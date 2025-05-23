import torch
import torch.nn as nn
from tools import builder
from utils import misc, dist_utils
import time
from utils.logger import *
from utils.AverageMeter import AverageMeter
from sklearn.svm import LinearSVC,SVC
import numpy as np
from torchvision import transforms
from datasets import data_transforms
from utils.grad_utils import IterativePercentile
from sklearn.linear_model import LogisticRegression
train_transforms = transforms.Compose(
    [
        data_transforms.PointcloudScaleAndTranslate(),
    ]
)

class Acc_Metric:
    def __init__(self, acc = 0.):
        if type(acc).__name__ == 'dict':
            self.acc = acc['acc']
        else:
            self.acc = acc

    def better_than(self, other):
        if self.acc > other.acc:
            return True
        else:
            return False

    def state_dict(self):
        _dict = dict()
        _dict['acc'] = self.acc
        return _dict




def evaluate_svm(train_features, train_labels, test_features, test_labels):
    #clf = LinearSVC()
    #clf = SVC(kernel='linear', C=1)
    #clf.fit(train_features, train_labels)
    #pred = clf.predict(test_features)
    # 初始化 Logistic Regression 模型
    log_reg = LogisticRegression(
        penalty='l2',            # 使用 L2 正则化
        solver='saga',           # saga 是适合大规模数据的优化器
        max_iter=500,            # 最大迭代次数
        C=1.0,                   # 正则化强度（较小的 C 值加强正则化）
        random_state=42          # 固定随机种子
    )
    # 训练模型
    log_reg.fit(train_features, train_labels)
    # 预测测试集
    pred = log_reg.predict(test_features)
    return np.sum(test_labels == pred) * 1. / pred.shape[0]

def run_net(args, config, train_writer=None, val_writer=None):
    logger = get_logger(args.log_name)
    # build dataset
    (train_sampler, train_dataloader), (_, test_dataloader),= builder.dataset_builder(args, config.dataset.train), \
                                                            builder.dataset_builder(args, config.dataset.val)
    (_, extra_train_dataloader)  = builder.dataset_builder(args, config.dataset.extra_train) if config.dataset.get('extra_train') else (None, None)
    # build model
    base_model = builder.model_builder(config.model)
    if args.use_gpu:
        base_model.to(args.local_rank)

    # from IPython import embed; embed()
    
    # parameter setting
    start_epoch = 0
    best_metrics = Acc_Metric(0.)
    metrics = Acc_Metric(0.)

    # resume ckpts
    if args.resume:
        start_epoch, best_metric = builder.resume_model(base_model, args, logger = logger)
        best_metrics = Acc_Metric(best_metric)
    elif args.start_ckpts is not None:
        builder.load_model(base_model, args.start_ckpts, logger = logger)

    # DDP
    if args.distributed:
        # Sync BN
        if args.sync_bn:
            base_model = torch.nn.SyncBatchNorm.convert_sync_batchnorm(base_model)
            print_log('Using Synchronized BatchNorm ...', logger = logger)
        base_model = nn.parallel.DistributedDataParallel(base_model, device_ids=[args.local_rank % torch.cuda.device_count()], find_unused_parameters=True)
        print_log('Using Distributed Data parallel ...' , logger = logger)
    else:
        print_log('Using Data parallel ...' , logger = logger)
        base_model = nn.DataParallel(base_model).cuda()
    # optimizer & scheduler
    optimizer, scheduler = builder.build_opti_sche(base_model, config)
    
    if args.resume:
        builder.resume_optimizer(optimizer, args, logger = logger)

    spike_idx = 0
    gradclip_percentile = config.get('gradclip_percentile', -1)
    if gradclip_percentile > 0:
        grad_history = IterativePercentile(p=gradclip_percentile)
    else:
        grad_history = None

    # trainval
    # training
    base_model.zero_grad()
    fixed_sample = None
    for epoch in range(start_epoch, config.max_epoch + 1):
        if args.distributed:
            train_sampler.set_epoch(epoch)
        #torch.cuda.empty_cache()
        base_model.train()

        epoch_start_time = time.time()
        batch_start_time = time.time()
        batch_time = AverageMeter()
        data_time = AverageMeter()
        losses = AverageMeter(['Loss1', 'Loss2'])
        acc_avg = AverageMeter()
        num_iter = 0
        grad_clip_val = config.grad_norm_clip

        base_model.train()  # set model to training mode
        n_batches = len(train_dataloader)
        for idx, (taxonomy_ids, model_ids, data) in enumerate(train_dataloader):
            num_iter += 1
            n_itr = epoch * n_batches + idx
            
            data_time.update(time.time() - batch_start_time)
            npoints = config.dataset.train.others.npoints
            dataset_name = config.dataset.train._base_.NAME
            if dataset_name == 'ShapeNet':
                points = data.cuda()
            elif dataset_name == 'ModelNet':
                points = data[0].cuda()
                points = misc.fps(points, npoints)
            elif dataset_name.startswith('ScanNet'):
                points = data.cuda()
            elif dataset_name== 'Amptdata':
                points = data.cuda()
            else:
                raise NotImplementedError(f'Train phase do not support {dataset_name}')

            assert points.size(1) == npoints
            points = train_transforms(points)

            if args.overfit_single_batch:
                if fixed_sample is None:
                    fixed_sample = points.clone()
                else:
                    points = fixed_sample.clone()

            loss_1, loss_2, acc = base_model(points)

            _loss = loss_1 + loss_2

            _loss.backward()

            # forward
            if num_iter == config.step_per_update:
                if config.get('grad_norm_clip') is not None:
                    if grad_history is not None:
                        grad_norm = torch.cat([p.grad.view(-1) for p in base_model.parameters() if p.grad is not None]).norm().item()
                        grad_clip_val = grad_history.add(grad_norm)
                        if train_writer is not None:
                            train_writer.add_scalar('Loss/Batch/clip_val', grad_clip_val, n_itr)
                    grad_norm = torch.nn.utils.clip_grad_norm_(base_model.parameters(), grad_clip_val, norm_type=2)
                    if train_writer is not None:
                        train_writer.add_scalar('Loss/grad_norm', grad_norm.item(), n_itr)
                num_iter = 0
                optimizer.step()
                base_model.zero_grad()

            if args.distributed:
                loss_1 = dist_utils.reduce_tensor(loss_1, args)
                loss_2 = dist_utils.reduce_tensor(loss_2, args)
                losses.update([loss_1.item(), loss_2.item()])
            else:
                losses.update([loss_1.item(), loss_2.item()])
                acc_avg.update(acc.item())

            if args.distributed:
                torch.cuda.synchronize()
                acc_avg.update(acc.item())


            if train_writer is not None:
                train_writer.add_scalar('Loss/Batch/Loss_1', loss_1.item(), n_itr)
                train_writer.add_scalar('Loss/Batch/Loss_2', loss_2.item(), n_itr)
                train_writer.add_scalar('Loss/Batch/Accuracy', acc.item(), n_itr)
                train_writer.add_scalar('Loss/Batch/LR', optimizer.param_groups[0]['lr'], n_itr)


            batch_time.update(time.time() - batch_start_time)
            batch_start_time = time.time()

            if idx % 20 == 0:
                print_log('[Epoch %d/%d][Batch %d/%d] BatchTime = %.3f (s) DataTime = %.3f (s) Losses = %s lr = %.6f acc = %s ' %
                            (epoch, config.max_epoch, idx + 1, n_batches, batch_time.val(), data_time.val(),
                            ['%.4f' % l for l in losses.val()], optimizer.param_groups[0]['lr'],acc_avg.val()), logger = logger)
        if isinstance(scheduler, list):
            for item in scheduler:
                item.step(epoch)
        else:
            scheduler.step(epoch)
        epoch_end_time = time.time()

        if train_writer is not None:
            train_writer.add_scalar('Loss/Epoch/Loss_1', losses.avg(0), epoch)
            train_writer.add_scalar('Loss/Epoch/Loss_2', losses.avg(1), epoch)
            train_writer.add_scalar('Loss/Epoch/Accuracy', acc_avg.avg(), epoch)


        print_log('[Training] EPOCH: %d EpochTime = %.3f (s) Losses = %s Acc = %s ' %
            (epoch,  epoch_end_time - epoch_start_time, ['%.4f' % l for l in losses.avg()],acc_avg.avg()), logger = logger)

        if epoch % args.val_freq == 0 and epoch != 0:
            # Validate the current model
            metrics = validate(base_model, extra_train_dataloader, test_dataloader, epoch, val_writer, args, config, logger=logger)

            # Save ckeckpoints
            if metrics.better_than(best_metrics):
                best_metrics = metrics
                builder.save_checkpoint(base_model, optimizer, epoch, metrics, best_metrics, 'ckpt-best', args, logger = logger)
        builder.save_checkpoint(base_model, optimizer, epoch, metrics, best_metrics, 'ckpt-last', args, logger = logger)
        if (config.max_epoch - epoch) < 10:
            builder.save_checkpoint(base_model, optimizer, epoch, metrics, best_metrics, f'ckpt-epoch-{epoch:03d}', args, logger = logger)     
    if train_writer is not None:
        train_writer.close()
    if val_writer is not None:
        val_writer.close()

#initial one
"""def validate(base_model, extra_train_dataloader, test_dataloader, epoch, val_writer, args, config, logger = None):
    print_log(f"[VALIDATION] Start validating epoch {epoch}", logger = logger)
    base_model.eval()  # set model to eval mode

    test_features = []
    test_label = []

    train_features = []
    train_label = []
    npoints = config.dataset.train.others.npoints
    with torch.no_grad():
        for idx, (taxonomy_ids, model_ids, data) in enumerate(extra_train_dataloader):
            points = data[0].cuda()
            label = data[1].cuda()

            points = misc.fps(points, npoints)

            assert points.size(1) == npoints
            feature = base_model(points, noaug=True)
            target = label.view(-1)

            train_features.append(feature.detach())
            train_label.append(target.detach())

        for idx, (taxonomy_ids, model_ids, data) in enumerate(test_dataloader):
            points = data[0].cuda()
            label = data[1].cuda()

            points = misc.fps(points, npoints)
            assert points.size(1) == npoints
            feature = base_model(points, noaug=True)
            target = label.view(-1)

            test_features.append(feature.detach())
            test_label.append(target.detach())


        train_features = torch.cat(train_features, dim=0)
        train_label = torch.cat(train_label, dim=0)
        test_features = torch.cat(test_features, dim=0)
        test_label = torch.cat(test_label, dim=0)

        if args.distributed:
            train_features = dist_utils.gather_tensor(train_features, args)
            train_label = dist_utils.gather_tensor(train_label, args)
            test_features = dist_utils.gather_tensor(test_features, args)
            test_label = dist_utils.gather_tensor(test_label, args)

        svm_acc = evaluate_svm(train_features.data.cpu().numpy(), train_label.data.cpu().numpy(), test_features.data.cpu().numpy(), test_label.data.cpu().numpy())

        print_log('[Validation] EPOCH: %d  acc = %.4f' % (epoch,svm_acc), logger=logger)

        if args.distributed:
            torch.cuda.synchronize()

    # Add testing results to TensorBoard
    if val_writer is not None:
        val_writer.add_scalar('Metric/ACC', svm_acc, epoch)

    return Acc_Metric(svm_acc)"""

def validate(base_model, extra_train_dataloader, test_dataloader, epoch, val_writer, args, config, logger=None):
    print_log(f"[VALIDATION] Start validating epoch {epoch}", logger=logger)
    base_model.eval()  # set model to eval mode

    test_features, test_labels = [], []
    train_features, train_labels = [], []
    npoints = config.dataset.train.others.npoints

    with torch.no_grad():
        for idx, (taxonomy_ids, model_ids, data) in enumerate(extra_train_dataloader):
            points = data.cuda()  # Directly use the points as ShapeNet format
            labels = torch.tensor([int(taxonomy_id) for taxonomy_id in taxonomy_ids]).cuda()  # Convert taxonomy_ids to tensor and move to GPU

            points = misc.fps(points, npoints)
            assert points.size(1) == npoints

            features = base_model(points, noaug=True)
            train_features.append(features.detach())
            train_labels.append(labels)  # Add the label tensor

        for idx, (taxonomy_ids, model_ids, data) in enumerate(test_dataloader):
            points = data.cuda()
            labels = torch.tensor([int(taxonomy_id) for taxonomy_id in taxonomy_ids]).cuda()  # Convert taxonomy_ids to tensor and move to GPU

            points = misc.fps(points, npoints)
            assert points.size(1) == npoints

            features = base_model(points, noaug=True)
            test_features.append(features.detach())
            test_labels.append(labels)

        train_features = torch.cat(train_features, dim=0)
        train_labels = torch.cat(train_labels, dim=0)  # Concatenate the label tensors
        test_features = torch.cat(test_features, dim=0)
        test_labels = torch.cat(test_labels, dim=0)

        if args.distributed:
            train_features = dist_utils.gather_tensor(train_features, args)
            train_labels = dist_utils.gather_tensor(train_labels, args)
            test_features = dist_utils.gather_tensor(test_features, args)
            test_labels = dist_utils.gather_tensor(test_labels, args)

        svm_acc = evaluate_svm(
            train_features.cpu().numpy(),
            train_labels.cpu().numpy(),
            test_features.cpu().numpy(),
            test_labels.cpu().numpy(),
        )

        print_log('[Validation] EPOCH: %d  acc = %.4f' % (epoch, svm_acc), logger=logger)

        if args.distributed:
            torch.cuda.synchronize()

    if val_writer is not None:
        val_writer.add_scalar('Metric/ACC', svm_acc, epoch)

    return Acc_Metric(svm_acc)


def test_net():
    pass
