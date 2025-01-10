import os

import torch
import torch.multiprocessing as mp
from torch import distributed as dist



def init_dist(launcher, backend='nccl', **kwargs):
    if mp.get_start_method(allow_none=True) is None:
        mp.set_start_method('spawn')
    if launcher == 'pytorch':
        _init_dist_pytorch(backend, **kwargs)
    else:
        raise ValueError(f'Invalid launcher type: {launcher}')


def _init_dist_pytorch(backend, **kwargs):
    # TODO: use local_rank instead of rank % num_gpus
    rank = int(os.environ['RANK'])
    num_gpus = torch.cuda.device_count()
    torch.cuda.set_device(rank % num_gpus)
    dist.init_process_group(backend=backend, **kwargs)
    print(f'init distributed in rank {torch.distributed.get_rank()}')


def get_dist_info():
    if dist.is_available():
        initialized = dist.is_initialized()
    else:
        initialized = False
    if initialized:
        rank = dist.get_rank()
        world_size = dist.get_world_size()
    else:
        rank = 0
        world_size = 1
    return rank, world_size


def reduce_tensor(tensor, args):
    '''
        for acc kind, get the mean in each gpu
    '''
    rt = tensor.clone()
    torch.distributed.all_reduce(rt, op=torch.distributed.ReduceOp.SUM)
    rt /= args.world_size
    return rt

def gather_tensor(tensor, args):
    output_tensors = [tensor.clone() for _ in range(args.world_size)]
    torch.distributed.all_gather(output_tensors, tensor)
    concat = torch.cat(output_tensors, dim=0)
    return concat


"""The **`/home/qiuyu/MaskPoint/utils/dist_utils.py`** file is designed to facilitate distributed computing using PyTorch's distributed package. Here's a breakdown of its purpose and functionality:

1. **Initialization of Distributed Environment:**
   - The `init_dist` function initializes the distributed environment based on the specified launcher type. It currently supports the 'pytorch' launcher, which uses PyTorch's native distributed package.
   - The `_init_dist_pytorch` function sets up the distributed process group using the specified backend (e.g., 'nccl' for NVIDIA Collective Communications Library). It determines the rank of the current process and sets the appropriate GPU device for the process.

2. **Distributed Information Retrieval:**
   - The `get_dist_info` function provides information about the current distributed environment, such as the rank of the current process and the total number of processes (world size). This is useful for coordinating tasks across multiple processes.

3. **Tensor Reduction and Gathering:**
   - The `reduce_tensor` function performs a reduction operation (specifically, a sum) across all processes for a given tensor. It then averages the result by dividing by the world size. This is commonly used to aggregate results like gradients or metrics across multiple GPUs.
   - The `gather_tensor` function collects tensors from all processes and concatenates them along a specified dimension. This is useful for gathering results from different processes into a single tensor for further processing or analysis.

Overall, this file provides utility functions to manage and operate in a distributed computing environment, enabling parallel processing and efficient use of multiple GPUs for tasks such as training machine learning models."""