import numpy as np
import os

def save_sorted_point_clouds(data_file, label_file, output_dir, output_txt):
    # open floder for save data
    os.makedirs(output_dir, exist_ok=True)
    
    # load pc data
    data = np.load(data_file)  # data shape: (num_point_clouds, 1, num_points, num_dimensions)
    labels = np.load(label_file)  # labels shape: (num_point_clouds, 1, one_hot_encode)
    
    # 检查数据格式
    if data.ndim != 4 or data.shape[1] != 1:
        raise ValueError("点云数据的维度应为 (num_point_clouds, 1, num_points, num_dimensions)")
    if labels.ndim != 2 :
        raise ValueError("标签数据的维度应为 (num_point_clouds, one_hot_encode)")
    if data.shape[0] != labels.shape[0]:
        raise ValueError("点云数据和标签数据的点云数量不一致")

    # 移除冗余维度
    data = data[:, 0, :, :]  # shape: (num_point_clouds, num_points, num_dimensions)
    labels = np.argmax(labels, axis=1)  # 转换为类别标签 (num_point_clouds,)
    
    # 按类别排序
    sorted_indices = np.argsort(labels)  # 获取按类别排序的索引
    data = data[sorted_indices]
    labels = labels[sorted_indices]
    
    # 创建一个字典存储每个类别的计数
    category_counts = {}
    
    # 打开 TXT 文件用于保存文件名
    with open(output_txt, 'a') as txt_file:
        for i, (point_cloud, label) in enumerate(zip(data, labels)):
            # 更新类别计数
            if label not in category_counts:
                category_counts[label] = 0
            category_counts[label] += 1
            
            # 生成紧凑的文件名，例如 1-0000001.npy
            file_name = f"{label+1}-{category_counts[label]:07d}.npy"
            file_path = os.path.join(output_dir, file_name)
            np.save(file_path, point_cloud)
            
            # 写入文件名到 TXT 文件
            txt_file.write(file_name + '\n')

    print(f"点云按类别保存完成，文件名记录在 {output_txt}")

# 示例使用
#data_file = "input_data.npy"  # 输入的点云数据文件路径
#label_file = "input_labels.npy"  # 输入的标签数据文件路径
#output_dir = "sorted_point_clouds"  # 输出文件夹
#output_txt = "sorted_point_cloud_files.txt"  # 输出的 TXT 文件

data_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch90_120_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_points.npy"
label_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch90_120_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch90_120"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/train_0mb_Nch90_120.txt"
save_sorted_point_clouds(data_file, label_file, output_dir, output_txt)

data_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch90_120_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_points.npy"
label_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch90_120_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch90_120"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/validation_0mb_Nch90_120.txt"
save_sorted_point_clouds(data_file2, label_file2, output_dir, output_txt)

