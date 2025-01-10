import numpy as np
import os
import re

def save_sorted_point_clouds(data_file, label_file, output_dir, output_txt):
    # 创建保存的文件夹
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取目录中所有现有的 .npy 文件
    existing_files = [f for f in os.listdir(output_dir) if f.endswith('.npy')]
    
    # 提取现有文件名中的编号部分
    existing_numbers = []
    for filename in existing_files:
        match = re.match(r"(\d+)-(\d+)\.npy", filename)
        if match:
            num = int(match.group(2))  # 获取文件名中的编号部分
            existing_numbers.append(num)
    
    # 获取下一个可用的编号
    next_number = max(existing_numbers, default=0) + 1  # 如果没有现有文件，默认从 1 开始
    
    # 加载点云数据和分类数据
    data = np.load(data_file)  # data shape: (num_point_clouds, 1, num_points, num_dimensions)
    labels = np.load(label_file)  # labels shape: (num_point_clouds, 2)
    
    # 检查数据格式
    if data.ndim != 4 or data.shape[1] != 1:
        raise ValueError("点云数据的维度应为 (num_point_clouds, 1, num_points, num_dimensions)")
    if labels.ndim != 2 or labels.shape[1] != 2:
        raise ValueError("标签数据的维度应为 (num_point_clouds, 2)")
    if data.shape[0] != labels.shape[0]:
        raise ValueError("点云数据和标签数据的点云数量不一致")
    
    # 移除冗余维度
    #data = data[:, 0, :, :]  # shape: (num_point_clouds, num_points, num_dimensions)
    data = data[:, 0, :,0:3]
    # 将 one-hot 编码转换为类别索引
    labels = np.argmax(labels, axis=1)  # shape: (num_point_clouds,)
    
    # 按类别排序
    sorted_indices = np.argsort(labels)  # 获取按类别排序的索引
    data = data[sorted_indices]
    labels = labels[sorted_indices]
    
    # 创建一个字典存储每个类别的计数
    category_counts = {}
    
    # 打开 TXT 文件用于保存文件名，使用追加模式 'a'
    with open(output_txt, 'a') as txt_file:
        for i, (point_cloud, label) in enumerate(zip(data, labels)):
            # 更新类别计数
            if label not in category_counts:
                category_counts[label] = 0
            category_counts[label] += 1
            
            # 生成文件名，例如 1-0000001.npy
            file_name = f"{label+1}-{next_number:07d}.npy"
            file_path = os.path.join(output_dir, file_name)
            np.save(file_path, point_cloud)
            
            # 写入文件名到 TXT 文件（追加模式不会覆盖现有内容）
            txt_file.write(file_name + '\n')
            
            # 更新下一个编号
            next_number += 1

    print(f"点云按类别保存完成，文件名记录在 {output_txt}")

# 示例使用
# 示例使用
"""data_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch120_150_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_points.npy"
label_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch120_150_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch120_150_3D"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/0mb_Nch120_150_3D/train.txt"
save_sorted_point_clouds(data_file, label_file, output_dir, output_txt)

data_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch120_150_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_points.npy"
label_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch120_150_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch120_150_3D"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/0mb_Nch120_150_3D/test.txt"
save_sorted_point_clouds(data_file2, label_file2, output_dir, output_txt)
"""

# 示例使用
data_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch150_185_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_points.npy"
label_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch150_185_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch150_185_3D"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/0mb_Nch150_185_3D/train.txt"
save_sorted_point_clouds(data_file, label_file, output_dir, output_txt)

data_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch150_185_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_points.npy"
label_file2 = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch150_185_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/validation_labels.npy"
output_dir = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch150_185_3D"
output_txt = "/home/qiuyu/MaskPoint/data/ampt_np/txt/0mb_Nch150_185_3D/test.txt"
save_sorted_point_clouds(data_file2, label_file2, output_dir, output_txt)
