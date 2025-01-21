import numpy as np
import random
import os
import re
import glob
from tqdm.contrib import tzip
from sklearn.neighbors import KDTree


class PointCloudProcessor:
    def __init__(self, threshold=150, use_normalized=True):
        """
        初始化处理类，设置默认的过滤阈值。
        
        参数：
        - threshold: int，粒子个数的默认阈值。
        """
        self.threshold = threshold
        self.use_normalized = use_normalized

    def _normalize(self, points):
        """
        使用横向动量 p_T 来归一化 px 和 py。
        
        参数：
        - points: np.ndarray，点云数据，shape: (num_points, num_dimensions)。
        
        返回：
        - normalized_points: np.ndarray，归一化后的点云数据，shape: (num_points, num_dimensions)。
        """
        px = points[:, 0]
        py = points[:, 1]

        # 计算横向动量 p_T
        p_T = np.sqrt(px**2 + py**2)

        # 避免除以 0，确保 p_T > 0
        p_T = np.maximum(p_T, 1e-6)  # 防止除以零

        # 使用 p_T 来归一化 px 和 py
        px_normalized = px / p_T
        py_normalized = py / p_T

        return np.vstack((px_normalized, py_normalized)).T

    def _rotate(self, points):
        """
        随机旋转点云，旋转角度在 [0, 2π] 范围内。
        
        参数：
        - points: np.ndarray，点云数据，shape: (num_points, num_dimensions)。
        
        返回：
        - rotated_points: np.ndarray，旋转后的点云数据，shape: (num_points, num_dimensions)。
        """
        theta = random.uniform(0, 2 * np.pi)  # 随机选择旋转角度

        # 旋转矩阵
        rotation_matrix = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), -np.cos(theta)]
        ])

        # 应用旋转矩阵到 px 和 py
        rotated_points = points[:, :2] @ rotation_matrix.T  # 只旋转前两维：px 和 py
        return np.hstack((rotated_points, np.zeros((rotated_points.shape[0], 1))))  # 将 pz 设置为 0

    def _get_sorted_files(self, file_pattern):
        """
        获取按规则排序的文件列表。
        
        参数：
        - file_pattern: str，文件路径的通配符模式（例如：'*-1.npy'）。
        
        返回：
        - sorted_files: list，按数字部分排序的文件列表。
        """
        files = glob.glob(file_pattern)
        sorted_files = sorted(files, key=lambda x: int(re.search(r'(\d+)-\d+\.npy', os.path.basename(x)).group(1)))
        return sorted_files

    def _filter_and_sample(self, data, indices):
        """
        删除零点并随机采样，使点数等于 threshold。
        
        参数：
        - data: np.ndarray，点云数据，shape: (num_point_clouds, num_points, num_dimensions)。
        - indices: np.ndarray，符合条件的点云索引。
        
        返回：
        - sampled_data: list，处理后的点云数据列表，每个元素是 shape (threshold, num_dimensions) 的数组。
        - valid_indices: list，成功处理的点云索引。
        """
        sampled_data = []
        valid_indices = []

        for i in indices:
            # 删除零点
            point_cloud = data[i, :, 0:3]
            non_zero_mask = np.any(point_cloud != 0, axis=1)
            filtered_points = point_cloud[non_zero_mask]

            # 如果剩余点数不足 threshold，则跳过该点云
            if filtered_points.shape[0] < self.threshold:
                print(f"点云索引 {i} 非零点数量不足 {self.threshold}，跳过。")
                continue

            # 随机采样
            sampled_indices = np.random.choice(filtered_points.shape[0], self.threshold, replace=False)
            sampled_points = filtered_points[sampled_indices]

            sampled_data.append(sampled_points)
            valid_indices.append(i)

        return sampled_data, valid_indices

    def process_and_save(self, data_pattern, label_pattern, output_dir, output_txt):
        """
        综合功能：筛选点云系统、删除零点、随机采样，并保存处理后的点云。
        
        参数：
        - data_pattern: str，点云数据文件路径的通配符模式（例如：'*-1.npy'）。
        - label_pattern: str，标签数据文件路径的通配符模式（例如：'*-3.npy'）。
        - output_dir: str，保存的目标目录。
        - output_txt: str，记录保存文件名的 TXT 文件路径。
        """
        # 获取按规则排序的数据文件和标签文件
        data_files = self._get_sorted_files(data_pattern)
        label_files = self._get_sorted_files(label_pattern)

        # 检查数据文件和标签文件数量是否一致
        if len(data_files) != len(label_files):
            raise ValueError("数据文件和标签文件的数量不一致，请检查文件命名和匹配规则！")

        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)

        # 获取当前目录中现有文件的编号，以确保编号唯一
        existing_files = [f for f in os.listdir(output_dir) if f.endswith('.npy')]
        existing_numbers = [
            int(re.search(r'(\d+)-(\d+)\.npy', f).group(2))
            for f in existing_files if re.search(r'(\d+)-(\d+)\.npy', f)
        ]
        next_number = max(existing_numbers, default=0) + 1

        # 打开输出 TXT 文件用于记录
        with open(output_txt, 'a') as txt_file:
            for data_file, label_file in tzip(data_files, label_files):
                #print(f"正在处理数据文件 {data_file} 和标签文件 {label_file}")
                
                # 加载点云和标签数据
                data = np.load(data_file)[:, 0, :, :]  # shape: (num_point_clouds, num_points, num_dimensions)
                labels = np.load(label_file)

                # 筛选符合条件的点云索引
                is_nonzero = np.any(data != 0, axis=-1)  # shape: (num_point_clouds, num_points)
                nonzero_counts = np.sum(is_nonzero, axis=1)
                indices = np.where(nonzero_counts > self.threshold)[0]

                # 删除零点并随机采样
                sampled_data, valid_indices = self._filter_and_sample(data, indices)


                # 保存处理后的点云和标签
                for point_cloud, idx in zip(sampled_data, valid_indices):
                    # normalize and rotate
                    if self.use_normalized:
                        point_cloud = self._normalize(point_cloud)
                        point_cloud = self._rotate(point_cloud)
                        point_cloud[:, 2] = 0 
                    label = np.argmax(labels[idx])
                    file_name = f"{label+1}-{next_number:07d}.npy"
                    file_path = os.path.join(output_dir, file_name)
                    np.save(file_path, point_cloud)
                    txt_file.write(file_name + '\n')
                    next_number += 1

        print(f"所有符合条件的点云已保存到 {output_dir}，文件名记录在 {output_txt}")


processor = PointCloudProcessor(threshold=128)
data_pattern = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data_processing/modeloutput_total/*-1.npy"
label_pattern = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data_processing/modeloutput_total/*-3.npy"
output_dir = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data/ampt_np_2D"
output_txt = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data/ampt_np_2D/txt/train.txt"

processor.process_and_save(data_pattern, label_pattern, output_dir, output_txt)


data_pattern = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data_processing/modeloutput_total/*-2.npy"
label_pattern = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data_processing/modeloutput_total/*-4.npy"
output_dir = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data/ampt_np_2D"
output_txt = "/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data/ampt_np_2D/txt/test.txt"

processor.process_and_save(data_pattern, label_pattern, output_dir, output_txt)