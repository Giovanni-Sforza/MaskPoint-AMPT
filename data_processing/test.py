import numpy as np

data_file = "/home/qiuyu/MaskPoint/data/ampt_np/0mb_Nch120_150_3D/1-0000023.npy"
label_file = "/home/qiuyu/MaskPoint/data_processing/model_output/NUM_POINTS_PbPb0mb_PbP0mb_Nch120_150_pT0.4_12.0_1events_4D_0.0001INITIAL_LR_0_2p4/train_points.npy"

data = np.load(data_file)  # data shape: (num_point_clouds, 1, num_points, num_dimensions)
labels = np.load(label_file)  # labels shape: (num_point_clouds, 1, one_hot_encode)
print(data[120,0:3])
print(labels.shape)

myarray = np.fromfile("data/ModelNet/modelnet40_normal_resampled/modelnet40_test_8192pts_fps.dat", dtype=float)
print("len(myarray)::", myarray.size)

def find_min_nonzero_points(data_file):
    """
    找到点云数据中每个模型包含的非零点的个数，并返回最小的非零点个数。
    
    参数：
    - data_file: str，点云数据的文件路径（NumPy 文件，shape 为 (num_point_clouds, 1, num_points, num_dimensions)）
    
    返回：
    - min_nonzero_points: int，点云数据中最小的非零点个数。
    """
    # 加载点云数据
    data = np.load(data_file)  # shape: (num_point_clouds, 1, num_points, num_dimensions)
    
    # 移除第二个维度
    data = data[:, 0, :, :]  # shape: (num_point_clouds, num_points, num_dimensions)
    
    # 计算每个点是否为零点（所有维度均为零）
    is_nonzero = np.any(data != 0, axis=-1)  # shape: (num_point_clouds, num_points)
    
    # 计算每个点云模型中非零点的个数
    nonzero_counts = np.sum(is_nonzero, axis=1)  # shape: (num_point_clouds,)
    
    # 找到最小的非零点个数
    #min_nonzero_points = np.min(nonzero_counts)
    
    return nonzero_counts

# 示例使用
#data_file = "input_data.npy"  # 输入点云数据文件路径
nonzero_counts = find_min_nonzero_points(label_file)
#nonzero_counts = np.average(nonzero_counts)
nonzero_counts = np.bincount(nonzero_counts)
print(f"最小的非零点个数是：{ nonzero_counts}")
#nonzero_counts = np.min(nonzero_counts)
print(f"最小的非零点个数是：{ nonzero_counts}")
