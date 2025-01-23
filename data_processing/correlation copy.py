import numpy as np
import os
import seaborn as sns
import pandas as pd
from tqdm.contrib import tzip
from itertools import combinations
import matplotlib.pyplot as plt
def calculate_observables(momentum_data, particle_mass=0.0):
    """
    计算单个事件的物理观测量，包括平均横动量 (pT), 赝快度 (eta), 椭圆流 (v2), 和三角流 (v3)
    参数:
        momentum_data: ndarray, 每行是一个粒子的动量 [p_x, p_y, p_z]
        particle_mass: float, 粒子的质量 (默认为 0)
    返回:
        observables: dict, 包含计算的物理观测量
    """
    px, py, pz = momentum_data[:, 0], momentum_data[:, 1], momentum_data[:, 2]
    pt = np.sqrt(px**2 + py**2)  # 横动量 pT
    phi = np.arctan2(py, px)     # 动量方位角 phi
    
    # 赝快度 eta
    # 粒子的纵向动量 pz 和横动量 pt
    
    theta = np.arctan2(pt, pz)  # 极角

    # 伪快度
    eta = -np.log(np.tan(theta / 2))
    
    """# 椭圆流 v2 和三角流 v3
    v2 = np.mean(np.cos(2 * phi))
    v3 = np.mean(np.cos(3 * phi))"""



    Qx2 = np.sum(pt * np.cos(2 * phi))
    Qy2 = np.sum(pt * np.sin(2 * phi))
    Psi_event2 = 0.5 * np.arctan2(Qy2, Qx2)  # 事件平面角
    delta_phi2 = phi - Psi_event2
    v2 = np.mean(np.cos(2 * delta_phi2))

    Qx3 = np.sum(pt * np.cos(3 * phi))
    Qy3 = np.sum(pt * np.sin(3 * phi))
    Psi_event3 = 1/3 * np.arctan2(Qy3, Qx3)  # 事件平面角
    delta_phi3 = phi - Psi_event3
    v3 = np.mean(np.cos(3 * delta_phi3))

    #print(f"v2: {Psi_event2}, v3: {Psi_event3},err:{abs(Psi_event2-Psi_event3)/Psi_event2}")

    # 计算结果
    observables = {
        "mean_pT": np.mean(pt),  # 平均横动量
        "mean_eta": np.mean(eta),  # 平均赝快度
        "v2": v2,
        "v3": v3
    }
    #return observables
    return np.array([observables["mean_pT"], observables["mean_eta"], observables["v2"], observables["v3"]])

# 加载降维后的特征数据
pca_data_path = f'/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data/features_pca.npz'
pca_data = np.load(pca_data_path, allow_pickle=True)
pca_features = pca_data['features']    # PCA 特征
labels = pca_data['labels']          # 标签
identifiers = pca_data['identifiers']  # 标识符

pca_features = pca_features[0:177714:80, :]
labels = labels[0:177714:80]
identifiers = identifiers[0:177714:80]

print(f"Loaded Labels: {labels.shape}, Identifiers: {identifiers.shape}")

# 定义原始文件路径
original_data_path = f'/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/data/ampt_np/total/'  # 假设原始文件存放在此目录
if not os.path.exists(original_data_path):
    raise FileNotFoundError(f"Original data directory does not exist: {original_data_path}")

# 存储每个事件的物理观测量
all_observables = []

# 遍历所有的标签和标识符，加载对应文件，计算物理量
for label, identifier in tzip(labels, identifiers):
    file_name = f"{label+1}-{identifier}.npy"  # 构造文件名
    file_path = os.path.join(original_data_path, file_name)
    
    if not os.path.exists(file_path):
        print(f"Warning: File not found {file_path}, skipping...")
        continue
    
    # 加载动量数据
    momentum_data = np.load(file_path)  # 每行是一个粒子的动量 [p_x, p_y, p_z]
    if momentum_data.shape[1] != 3:
        print(f"Invalid shape in file {file_path}: expected Nx3, got {momentum_data.shape}")
        continue
    
    # 计算物理观测量
    observables = calculate_observables(momentum_data)
    all_observables.append(observables)

# 将所有事件的观测量转换为结构化数组
observables_array = np.array(all_observables)
print("Shape of observables array:", observables_array.shape)

# 保存为 .npz 文件
observables_save_path = f"/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data/event_observables.npz"
np.savez(observables_save_path, observables=observables_array, labels=labels, identifiers=identifiers)
print(f"Event observables saved to: {observables_save_path}")


def plot_correlation_heatmap(pca_features, observables, pca_names, observable_names):

    #绘制 PCA 特征与物理观测量的相关性热力图

    if observables.ndim == 1:
        observables = observables.reshape(-1, 1)  # 将 1D 数组扩展为 2D 数组

    # 检查 pca_features 和 observables 的形状
    print(f"Shape of pca_features: {pca_features.shape}")
    print(f"Shape of observables: {observables.shape}")

    # 合并 PCA 特征和物理观测量到一个 DataFrame
    data = np.hstack([observables, pca_features])  # 横向拼接
    column_names = observable_names + pca_names
    df = pd.DataFrame(data, columns=column_names)

    # 计算相关性矩阵
    corr_matrix = df.corr()
    corr_matrix = abs(corr_matrix)  # 取绝对值
    plt.rcParams.update({
        "font.family": "serif",      # 使用衬线字体
        "font.size": 12,             # 全局字体大小
        "axes.titlesize": 14,        # 坐标轴标题字体大小
        "axes.labelsize": 12,        # 坐标轴标签字体大小
        "xtick.labelsize": 10,       # X轴刻度字体大小
        "ytick.labelsize": 10,       # Y轴刻度字体大小
        "legend.fontsize": 12,       # 图例字体大小
        "figure.dpi": 300            # 默认图像分辨率
    })
    # 绘制热力图
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix.loc[pca_names, observable_names], annot=True, cmap='mako', fmt=".2f", linewidths=0.5, cbar_kws={'shrink': .8}, square=True)
    plt.title("Correlation Heatmap: Observables vs PCA Features", fontsize=18, weight='bold')
    plt.xlabel("Observables", weight='bold')
    plt.ylabel("PCA Features",  weight='bold')
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=0)
    plt.tight_layout()
    output_path_png = f"/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data/correlation.png"
    plt.savefig(output_path_png, format='png', dpi=300, bbox_inches='tight')
    print(f"Cumulative variance plot saved to: {output_path_png}")


# 示例调用
pca_names = [f"PC_{i+1}" for i in range(pca_features.shape[1])]
observable_names = ["mean_pT", "mean_eta", "v2", "v3"]
plot_correlation_heatmap(pca_features, observables_array, pca_names, observable_names)
