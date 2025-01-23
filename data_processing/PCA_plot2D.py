import argparse
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from tqdm.contrib import tzip
# 参数解析
parser = argparse.ArgumentParser(description="PCA Analysis for Extracted Features")
parser.add_argument("--experiment_path", type=str, default="/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data/",
                    help="Path to save or load extracted features")
args = parser.parse_args()

# 确保路径存在
if not os.path.exists(args.experiment_path):
    os.makedirs(args.experiment_path)

# 加载 .npz 文件
data_path = os.path.join(args.experiment_path, "features_pca2D.npz")
data = np.load(data_path, allow_pickle=True)


# 提取特征和标签
features = data['features']  # 特征矩阵 (num_samples, classify(1)+group_num, feature_dim)

labels = data['labels']      # 标签数组 (num_samples,)

features = features[0:177714:70, :]
labels = labels[0:177714:70]
# 配色方案和标记
colors = ["#e64532", "#3c7fb1"] 
markers = ["o", "s"]  
plt.style.use('seaborn-paper')  # 或 'ggplot', 'seaborn-whitegrid'
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
for data,lab in tzip(features,labels):
    plt.scatter(data[0],data[1],color=colors[lab],marker=markers[lab],s=10,alpha=0.85, edgecolors='w', linewidth=0.5)

plt.scatter([], [], color=colors[0], marker=markers[0], label="Pbp")
plt.scatter([], [], color=colors[1], marker=markers[1], label="PbPb")
plt.legend(loc='best', frameon=True, shadow=True, fancybox=True)

# 添加标签和标题
plt.xlabel("Component 1")
plt.ylabel("Component 2")
plt.title("PCA 2D Plot")
plt.savefig("/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data/PCA2D.png",dpi=300)