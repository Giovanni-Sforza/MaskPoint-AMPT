import argparse
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# 参数解析
parser = argparse.ArgumentParser(description="PCA Analysis for Extracted Features")
parser.add_argument("--experiment_path", type=str, default="/storage/fdunphome/zhangjingzong/MaskPoint-AMPT/experiments/finetune_ampt_data/test_fineturn_ampt_data",
                    help="Path to save or load extracted features")
args = parser.parse_args()

# 确保路径存在
if not os.path.exists(args.experiment_path):
    os.makedirs(args.experiment_path)

# 加载 .npz 文件
data_path = os.path.join(args.experiment_path, "features.npz")
data = np.load(data_path, allow_pickle=True)


# 提取特征和标签
features = data['features']  # 特征矩阵 (num_samples, classify(1)+group_num, feature_dim)
features = features[:, 0, :]

labels = data['labels']      # 标签数组 (num_samples,)
identifiers = data['identifiers']  # 标识符 (num_samples,)
print("Shape of features:", features.shape)
#print(features[0:10,1,1])
# 标准化特征
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features)

# 执行 PCA
pca = PCA()
pca.fit(features_scaled)

def plot_PCA(_pca):
    # 设置全局样式
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

    # 创建累计方差图
    plt.figure(figsize=(6, 4))  # 图像大小 (宽, 高)
    cumulative_variance = np.cumsum(_pca.explained_variance_ratio_)
    plt.plot(cumulative_variance, marker='o', linestyle='-', linewidth=1.0, color='#3c7fb1', label='Cumulative Explained Variance')

    # 添加两条基准线
    plt.axhline(y=0.95, color='#e64532', linestyle='--', linewidth=1, label='95% Variance Threshold')
    plt.axhline(y=0.99, color='#b53289', linestyle='--', linewidth=1, label='99% Variance Threshold')

    # 添加关键点标记
    optimal_components_95 = np.argmax(cumulative_variance >= 0.95)   # 找到达到95%方差的主成分数量
    optimal_components_99 = np.argmax(cumulative_variance >= 0.99)   # 找到达到99%方差的主成分数量
    #plt.axvline(x=optimal_components_95 - 1, color='g', linestyle='--', linewidth=1, label=f'95%: {optimal_components_95} Components')
    #plt.axvline(x=optimal_components_99 - 1, color='orange', linestyle='--', linewidth=1, label=f'99%: {optimal_components_99} Components')

    # 在曲线上标注关键点
    plt.scatter(optimal_components_95, cumulative_variance[optimal_components_95], color='#e64532', zorder=5)
    plt.text(optimal_components_95, cumulative_variance[optimal_components_95], f' Component {optimal_components_95+1}', color='#e64532')

    plt.scatter(optimal_components_99, cumulative_variance[optimal_components_99], color='#b53289', zorder=5)
    plt.text(optimal_components_99, cumulative_variance[optimal_components_99], f' Component {optimal_components_99+1}', color='#b53289')

    # 添加标签和标题
    plt.title("Cumulative Explained Variance by PCA")
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    #plt.xticks(range(1, 20))  # 自定义 X 轴刻度
    plt.yticks(np.arange(0.8, 1.01, 0.02))  # 自定义 Y 轴刻度
    plt.legend(loc='lower right')  # 图例位置
    plt.grid(alpha=0.5)  # 网格透明度

    """# 保存图像为高分辨率 PDF 或 PNG
    output_path = f"{args.experiment_path}/pca_cumulative_variance.pdf"
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    print(f"Cumulative variance plot saved to: {output_path}")"""

    # 如果需要导出 PNG 格式
    output_path_png = f"{args.experiment_path}/pca_cumulative_variance.png"
    plt.savefig(output_path_png, format='png', dpi=300, bbox_inches='tight')
    print(f"Cumulative variance plot saved to: {output_path_png}") 
plot_PCA(pca)
# 打印前几个主成分的重要性
explained_variance = pca.explained_variance_ratio_
print("Explained variance by each component:")
for i, var in enumerate(explained_variance[:15]):  # 只打印前10个
    print(f"Component {i+1}: {var:.4f}")

# 设置降维目标为95%的累计方差
n_components = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.95) + 1
print(f"Number of components to explain 95% variance: {n_components}")

# 设置降维目标为95%的累计方差
n_components = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.99) + 1
print(f"Number of components to explain 99% variance: {n_components}")

# 降维
pca = PCA(n_components=2)
features_reduced = pca.fit_transform(features_scaled)

"""# 保存降维后的特征
np.savez(
    f'{args.experiment_path}/features_pca',
    features=features_reduced,  # 降维后的特征
    labels=labels,  # 原始标签
    identifiers=identifiers  # 原始标识符
)"""
print(f"Reduced features saved to {args.experiment_path}/features_pca.npz")
