from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# 假设 data 是 N x 96 的特征矩阵
data = np.random.rand(100, 96)  # 示例数据
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data)

# 执行 PCA
pca = PCA()
pca.fit(data_scaled)

# 绘制累计方差图
plt.figure(figsize=(8, 5))
plt.plot(np.cumsum(pca.explained_variance_ratio_), marker='o')
plt.xlabel('Number of Components')
plt.ylabel('Cumulative Explained Variance')
plt.title('Cumulative Explained Variance by PCA')
plt.grid(True)
plt.show()
