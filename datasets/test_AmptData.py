import os
import yaml
import torch
from your_module import Amptdata  # 替换成实际的模块路径

# 示例配置类
class Config:
    def __init__(self):
        self.DATA_PATH = "data/ampt_np/0mb_Nch90_120_3D"  # 替换为实际路径
        self.PC_PATH = "/path/to/pc_data"     # 替换为实际路径
        self.subset = "train"                 # 替换为实际子集名，例如 train/test/val
        self.N_POINTS = 2048                  # 替换为实际点数
        self.whole = False                    # 替换为实际配置

# 测试代码
def test_dataset():
    # 加载配置
    config = Config()
    
    # 创建数据集实例
    dataset = Amptdata(config)
    
    # 测试数据集长度
    print(f"Dataset size: {len(dataset)}")
    
    # 随机抽样测试一个数据
    idx = 0  # 你可以改成随机数，例如 idx = torch.randint(0, len(dataset), (1,)).item()
    taxonomy_id, model_id, data = dataset[idx]
    
    # 打印测试结果
    print("Sample output:")
    print(f"  Taxonomy ID: {taxonomy_id}")
    print(f"  Model ID: {model_id}")
    print(f"  Data shape: {data.shape}")
    print(f"  Data sample (first 5 points):\n{data[:5]}")

# 运行测试
if __name__ == "__main__":
    test_dataset()
