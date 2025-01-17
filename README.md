# MaskPoint-AMPT

MaskPoint-AMPT is a project designed for point cloud data processing, based on the framework described in the paper **"Masked Discrimination for Self-Supervised Learning on Point Clouds"** by Liu, Haotian et al. (ECCV 2022). This project integrates deep learning techniques with high-energy physics research, focusing on analyzing AMPT (A Multi-Phase Transport) simulation data.

## Objectives

The primary goals of the MaskPoint-AMPT project are:

1. **Final State Particle Reconstruction**: Employ self-supervised learning to reconstruct final-state particles in high-energy heavy-ion collision data.
2. **Feature Extraction and Correlation Analysis**: Extract high-dimensional features from AMPT data that reflect underlying physical phenomena, and investigate their correlations with traditional physical observables, such as collision system types, average transverse momentum, elliptic flow, and triangular flow.
3. **Interpretability Analysis**: Enhance the interpretability of deep learning models in high-energy physics research by analyzing the features extracted by the autoencoder, shedding light on potential patterns in small and large collision systems (e.g., Pb-P vs. Pb-Pb systems).

## Project Structure

```
├── cfgs/                 # Configuration files for experiments and datasets
├── data/                 # Symlinked folder for storing input data
├── data_processing/      # Scripts for preprocessing AMPT data
├── datasets/             # Dataset loaders and transformation utilities
├── experiments/          # Logs and results for experiments
├── extensions/           # PointNet++ and custom point cloud operators
├── figure/               # Visual resources such as diagrams
├── models/               # Core models, including MaskPoint
├── scripts/              # Helper scripts for testing
├── tools/                # Training and pretraining utilities
├── utils/                # Helper utilities for metrics, logging, and configuration
├── install.sh            # Installation script for environment setup
├── main.py               # Entry point for training and evaluation
├── requirements.txt      # Python dependencies
├── LICENSE               # License information
└── README.md             # Project documentation (this file)
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Giovanni-Sforza/MaskPoint-AMPT.git
   cd MaskPoint-AMPT
   ```

2. Install dependencies:
   Ensure your system meets the following requirements:
   - PyTorch >= 1.7.0
   - Python >= 3.7
   - CUDA >= 9.0
   - GCC >= 4.9
   - torchvision

   if your meet the trouble that `version `GLIBC_2.18' not found`, you can try to install the open3d == 0.9.0 in env python == 3.7:
   ```bash
   conda search -c open3d-admin open3d==0.9.0
   ```

   Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

   if you find that pip have not install these pkgs in the envs you need:
   ```bash
   vim ~/anaconda3/envs/MP10/lib/python3.7/site.py
   ```

   and change `USER_SITE`, `USER_BASE`:
   ```bash
   USER_SITE = /storage/fdunphome/zhangjingzong/anaconda3/envs/MP10/lib/python3.7/site-packages
   USER_BASE = /storage/fdunphome/zhangjingzong/anaconda3/envs/MP10
   ```

3. Run the installation script:
   ```bash
   bash install.sh
   ```
please make sure that gcc has been updated. you can use conda to install some new version of gcc, and add the PATH in ~/.bashrc

---
**Important Notes**:
   - If your system already includes GCC version >= 4.9, running `bash install.sh` should work fine.
   - If GCC version is lower than 4.9, you will need to install GCC 9.4.0 and G++ 9.4.0 in the appropriate virtual environment using:
     ```bash
     conda install -c conda-forge gcc=9.4.0 g++=9.4.0
     ```
     Then, set the environment variables as follows:
     ```bash
     export PATH="/storage/fdunphome/zhangjingzong/anaconda3/envs/base_gcc_9_4_0/bin:$PATH"
     export CC=/storage/fdunphome/zhangjingzong/anaconda3/envs/base_gcc_9_4_0/bin/gcc
     export CXX=/storage/fdunphome/zhangjingzong/anaconda3/envs/base_gcc_9_4_0/bin/g++
     export LD_LIBRARY_PATH=/storage/fdunphome/zhangjingzong/anaconda3/envs/base_gcc_9_4_0/lib:$LD_LIBRARY_PATH
     ```
     If you have multiple environments with GCC and G++, ensure the environment variables point to the folder of the virtual environment where the program will run.


## Usage

### Pretraining
To start pretraining with the AMPT dataset, run the following command:
```bash
python main.py --config cfgs/pretrain_ampt_data.yaml \
    --exp_name ampt \
    --val_freq 10
```

- `--config`: Path to the configuration file
- `--exp_name`: Name of the experiment
- `--val_freq`: Validation frequency during training

### Fineturning

Coming soon !

## License
This project is licensed under the terms of the [LICENSE](LICENSE) file.

## Acknowledgments
Special thanks to the developers of PointNet++ and the authors of related datasets.

