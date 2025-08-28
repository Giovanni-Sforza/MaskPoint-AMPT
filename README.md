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
├── figure/               # Visual resources
├── models/               # Core models, including MaskPoint
├── scripts/              # Helper scripts for testing
├── tools/                # Training and pretraining utilities
├── utils/                # Helper utilities for metrics, logging, and configuration
├── install.sh            # Installation script for environment setup
├── main.py               # Entry point for training and evaluation
├── requirements.txt      # Python dependencies
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

   if your meet the trouble that `version `GLIBC_2.18' not found`, which may be caused by open3d. so you can try to install the open3d == 0.9.0 in env python == 3.7:
   ```bash
   conda search -c open3d-admin open3d==0.9.0
   ```

   then, install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

   if you find that pip have not install these pkgs in the envs you need:
   ```bash
   vim ~/anaconda3/envs/MP10/lib/python3.7/site.py
   ```

   and change `USER_SITE`, `USER_BASE`:
   ```bash
   USER_SITE = ～/anaconda3/envs/MP10/lib/python3.7/site-packages
   USER_BASE = ～/anaconda3/envs/MP10
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
     export PATH="path_to_anaconda3_envs/bin:$PATH"
     export CC="path_to_anaconda3_envs/bin/gcc"
     export CXX="path_to_anaconda3_envs/bin/g++"
     export LD_LIBRARY_PATH="path_to_anaconda3_envs/lib:$LD_LIBRARY_PATH"
     ```
     If you have multiple environments with GCC and G++, or you are using gcc in a virtual env, as shown above, ensure the environment variables point to the folder of the virtual environment where the program will run, since pkgs in extensions will use gcc and be installed in the same env as the gcc you are using. 


## Usage

### Pretraining
To start pretraining with the AMPT dataset, run the following command:
```bash
python main.py --config cfgs/pretrain_ampt_data.yaml --exp_name ampt --val_freq 10 --gpu 0
```

- `--config`: Path to the configuration file
- `--exp_name`: Name of the experiment
- `--val_freq`: Validation frequency during training
- `--gpu`: Gpu device number
### Fineturning

To finetune a pre-trained MaskPoint-AMPT model, simply run like:
```bash
python main.py --config cfgs/finetune_ampt_data.yaml --finetune_model --ckpts experiments/pretrain_ampt_data/ampt/ckpt-last.pth --exp_name ampt_data --gpu 0
```

## License
coming soon!

## Acknowledgments
Special thanks to the developers of PointNet++, MaskPoint and the authors of related datasets. and the original MaskPoint can be found in https://github.com/WisconsinAIVision/MaskPoint
@article{liu2022masked,
  title={Masked Discrimination for Self-Supervised Learning on Point Clouds},
  author={Liu, Haotian and Cai, Mu and Lee, Yong Jae},
  journal={Proceedings of the European Conference on Computer Vision (ECCV)},
  year={2022}
}
