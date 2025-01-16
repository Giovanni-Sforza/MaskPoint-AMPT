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

   Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the installation script:
   ```bash
   bash install.sh
   ```

4. Verify the installation by running a test script:
   ```bash
   python test.py
   ```

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

