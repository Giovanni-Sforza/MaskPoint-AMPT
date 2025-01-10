
# main
## import 
### tools
- tools/runner_pretrain -> run_net == pretrain_run_net
- tools/runner_finetune -> run_net == finetune_run_net
- tools/runner_finetune -> test_net == test_run_net
- tools/runner_finetune -> test_net_extract_feat 


### utils
- utils -> parse (get arg at terminal), dist_utils (train in distributed GPU, which i did not need)
- utils -> misc(fps seprate_point_cloud)
- utils -> logger(get log by logging in the terminal)
- utils -> config ***(get config from .yaml)*** (functions to manage, log, and manipulate configuration data, which is crucial for setting up experiments, managing hyperparameters, and ensuring reproducibility)

----
## prepare
- use distributed GPU ?
- logger
- tensorboard
- config
- batch size
- log_args/config_to_file
----
## run
- l92 -> 112 ( tools/runner_pretrain, tools/runner_finetune)

----

# tools/runner_pretrain
## import
### utils 
- utils -> misc, dist_utils, logger, AverageMeter
- datasets -> data_transforms




-------
# How the Project Loads Data
The project "MaskPoint" is designed to work with 3D point cloud data, and it uses a self-supervised learning approach for tasks like 3D object detection and classification. The project is structured to handle datasets like ShapeNet and ModelNet, as indicated in the log file and the provided code snippets.



1. **Configuration Files**: The project uses configuration files (e.g., `cfgs/pretrain_shapenet.yaml`) to specify dataset paths, number of points, and other parameters. These configurations are logged at the start of the training process, as seen in the log file.

2. **Dataset Classes**: The project has specific Python classes for handling different datasets. For example:
   - `ModelNetDataset.py` and `ModelNetDatasetFewShot.py` handle the ModelNet dataset.
   - `ShapeNet55Dataset.py` handles the ShapeNet dataset.

3. **Data Loading Process**:
   - Each dataset class initializes with a configuration object that provides paths and parameters.
   - The classes read data files (e.g., `.txt` or `.pkl` files) to load point cloud data and labels.
   - The data is processed and stored in a format suitable for training, such as numpy arrays or PyTorch tensors.

4. **Data Preprocessing**: The project may perform preprocessing steps like sampling a fixed number of points from each point cloud, normalizing the data, or applying transformations.

5. **DataLoader**: The project likely uses PyTorch's `DataLoader` to handle batching and shuffling of data during training and evaluation.

### Changing the Dataset

To use a different dataset with this project, you would need to follow these steps:

1. **Prepare the Dataset**: Ensure your dataset is in a format that can be easily read and processed. This might involve converting your data into a format similar to the existing datasets (e.g., point cloud files with labels).

2. **Modify Configuration Files**: Update the configuration files to point to your new dataset paths and adjust any dataset-specific parameters (e.g., number of points, categories).

3. **Create a New Dataset Class**: If your dataset has a different structure, you might need to create a new dataset class similar to `ModelNetDataset.py` or `ShapeNet55Dataset.py`. This class should handle loading and preprocessing your dataset.

4. **Update Data Loading Logic**: Ensure that the data loading logic in the training scripts (e.g., `main.py`, `runner_finetune.py`) uses your new dataset class.

5. **Test the Integration**: Run the training or evaluation scripts to ensure that your dataset is being loaded correctly and that the model can train on it.

By following these steps, you can adapt the "MaskPoint" project to work with a different 3D point cloud dataset.

# Example: ShapeNet
The `ShapeNet` class in the provided code snippet is designed to load and process point cloud data from the ShapeNet dataset. The class inherits from `torch.utils.data.Dataset`, which is a common interface for working with datasets in PyTorch.

Here's a breakdown of how the `ShapeNet` class works:

1. Initialization: The `__init__` method initializes the class with the provided configuration parameters. It sets the data root directory, point cloud path, subset (train or test), number of points, and other relevant attributes.

2. Loading the file list: The class reads the list of files for the specified subset from a text file (`self.data_list_file`). It also opens a separate file for the test data list (`test_data_list_file`) and appends the test data to the file list if the `whole` attribute is set to `True`.

3. Processing each instance: The `__getitem__` method is responsible for loading and processing each instance in the dataset. It retrieves the taxonomy ID, model ID, and file path for the current instance. It then reads the point cloud data from the specified file path, performs random sampling to reduce the number of points, normalizes the point cloud, and converts the data to a PyTorch tensor.

4. Returning the data: The `__getitem__` method returns the taxonomy ID, model ID, and the processed point cloud data as a tuple.

5. Length of the dataset: The `__len__` method returns the total number of instances in the dataset.

Regarding the identification of the label of the points in the ShapeNet dataset, the provided code snippet does not explicitly handle labeling. The ShapeNet dataset is a collection of 3D models, and the labels are typically associated with the categories or classes of the models, not the individual points within the models.

If you need to work with labeled point cloud data from the ShapeNet dataset, you would need to find additional datasets or annotations that provide point-level labels. The ShapeNet dataset itself does not provide point-level labels.

If you are interested in working with labeled point cloud data, you might consider using datasets like the ModelNet40 or ScanObjectNN, which do provide point-level labels. These datasets are often used for benchmarking point cloud classification and segmentation tasks.