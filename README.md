# Reproducibility study for the [paper](https://ieeexplore.ieee.org/document/9101598) "SuRF: Identification of Interesting Data Regions with Surrogate Models"

experiments, from the original SuRF [repository](https://github.com/Skeftical/SuRF-Reproducibility). This is done to ensure a straightforward comparison of our results with the original ones. Therefore, repeating our experiments requires cloning this repository and installing the necessary dependencies. Copy the files from the current repository into the same folder. You will also need to have [R](https://cran.r-project.org/) and its package ['prim'](https://cran.r-project.org/web/packages/prim/index.html) installed.

Perform the following steps:

To generate datasets and models
- Run the Synthetic Data Generation notebook. This will create the base synthetic datasets in input
- Run the python codeabase/query_generation.py script. To generate workloads with past region evaluations
- Run "python codebase/model_training.py" to train models on queries.

For accuracy experiments and plots
- Run "acc_experiments.R" to get and store the boxes found by PRIM.
- Run "Accuracy_Experiments.ipynb".

For performance experiments and plots
- Run Performance_Experiments.ipynb

For the experiments demonstrating PRIM interactivity
- Run Interactivity.ipynb

