# Reproducibility study for the [paper](https://ieeexplore.ieee.org/document/9101598) "SuRF: Identification of Interesting Data Regions with Surrogate Models"

This reproducibility study inherits some elements: dataset generation, SuRF-related experiments, from the original SuRF [repository](https://github.com/Skeftical/SuRF-Reproducibility). This is done to ensure a straightforward comparison of our results with the original.
Repeating our experiments requires cloning this repository and installing dependencies. Copy the files from the current repository into the same folder. You will also need to have [R](https://cran.r-project.org/) and its packages [prim](https://cran.r-project.org/web/packages/prim/index.html), [FNN](https://cran.r-project.org/web/packages/FNN/index.html), and [ks](https://cran.r-project.org/web/packages/ks/index.html) installed.

Perform the following steps:

Create folders required by the code from SuRF repository:
 - Run R script `create_folders.R` 

Generate datasets and models:
- Run the `Synthetic Data Generation-Copy1.ipynb` notebook. This will create the basic synthetic datasets in input folder.
- Run the python script `codeabase/query_generation.py`. To create workloads with past region evaluations.
- Run the python script `codebase/model_training.py` to train models on queries.

For accuracy experiments and plots
- Run `acc_experiments.R`
- Run `Accuracy_Experiments.ipynb`.

For performance experiments and plots
- Run `Performance_Experiments.ipynb`

For the experiments that demonstrate PRIM interactivity
- Run `Interactivity.ipynb`