# bmDA

The implementation of benchmarking the performance of differential abundance (DA) testing methods including 5 clustering-free methods:

1. [Testing for differential abundance in mass cytometry data (Cydar)](https://www.nature.com/articles/nmeth.4295)
2. [Detection of differentially abundant cell subpopulations in scRNA-seq data (DAseq)](https://www.pnas.org/doi/abs/10.1073/pnas.2100293118)
3. [Quantifying the effect of experimental perturbations at single-cell resolution (MELD)](https://www.nature.com/articles/s41587-020-00803-5)
4. [Differential abundance testing on single-cell data using k-nearest neighbor graphs (Milo)](https://www.nature.com/articles/s41587-021-01033-z)
5. [Co-varying neighborhood analysis identifies cell populations associated with phenotypes of interest from single-cell transcriptomics (CNA)](https://www.nature.com/articles/s41587-021-01066-4),

and a clustering-based method, *Louvain*.

## Dependencies

To run this benchmarking codes, it needs to install a list of R and Python packages. The R packages needed are:

- argparse
- SingleCellExperiment
- scran
- [DAseq](https://github.com/KlugerLab/DAseq)
- [miloR](https://github.com/MarioniLab/miloR)
- tibble
- dplyr
- tidyverse
- igraph
- [cydar](http://bioconductor.org/packages/cydar)
- pdist
- reshape2

The Python packages needed are

- [MELD](https://github.com/KrishnaswamyLab/MELD)
- [cna](https://github.com/immunogenomics/cna)
- scanpy
- [graphtools](https://github.com/KrishnaswamyLab/graphtools)
- scikit-learn
- multianndata

The python packages list are saved in the differential_abundance_env_list.yml for building environment.


The R packages list are saved in the renv.lock.
- Use module load R/4.3.1-gfbf-2022b first and then open R to create the renv environemnt by renv::restore()

## Data

- Synthetic datasets and BCR-XL dataset are available under the `data` directory.
- The COVID-19 PBMC dataset is available at https://www.covid19cellatlas.org/#wilk20.

But for now, all datasets are labeled with the link where they saved directly in bash scripts.
- For finding datasets:
   1. Please open bin/bm_syn_python_real_pca_final.sh or bin/bm_syn_python_synthetic_pca_final.sh for finding datasets with PCA layer saved in h5ad format.
    2. Please open bin/bm_syn_python_real_dm_final.sh or bin/bm_syn_python_synthetic_dm_final.sh for finding datasets with Diffusion Map layer saved in h5ad format.
    3. Please open bin/R_method_evaluation_dm.sh for finding datasets with Diffusion Map layer saved in RDS format.
    4. Please open bin/R_method_evaluation_pca.sh for finding datasets with PCA layer saved in RDS format.

For those datasets which has "dm", it means that its X_pca or PCA has been replaced by diffusion map values.

More note about diffusion map datasets
- for synthetic datasets, if want dm=30 instead of 10, please delete "_10" in the file name
- for real datasets, the dm = 30 (only limit to single cell rna seq datasets)
- for CYTOF datasets, the number of features are too small so there is no dimensional reduction applied now. We can put a small diffusion map component number on it later.

## Usage

*Note:* Our implementation can only be used on a cluster with Slurm job scheduler since we need to run thousands of jobs.

The benchmarking scripts are all located in the `bin` drectory.

```text
Original scripts from benchmarkDA for benchmarking has been saved at bin/original_bash_scripts_from_benchmarkDA
```

To run a benchmarking job, use the following command:

### Generate synthetic labels
PCA
```sh
bash bin/modified_benchmarkda.sh $1 $2 $3
```
- $1 :dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: for indicating whether analysis is on PCA or Diffuson map, can be any name, but just need to be consistent with further analysis
- $3: for showing whether the ground truth is balanced or not. Using "No" now.

Example:
```sh
bash bin/modified_benchmarkda.sh linear pca No
```
Diffusion Map

```sh
bash bin/modified_benchmarkda_dm_all.sh $1 $2 $3
```
- Input similar as PCA

Example:
```sh
bash bin/modified_benchmarkda_dm_all.sh linear dm No
```

### Running benchmarking
PCA

Python packages
- For synthetic datasets
```sh
bash bin/bm_syn_python_synthetic_pca_final.sh $1 $2 $3 $4 $5 $6 $7 $8
```
- For real datasets
```sh
bash bin/bm_syn_python_real_pca_final.sh $1 $2 $3 $4 $5 $6 $7 $8
```

R packages
- For all datasets
```sh
bash bin/R_method_evaluation_pca.sh $1 $2 $3 $4 $5 $6 $7 $8
```

- $1 :dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: for indicating whether analysis is on PCA or Diffuson map, can be any name, need to match with what set in the "Generate synthetic labels" step.
- $3: iteration number : number of iterations if we select centroid randomly, now just use 0
- $4: for showing whether the ground truth is balanced or not. Using "No" now.
- $5: mellon d method: Use "fractal"
- $6: whether Mellon density is normalized or not: Use "No", un-normalized
- $7: whether apply correction on Mellon density log fold change: Use "No", no correction
- $8: whether Mellon parameters for density estimation are synchronized: Use "Yes", synchronized.

Example:
```sh
bash bin/bm_syn_python_synthetic_pca_final.sh linear pca 0 No fractal No No Yes
```

```sh
bash bin/bm_syn_python_real_pca_final.sh covid19-pbmc pca 0 No fractal No No Yes
```

```sh
bash bin/R_method_evaluation_pca.sh linear pca 0 No fractal No No Yes
```


Diffusion Map


Python packages
- For synthetic datasets
```sh
bash bin/bm_syn_python_synthetic_dm.sh $1 $2 $3 $4 $5 $6 $7 $8
```
- For real datasets
```sh
bash bin/bm_syn_python_real_pca_dm $1 $2 $3 $4 $5 $6 $7 $8
```

R packages
- For all datasets
```sh
bash bin/R_method_evaluation_dm.sh $1 $2 $3 $4 $5 $6 $7 $8
```
The parameter settings are same. Only path for saving results are different. 

For running bash scripts, please make sure that the jobid match (${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}) with where we save ground truth.

### Helper notebooks

Notebooks for calculating diffusion map, converting anndata and rds, and perform evaluation are saved in the notebook/.

### Codes which will be updated soon:
- Code for generating the simulated batch effect based on the shape scale of embedding layer
- Code for Meld beta parameter optimization

## Acknowledgement

Our implementation is inspired by the repo https://github.com/MarioniLab/milo_analysis_2020.
