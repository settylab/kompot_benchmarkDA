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

- Synthetic datasets are available under the `data` directory.
- Real datasets are available at:https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH

## Note about dataset embedding layer information

More note about diffusion map datasets
- for synthetic datasets, I used dm = 10 for running scripts
- for real datasets, I used the dm = 30 (only limit to single cell rna seq datasets)
- for CYTOF datasets, the number of features are too small so there is no dimensional reduction applied now. We can put a small diffusion map component number on it later.

## Usage

*Note:* Our implementation can only be used on a cluster with Slurm job scheduler since we need to run thousands of jobs.

The benchmarking scripts are all located in the `bin` drectory.

```text
Original scripts from benchmarkDA for benchmarking has been saved at bin/original_bash_scripts_from_benchmarkDA
```

To run a benchmarking job, use the following command:

### Data preprocessing

#### For transfering from rds to anndata:

```text
Rscript scripts/convert_seurat.R input_path_rds output_path_h5ad
```

The name of input_path_rds and output_path_h5ad should be the same.

#### For data preprocessing

```sh
bash bin/dataset_preprocessing.sh $1 $2 $3 $4
```
- $1 :dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: name of embedding layer in anndata obsm, if the name of pca layer is X_pca, enter X_pca. 
- $3: n_dm, number of diffusion map component, if benchmarking on pca, enter 0
- $4: mode_embedding, if want to benchmarking on pca, enter "PCA"; if want to benchmarking on diffusion map, enter "DM"

Example:
```sh
bash bin/dataset_preprocessing.sh linear X_pca 10 DM
```

### Generate synthetic labels
PCA
```sh
bash bin/modified_benchmarkda.sh $1 $2 $3 $4 $5
```
- $1 :dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: for indicating whether analysis is on PCA or Diffuson map, can be any name, but just need to be consistent with further analysis
- $3: for showing whether the ground truth is balanced or not. Using "No" now.
- $4: mode_embedding, if want to benchmarking on pca, enter "PCA"; if want to benchmarking on diffusion map, enter "DM"
- $5: n_dm, number of diffusion map component, if benchmarking on pca, enter 0

Example:
```sh
bash bin/modified_benchmarkda.sh linear pca No PCA 0
```
Diffusion Map

```sh
bash bin/modified_benchmarkda_dm_all.sh $1 $2 $3 $4 $5
```
- Input similar as PCA

Example:
```sh
bash bin/modified_benchmarkda_dm_all.sh linear dm No DM 10
```

### Running benchmarking
PCA

Python packages
- For synthetic datasets
```sh
bash bin/bm_syn_python_synthetic_pca_final.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```
- For real datasets
```sh
bash bin/bm_syn_python_real_pca_final.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```

R packages
- For all datasets
```sh
bash bin/R_method_evaluation_pca.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```

- $1 :dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: for indicating whether analysis is on PCA or Diffuson map, can be any name, need to match with what set in the "Generate synthetic labels" step.
- $3: iteration number : number of iterations if we select centroid randomly, now just use 0
- $4: for showing whether the ground truth is balanced or not. Using "No" now.
- $5: mellon d method: Use "fractal"
- $6: whether Mellon density is normalized or not: Use "No", un-normalized
- $7: whether apply correction on Mellon density log fold change: Use "No", no correction
- $8: whether Mellon parameters for density estimation are synchronized: Use "Yes", synchronized.
- $9: number of diffusion map components, if pca, n_dm = 0; if dm, the n_dm should match with previous dm
- $10: if benchmarking on pca, enter "PCA"; if benchmarking on dm, enter "DM".

Example:
```sh
bash bin/bm_syn_python_synthetic_pca_final.sh linear pca 0 No fractal No No Yes 0 PCA
```

```sh
bash bin/bm_syn_python_real_pca_final.sh covid19-pbmc pca 0 No fractal No No Yes 0 PCA
```

```sh
bash bin/R_method_evaluation_pca.sh linear pca 0 No fractal No No Yes 0 PCA
```


Diffusion Map


Python packages
- For synthetic datasets
```sh
bash bin/bm_syn_python_synthetic_dm.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```
- For real datasets
```sh
bash bin/bm_syn_python_real_pca_dm.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```

R packages
- For all datasets
```sh
bash bin/R_method_evaluation_dm.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
```

Example:
```sh
bash bin/bm_syn_python_synthetic_dm.sh linear dm 0 No fractal No No Yes 10 DM
```


The parameter settings are same. Only path for saving results are different. 

For running bash scripts, please make sure that the jobid match (${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}) with where we save ground truth.

### Helper notebooks

Notebooks for calculating diffusion map, converting anndata and rds, and perform evaluation are saved in the notebook/.

### Codes which will be updated soon:
- Code for generating the simulated batch effect based on the added scaler on batch standard deviation
- Code for Meld beta parameter optimization

## Acknowledgement

Our implementation is inspired by the repo https://github.com/MarioniLab/milo_analysis_2020.
