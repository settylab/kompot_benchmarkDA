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

## Environment Setup

### Using Micromamba (Recommended)
The project now supports using micromamba instead of conda for environment management. Micromamba is a tiny version of mamba, which is a fast, robust package manager compatible with conda packages.

1. Install micromamba:
   ```
   curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
   ```

2. Create the environment:
   ```
   eval "$(micromamba shell hook --shell bash)"
   micromamba create -f differential_abundance_env_list.yml -n DiffAbundance
   ```

3. Activate the environment:
   ```
   micromamba activate DiffAbundance
   ```

Note: When running jobs on a cluster, make sure to load required modules before activating the micromamba environment, e.g.:
```
module load cuDNN/8.4.1.50-CUDA-11.7.0
```

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

### Architecture

The benchmarking implementation has been refactored for improved readability, reliability, and reduced redundancy:

1. **Configuration System**:
   - `config/dataset_config.py`: Central repository for dataset parameters
   - `config/method_config.py`: Method-specific parameters and command generation

2. **Unified Data Loading**:
   - `python_method/data_loader.py`: Common data loading patterns for all methods

3. **Standardized Method Implementation**:
   - Consistent Python method interfaces for Mellon, MELD, and CNA
   - Parameterized method configuration to reduce duplication

4. **Unified Script Generation**:
   - `bin/run_benchmark.py`: Single script generator for all benchmark runs

### Data preprocessing

#### For transfering from rds to anndata:

```text
Rscript scripts/convert_seurat.R input_file_rds output_file_h5ad
```
Or using convert2anndata directly:

```text
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input_file.rds -o /path/to/output_file.h5ad
```

NOTE: The name and path of input_file_rds and output_file_h5ad should be the same!

For covid19-pbmc dataset: Needs to use UpdateSeuratObject(seurat_obj) before sce <- convert_seurat_to_sce(seurat_obj)

#### For data preprocessing

```sh
bash bin/dataset_preprocessing.sh $1 $2 $3 $4
```
- $1: dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: name of embedding layer in anndata obsm, if the name of pca layer is X_pca, enter X_pca. 
- $3: n_dm, number of diffusion map component, if benchmarking on pca, enter 0
- $4: mode_embedding, if want to benchmarking on pca, enter "PCA"; if want to benchmarking on diffusion map, enter "DM"

Example:
```sh
bash bin/dataset_preprocessing.sh linear X_pca 10 DM
```

### Generate synthetic labels

```sh
bash bin/modified_benchmarkda.sh $1 $2 $3 $4 $5
```
- $1: dataset name (can be linear, branch, cluster, covid19-pbmc, bcr-xl, pancreas, aging, levine32)
- $2: for indicating whether analysis is on PCA or Diffuson map, can be any name, but just need to be consistent with further analysis
- $3: for showing whether the ground truth is balanced or not. Using "No" now.
- $4: mode_embedding, if want to benchmarking on pca, enter "PCA"; if want to benchmarking on diffusion map, enter "DM"
- $5: n_dm, number of diffusion map component, if benchmarking on pca, enter 0

Example:
```sh
# For PCA
bash bin/modified_benchmarkda.sh linear pca No PCA 0

# For Diffusion Map
bash bin/modified_benchmarkda_dm_all.sh linear dm No DM 10
```

### Running benchmarks

The benchmarking process has been simplified using a unified script generator. This allows running benchmark analysis for any combination of dataset, methods, and parameters.

#### Generate a benchmark script

```sh
python bin/run_benchmark.py --dataset linear --method_type python --mode_embedding PCA
```

Parameters:
- `--dataset`: Dataset name (e.g., linear, branch, cluster, covid19-pbmc)
- `--method_type`: Type of methods to run (python or r)
- `--analysis_layer`: Analysis layer name (default: pca)
- `--iteration_num`: Number of iterations (default: 0)
- `--balance`: Balance flag (default: No)
- `--n_dm`: Number of diffusion map components (default: 0)
- `--mode_embedding`: Embedding mode (PCA or DM, default: PCA)
- `--methods`: Specific methods to run (space-separated list)
- `--output`: Custom output script path

#### Examples

Generate a script for running Python methods on the linear dataset with PCA embedding:
```sh
python bin/run_benchmark.py --dataset linear --method_type python --mode_embedding PCA
```

Generate a script for running R methods on the covid19-pbmc dataset with diffusion map embedding:
```sh
python bin/run_benchmark.py --dataset covid19-pbmc --method_type r --mode_embedding DM --n_dm 30
```

Generate a script for running only the Mellon and MELD methods:
```sh
python bin/run_benchmark.py --dataset branch --method_type python --methods mellon meld
```

#### Run the generated script

The script generator will create an executable shell script that can be run directly:
```sh
bash run_benchmark_linear_pca_python.sh
```

For all benchmark scripts, ensure that the jobid matches (${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}) with where the ground truth was saved.

### Legacy Scripts

Original scripts are still available in the `bin` directory for backward compatibility:

```text
Original scripts from benchmarkDA for benchmarking have been saved at bin/original_bash_scripts_from_benchmarkDA
```

### Helper notebooks

Notebooks for calculating diffusion map, converting anndata and rds, and perform evaluation are saved in the notebook/.

### Codes which will be updated soon:
- Code for generating the simulated batch effect based on the added scaler on batch standard deviation
- Code for Meld beta parameter optimization

## Acknowledgement

Our implementation is inspired by the repo https://github.com/MarioniLab/milo_analysis_2020.
