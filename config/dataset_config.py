"""
Dataset configuration for the benchmarkDA project.
Contains dataset-specific parameters used across different scripts.
"""

# Dataset configurations
DATASET_CONFIGS = {
    "cluster": {
        "pops": ["M1", "M2", "M3"],
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 0.2,
        "beta": 33,
        "downsample": 3,
        "pop_col": "celltype",
        "n_dm": 10,  # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use
    },
    "linear": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7"],
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 1,
        "beta": 71,
        "downsample": 3,
        "pop_col": "celltype",
        "n_dm": 10,  # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use
    },
    "branch": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"],
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 1,
        "beta": 65,
        "downsample": 3,
        "pop_col": "celltype",
        "n_dm": 10,  # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use
    },
    "covid19-pbmc": {
        "pops": ["PB", "CD14_Monocyte", "CD8_T", "CD4_T", "Platelet", "NK", "Granulocyte",
                "CD16_Monocyte", "gd_T", "pDC", "DC"],
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 0.5,
        "beta": 40,
        "downsample": 3,
        "pop_col": "cell.type.coarse",
        "n_dm": 30,  # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use
    },
    "bcr-xl": {
        "pops": ["CD4_T-cells", "CD8_T-cells", "NK_cells", "B-cells_IgM+", "B-cells_IgM-",
                "monocytes", "DC"],
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 0.5,
        "beta": 40,
        "downsample": 3,
        "pop_col": "cell_type",
        "n_dm": 5,   # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use (preprocessing creates X_pca from X_x for CyTOF)
    },
    "levine32": {
        "pops": ["pDCs", "CD4_T_cells", "CD8_T_cells", "Pre_B_cells", "Mature_B_cells", "Monocytes", "Basophils"],  # Major immune cell populations (7 total)
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 0.5,
        "beta": 40,
        "downsample": 3,
        "pop_col": "cell_type",
        "n_dm": 5,   # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use (preprocessing creates X_pca from X_x for CyTOF)
    },
    "pancreas": {
        "pops": ["alpha_cell", "beta_cell", "delta_cell", "gamma_cell", "ductal_cell",
                "acinar_cell", "endothelial_cell", "PSC_cell"],  # Major pancreatic cell types
        "batch_vec": [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],  # Variance-based batch effect strengths
        "k": 30,
        "resolution": 0.5,
        "beta": 40,
        "downsample": 3,
        "pop_col": "Factor.Value.inferred.cell.type...authors.labels.",
        "n_dm": 30,  # Diffusion map components for this dataset
        "layer_embedding": "X_pca"  # Embedding layer to use
    }
}

# Common parameter sets
SEEDS = [43, 44, 45]
ENRICHMENT_VALUES = [0.75, 0.85, 0.95]

# Helper functions
def get_n_dm_for_dataset(dataset):
    """Get the number of diffusion map components for a specific dataset."""
    if dataset in DATASET_CONFIGS:
        return DATASET_CONFIGS[dataset]["n_dm"]
    else:
        # Default fallback
        return 10


# Path templates
def get_data_file_path(root, data_id, mode_embedding, n_dm):
    """Get the path to the data file."""
    if data_id in ["linear", "branch", "cluster"]:
        return f"{root}/data/synthetic/{data_id}/{data_id}_{mode_embedding}_{n_dm}.h5ad"
    else:
        return f"{root}/data/{data_id}/{data_id}_{mode_embedding}_{n_dm}.h5ad"

def get_job_id(data_id, pop, enr, seed, batch_sd_num, balance_bool, analysis_layer):
    """Generate a consistent job ID across scripts."""
    return f"{data_id}-{pop}-{enr}-{seed}-{batch_sd_num}-{balance_bool}-{analysis_layer}"

def get_input_dir(root, data_id, job_id):
    """Get the input directory path."""
    if data_id in ["linear", "branch", "cluster"]:
        return f"{root}/data/synthetic/{data_id}/{job_id}/"
    else:
        return f"{root}/data/{data_id}/{job_id}/"

def get_output_dir(root, data_id, job_id, mode):
    """Get the output directory path based on mode (pca or dm)."""
    if mode.upper() == "PCA":
        return f"{root}/benchmark_pca/{'synthetic' if data_id in ['linear', 'branch', 'cluster'] else 'real'}/{data_id}/{job_id}"
    else:
        return f"{root}/benchmark_dm/{'synthetic' if data_id in ['linear', 'branch', 'cluster'] else 'real'}/{data_id}/{job_id}"

# Default SLURM options (can be overridden via CLI with --sbatch-options)
# These are reasonable defaults that work across most SLURM configurations
DEFAULT_SLURM_OPTIONS = {
    "time": "6:00:00",
    "mem": "32G"
}