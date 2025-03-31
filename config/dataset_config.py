"""
Dataset configuration for the benchmarkDA project.
Contains dataset-specific parameters used across different scripts.
"""

# Dataset configurations
DATASET_CONFIGS_PCA = {
    "cluster": {
        "pops": ["M1", "M2", "M3"],
        "batch_vec_orig": [0, 0.75, 1, 1.25, 1.5],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 0.2,
        "beta": 33,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "linear": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7"],
        #"pops": ["M1", "M2"],
        "batch_vec_orig": [0, 0.75, 1, 1.25, 1.5],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 1,
        "beta": 71,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "branch": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"],
        "batch_vec_orig": [0, 0.75, 1, 1.25, 1.5],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 1,
        "beta": 65,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "covid19-pbmc": {
        "pops": ["RBC","B","PB", "CD14_Monocyte", "CD8_T", "CD4_T", "Platelet", "NK", "Granulocyte", 
                "CD16_Monocyte", "gd_T", "pDC", "DC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30, 
        "resolution": 0.5,
        "beta": 25,
        "downsample": 3,
        "pop_col": "cell.type.coarse"
    },
    "bcr-xl": {
        "pops": ["CD4_T-cells", "NK_cells", "CD8_T-cells", "B-cells_IgM+", "monocytes", "surface-", "B-cells_IgM-", "DC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 0.6,
        "beta": 23,
        "downsample": 10,
        "pop_col": "cell_type"
    },
    "pancreas":{
        "pops": ["delta_cell", "alpha_cell", "gamma_cell", "acinar_cell", "beta_cell", "ductal_cell", "epsilon_cell"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 1.2,
        "beta": 80,
        "downsample": 3,
        "pop_col": "Factor.Value.inferred.cell.type...authors.labels."
    },
    "levine32":{
        "pops": ["pDCs" ,"CD4_T_cells" ,"CD8_T_cells","Pre_B_cells","Mature_B_cells","Monocytes", "Basophils"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 0.6,
        "beta": 36,
        "downsample": 25,
        "pop_col": "cell_type"
    },
    "aging":{
        "pops": ["CLP","Ery_P","HSC","ILC","Immature_B_cell","LMPP", "MBE","MKP","Mature_B_cell" ,"Mono_P" ,"Monocyte" ,"Myelo_P" ,"NK","Neutrophil","Pre-B_cell","T_cell","Treg","cDC", "pDC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.07, 0.1, 0.14, 0.19, 0.27, 0.38, 0.53, 0.74, 1.03],
        "k": 30,
        "resolution": 1,
        "beta": 64,
        "downsample": 3,
        "pop_col": "midres_celltype_benchmarking"
    }


}

DATASET_CONFIGS_DM = {
    "cluster": {
        "pops": ["M1", "M2", "M3"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 0.2,
        "beta": 33,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "linear": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7"],
        #"pops": ["M1", "M2"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 1,
        "beta": 71,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "branch": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 1,
        "beta": 65,
        "downsample": 3,
        "pop_col": "celltype"
    },
    "covid19-pbmc": {
        "pops": ["RBC","B","PB", "CD14_Monocyte", "CD8_T", "CD4_T", "Platelet", "NK", "Granulocyte", 
                "CD16_Monocyte", "gd_T", "pDC", "DC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30, 
        "resolution": 0.5,
        "beta": 25,
        "downsample": 3,
        "pop_col": "cell.type.coarse"
    },
    "bcr-xl": {
        "pops": ["CD4_T-cells", "NK_cells", "CD8_T-cells", "B-cells_IgM+", "monocytes", "surface-", "B-cells_IgM-", "DC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 0.6,
        "beta": 23,
        "downsample": 10,
        "pop_col": "cell_type"
    },
    "pancreas":{
        "pops": ["delta_cell", "alpha_cell", "gamma_cell", "acinar_cell", "beta_cell", "ductal_cell", "epsilon_cell"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 1.2,
        "beta": 80,
        "downsample": 3,
        "pop_col": "Factor.Value.inferred.cell.type...authors.labels."
    },
    "levine32":{
        "pops": ["pDCs" ,"CD4_T_cells" ,"CD8_T_cells","Pre_B_cells","Mature_B_cells","Monocytes", "Basophils"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 0.6,
        "beta": 36,
        "downsample": 25,
        "pop_col": "cell_type"
    },
    "aging":{
        "pops": ["CLP","Ery_P","HSC","ILC","Immature_B_cell","LMPP", "MBE","MKP","Mature_B_cell" ,"Mono_P" ,"Monocyte" ,"Myelo_P" ,"NK","Neutrophil","Pre-B_cell","T_cell","Treg","cDC", "pDC"],
        "batch_vec_orig": [0],
        "batch_vec_modified": [0, 0.05, 0.06, 0.08, 0.11, 0.14, 0.19, 0.24, 0.31, 0.41, 0.53],
        "k": 30,
        "resolution": 1,
        "beta": 64,
        "downsample": 3,
        "pop_col": "midres_celltype_benchmarking"
    }
}

# Common parameter sets
SEEDS = [43, 44, 45]
ENRICHMENT_VALUES = [0.75, 0.85, 0.95]

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