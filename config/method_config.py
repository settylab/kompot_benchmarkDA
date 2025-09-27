"""
Method configuration for the benchmarkDA project.
Contains method-specific parameters and command templates.
"""

# Python method configurations
PYTHON_METHODS = {
    "mellon": {
        "script": "Mellon_bm.py",
        "description": "Mellon with standard fractal density estimation and hyperparameter sync",
        "params": {
            "mellon_d_method": "fractal",
            "norm_density": "No",
            "hyperparameter": "Yes",
            "corrected": "No",
            "ls_factor": 1.5
        }
    },
    "mellon_noSync": {
        "script": "Mellon_bm.py",
        "description": "Mellon without hyperparameter synchronization",
        "params": {
            "mellon_d_method": "fractal",
            "norm_density": "No",
            "hyperparameter": "No",
            "corrected": "No",
            "ls_factor": 1.5
        }
    },
    "mellon_corr": {
        "script": "Mellon_bm.py",
        "description": "Mellon with batch correction",
        "params": {
            "mellon_d_method": "fractal",
            "norm_density": "No",
            "hyperparameter": "Yes",
            "corrected": "Yes",
            "ls_factor": 1.5
        }
    },
    "meld": {
        "script": "meld_bm.py",
        "description": "MELD with dataset-specific parameters",
        "params": {}
    },
    "meld_default": {
        "script": "meld_bm.py",
        "description": "MELD with fixed beta=40 parameter",
        "params": {
            "beta": 40
        }
    },
    # CNA methods temporarily disabled due to NumPy 2.0 compatibility issues
    # "cna": {
    #     "script": "CNA_bm.py",
    #     "description": "Conditional Neighborhood Analysis",
    #     "params": {}
    # },
    "kompot": {
        "script": "kompot_bm.py",
        "description": "Kompot differential abundance testing",
        "params": {
            "ls_factor": 10.0,
            "n_landmarks": None,
            "log_fold_change_threshold": 1.0,
            "pvalue_threshold": 0.05
        }
    },
    "kompot_pca": {
        "script": "kompot_bm.py",
        "description": "Kompot using PCA embedding",
        "params": {
            "ls_factor": 10.0,
            "n_landmarks": None,
            "log_fold_change_threshold": 1.0,
            "pvalue_threshold": 0.05,
            "force_pca": True
        }
    },
    "mellon_pca": {
        "script": "Mellon_bm.py",
        "description": "Mellon using PCA embedding",
        "params": {
            "mellon_d_method": "fractal",
            "norm_density": "No",
            "hyperparameter": "Yes",
            "corrected": "No",
            "ls_factor": 1.5,
            "force_pca": True
        }
    },
    "meld_pca": {
        "script": "meld_bm.py",
        "description": "MELD using PCA embedding",
        "params": {
            "force_pca": True
        }
    }
    # CNA PCA method temporarily disabled due to NumPy 2.0 compatibility issues
    # "cna_pca": {
    #     "script": "CNA_bm.py",
    #     "description": "CNA using PCA embedding",
    #     "params": {
    #         "force_pca": True
    #     }
    # }
}

# R method configurations
R_METHODS = {
    "milo": {
        "method_name": "milo",
        "description": "Neighborhood-based differential abundance testing",
        "params": {}
    },
    "daseq": {
        "method_name": "daseq",
        "description": "Differential abundance region detection",
        "params": {}
    },
    "cydar": {
        "method_name": "cydar",
        "description": "Hypersphere-based differential abundance testing",
        "params": {}
    },
    "louvain": {
        "method_name": "louvain",
        "description": "Clustering-based differential abundance testing",
        "params": {}
    }
}

# Common method parameters
COMMON_PARAMS = {
    "layer_embedding_pca": "X_pca",
    "layer_embedding_dm": "DM_EigenVectors",
}

def get_method_description(method_name):
    """Get description for a specific method."""
    if method_name in PYTHON_METHODS:
        return PYTHON_METHODS[method_name].get("description", "No description available")
    elif method_name in R_METHODS:
        return R_METHODS[method_name].get("description", "No description available")
    else:
        return "Unknown method"

def list_methods_with_descriptions():
    """Get formatted list of all methods with descriptions."""
    output = []

    output.append("PYTHON METHODS:")
    for name, config in PYTHON_METHODS.items():
        desc = config.get("description", "No description")
        output.append(f"  {name:<15} - {desc}")

    output.append("\nR METHODS:")
    for name, config in R_METHODS.items():
        desc = config.get("description", "No description")
        output.append(f"  {name:<15} - {desc}")

    return "\n".join(output)

def get_python_method_cmd(method, file_path, pop, pop_enr, pop_col, ds_type, batch_sd,
                         input_file, seed, layer_embedding, output_dir, k=30, beta=None,
                         n_dm=10, ls_mode="DM"):
    """Generate a Python method command with appropriate parameters."""
    
    method_config = PYTHON_METHODS[method]
    script = method_config["script"]
    method_params = method_config["params"]
    
    cmd = f"python {script} \\\n"
    cmd += f"    --file_path {file_path} \\\n"
    cmd += f"    --pop {pop} \\\n"
    cmd += f"    --pop_enr {pop_enr} \\\n"
    cmd += f"    --pop_column {pop_col} \\\n"
    cmd += f"    --ds_type {ds_type} \\\n"
    cmd += f"    --batch_sd {batch_sd} \\\n"
    cmd += f"    --input_file {input_file} \\\n"
    cmd += f"    --package {method} \\\n"
    cmd += f"    --seed {seed} \\\n"
    cmd += f"    --layer_embedding {layer_embedding} \\\n"
    
    # Add method-specific parameters
    if script == "Mellon_bm.py":
        # Force PCA mode if specified
        if method_params.get('force_pca'):
            n_dm = 0
        cmd += f"    --n_dm {n_dm} \\\n"
        cmd += f"    --mellon_d_method \"{method_params.get('mellon_d_method', 'fractal')}\" \\\n"
        cmd += f"    --norm_density \"{method_params.get('norm_density', 'No')}\" \\\n"
        cmd += f"    --hyperparameter \"{method_params.get('hyperparameter', 'Yes')}\" \\\n"
        cmd += f"    --corrected \"{method_params.get('corrected', 'No')}\" \\\n"
        cmd += f"    --ls_factor {method_params.get('ls_factor', 1.5)} \\\n"
        cmd += f"    --ls_mode {ls_mode} \\\n"
    elif script == "meld_bm.py":
        # Force PCA mode if specified
        if method_params.get('force_pca'):
            n_dm = 0
        cmd += f"    --beta {beta if beta else method_params.get('beta', 40)} \\\n"
        cmd += f"    --k_meld {k} \\\n"
        cmd += f"    --n_dm {n_dm} \\\n"
    elif script == "CNA_bm.py":
        # Force PCA mode if specified
        if method_params.get('force_pca'):
            n_dm = 0
        cmd += f"    --k_cna {k} \\\n"
        cmd += f"    --n_dm {n_dm} \\\n"
    elif script == "kompot_bm.py":
        cmd += f"    --n_dm {n_dm} \\\n"
        cmd += f"    --ls_factor {method_params.get('ls_factor', 10.0)} \\\n"
        if method_params.get('n_landmarks') is not None:
            cmd += f"    --n_landmarks {method_params.get('n_landmarks')} \\\n"
        cmd += f"    --log_fold_change_threshold {method_params.get('log_fold_change_threshold', 1.0)} \\\n"
        cmd += f"    --pvalue_threshold {method_params.get('pvalue_threshold', 0.05)} \\\n"
        if method_params.get('force_pca'):
            cmd += f"    --force_pca \\\n"

    cmd += f"    --output_dir {output_dir}/"

    return cmd

def get_r_method_cmd(method, data_file, pop, pop_enr, pop_col, ds_type, batch_sd,
                    input_file, seed, layer_embedding, output_dir, k=30, resolution=0.5,
                    n_dm=10, scripts_dir="${root}/scripts"):
    """Generate an R method command with appropriate parameters."""

    method_config = R_METHODS[method]
    method_name = method_config["method_name"]

    cmd = f"Rscript {scripts_dir}/run_DA.r \\\n"
    cmd += f"    --file_path {data_file} \\\n"
    cmd += f"    --pop {pop} \\\n"
    cmd += f"    --pop_enr {pop_enr} \\\n"
    cmd += f"    --pop_column {pop_col} \\\n"
    cmd += f"    --ds_type {ds_type} \\\n"
    cmd += f"    --batch_sd {batch_sd} \\\n"
    cmd += f"    --input_file {input_file} \\\n"
    cmd += f"    --package {method_name} \\\n"
    cmd += f"    --seed {seed} \\\n"
    cmd += f"    --layer_embedding {layer_embedding} \\\n"
    cmd += f"    --k {k} \\\n"
    cmd += f"    --resolution {resolution} \\\n"
    cmd += f"    --n_dm {n_dm} \\\n"
    cmd += f"    --output_dir {output_dir}/"
    
    return cmd