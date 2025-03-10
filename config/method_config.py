"""
Method configuration for the benchmarkDA project.
Contains method-specific parameters and command templates.
"""

# Python method configurations
PYTHON_METHODS = {
    "mellon": {
        "script": "Mellon_bm.py",
        "params": {
            "mellon_d_method": "fractal",
            "norm_density": "No",
            "hyperparameter": "Yes",
            "corrected": "No",
            "ls_factor": 1.5
        }
    },
    "mellon_noNorm": {
        "script": "Mellon_bm.py",
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
        "params": {}
    },
    "meld_default": {
        "script": "meld_bm.py",
        "params": {
            "beta": 40
        }
    },
    "cna": {
        "script": "CNA_bm.py",
        "params": {}
    }
}

# R method configurations
R_METHODS = {
    "milo": {
        "method_name": "milo",
        "params": {}
    },
    "daseq": {
        "method_name": "daseq",
        "params": {}
    },
    "cydar": {
        "method_name": "cydar",
        "params": {}
    },
    "louvain": {
        "method_name": "louvain",
        "params": {}
    }
}

# Common method parameters
COMMON_PARAMS = {
    "layer_embedding_pca": "X_pca",
    "layer_embedding_dm": "DM_EigenVectors",
}

def get_python_method_cmd(method, file_path, pop, pop_enr, pop_col, ds_type, batch_sd, 
                         input_file, seed, layer_embedding, output_dir, k=30, beta=None, 
                         n_dm=0, ls_mode="PCA"):
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
        cmd += f"    --n_dm {n_dm} \\\n"
        cmd += f"    --mellon_d_method \"{method_params.get('mellon_d_method', 'fractal')}\" \\\n"
        cmd += f"    --norm_density \"{method_params.get('norm_density', 'No')}\" \\\n"
        cmd += f"    --hyperparameter \"{method_params.get('hyperparameter', 'Yes')}\" \\\n"
        cmd += f"    --corrected \"{method_params.get('corrected', 'No')}\" \\\n"
        cmd += f"    --ls_factor {method_params.get('ls_factor', 1.5)} \\\n"
        cmd += f"    --ls_mode {ls_mode} \\\n"
    elif script == "meld_bm.py":
        cmd += f"    --beta {beta if beta else method_params.get('beta', 40)} \\\n"
        cmd += f"    --k_meld {k} \\\n"
    elif script == "CNA_bm.py":
        cmd += f"    --k_cna {k} \\\n"
    
    cmd += f"    --output_dir {output_dir}/"
    
    return cmd

def get_r_method_cmd(method, data_file, pop, pop_enr, pop_col, ds_type, batch_sd,
                    input_file, seed, layer_embedding, output_dir, k=30, resolution=0.5):
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
    cmd += f"    --output_dir {output_dir}/"
    
    return cmd