import sys

import os


import argparse


def main():

    parser = argparse.ArgumentParser(description='Your script description.')
    parser.add_argument('--input_file_path', type=str, help='Path to the input anndata file')
    parser.add_argument('--output_file_path', type=str, help='Path to the output rds file')

    args = parser.parse_args()
    input_file_path = args.input_file_path
    output_file_path = args.output_file_path


    if "/app/software/R" in os.environ["PATH"]:
            print(
            "R is loaded in the environment. Please make sure this is not part of "
            "your jupyter server start script and restart it from a clean command "
            "line environment."
        )
    else:
            print("R does not seem to be loaded and you are good to go.")
            


            # Make sure no R module is loaded
    if "/app/software/R" in os.environ["PATH"]:
            raise Exception("An R module seems to be loaded.")

            # Get the path to the python executable
    python_executable_path = sys.executable

            # Extract the path to the environment from the path to the python executable
    env_path = os.path.dirname(os.path.dirname(python_executable_path))

    print(
            f"Micromamba env path: {env_path}\n"
            "Please make sure you have R installed in the micromamba environment."
        )
    print(env_path)
    os.environ['R_HOME'] = os.path.join(env_path, 'lib', 'R')
        
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, r, ListVector

        #import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter


    # Function to get the library paths used by R
    def get_r_lib_paths():
        r_command = """
        .libPaths()
        """
        return robjects.r(r_command)

        # Execute the function and store the result
    r_lib_paths = get_r_lib_paths()

        # Print the first library path used by R (assuming there's at least one path)
    if r_lib_paths is not None:
        print("R is using packages from " + r_lib_paths[0])
    else:
        print("No library paths found.")
        
    import anndata
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    anndata2ri.activate()

    adata = anndata.read_h5ad(input_file_path)

    r_adata = ro.conversion.py2rpy(adata)
    ro.globalenv["r_adata"] = r_adata
    r_command = f'saveRDS(r_adata, file="{output_file_path}")'
    ro.r(r_command)
    print("Successfully save new rds file")
if __name__ == "__main__":
    main()