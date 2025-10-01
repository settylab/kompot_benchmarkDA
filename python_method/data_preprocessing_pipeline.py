import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from importlib import reload
import sys
from pathlib import Path



import read_file


import calculate_diffusion_map
import preprocessing

import argparse
import numpy as np


import sys
import os



def main():

    parser = argparse.ArgumentParser(description='Your script description.')
    parser.add_argument('--file_path', type=str, help='Path to the input anndata file')
    parser.add_argument('--embedding_layer', type=str, help='name of the embedding layer')
    parser.add_argument('--n_dm', type=int, help='number of diffusion map components number')
    parser.add_argument('--mode_embedding', type=str, help='DM or PCA')
    parser.add_argument('--pop_col', type=str, help='column for celltypes')
    parser.add_argument('--output_dir', type=str, help='where to save the output anndata')


    args = parser.parse_args()
    file_path = args.file_path
    embedding_layer = args.embedding_layer
    n_dm = args.n_dm
    mode_embedding = args.mode_embedding
    pop_col = args.pop_col
    output_dir = args.output_dir
    


    # Load dataset
    print(f"  Loading: {file_path}")
    adata = read_file.read_dataset(file_path, embedding_layer)
    print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} features")

    ## detecting whether there is the space inside the cell names string, replace the space with the underline.
    ## Only do this if pop_col exists in the data (preprocessing embeddings doesn't require this)
    if pop_col in adata.obs.columns:
        adata = preprocessing.replace_space_in_string(adata, pop_col)
        print(f"  Cell type column: '{pop_col}' (spaces replaced with underscores)")

    if mode_embedding == "DM":
        print(f"  Computing diffusion map ({n_dm} components) from {embedding_layer}...")
        adata = calculate_diffusion_map.calculate_dm(adata, embedding_layer, n_dm)
        print(f"  ✓ DM_EigenVectors computed")
    else:
        print(f"  Embedding layer '{embedding_layer}' prepared")

    # Check for duplicates
    if preprocessing.has_duplicate_rows(adata):
        print("  Removing duplicate cells...")
        adata = preprocessing.remove_duplicate_cells(adata, keep="first")
        print(f"  ✓ Duplicates removed, {adata.shape[0]} unique cells remain")
    else:
        print("  ✓ No duplicate cells found")

    print(f"  Writing output to: {output_dir}")
    adata.write(output_dir)
    print(f"  ✓ Saved successfully")


if __name__ == "__main__":
    main()