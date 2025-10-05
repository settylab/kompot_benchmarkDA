#!/usr/bin/env python3
"""
Simple script to convert RDS files to H5AD using rpy2
"""
import argparse
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

def convert_rds_to_h5ad(input_path, output_path):
    """Convert RDS file to H5AD using rpy2"""
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr

        # Activate automatic conversion
        pandas2ri.activate()
        numpy2ri.activate()

        # Import R libraries
        base = importr('base')
        sce_r = importr('SingleCellExperiment')

        print(f"Loading RDS file: {input_path}")

        # Load the RDS file
        r_readRDS = robjects.r['readRDS']
        sce = r_readRDS(str(input_path))

        print(f"RDS file loaded successfully")
        print(f"Object class: {robjects.r['class'](sce)[0]}")

        # Get basic info
        r_nrow = robjects.r['nrow']
        r_ncol = robjects.r['ncol']
        r_colnames = robjects.r['colnames']
        r_rownames = robjects.r['rownames']

        n_genes = int(r_nrow(sce)[0])
        n_cells = int(r_ncol(sce)[0])

        print(f"Dimensions: {n_genes} genes x {n_cells} cells")

        # Get counts matrix
        r_assay = robjects.r['assay']
        counts_matrix = r_assay(sce, 'counts')
        counts_np = np.array(counts_matrix).T  # Transpose to get cells x genes

        # Get cell metadata (colData)
        r_colData = robjects.r['colData']
        cell_metadata = r_colData(sce)
        obs_df = pandas2ri.rpy2py(cell_metadata)

        # Get gene names
        gene_names = list(r_rownames(sce))
        cell_names = list(r_colnames(sce))

        # Create AnnData object
        adata = ad.AnnData(
            X=counts_np,
            obs=obs_df,
            var=pd.DataFrame(index=gene_names)
        )
        adata.obs_names = cell_names

        # Try to get PCA if available
        try:
            r_reducedDim = robjects.r['reducedDim']
            pca_matrix = r_reducedDim(sce, 'PCA')
            adata.obsm['X_pca'] = np.array(pca_matrix)
            print("Added PCA coordinates")
        except:
            print("No PCA coordinates found")

        # Try to get UMAP if available
        try:
            umap_matrix = r_reducedDim(sce, 'UMAP')
            adata.obsm['X_umap'] = np.array(umap_matrix)
            print("Added UMAP coordinates")
        except:
            print("No UMAP coordinates found")

        # Save as H5AD
        print(f"Saving to: {output_path}")
        adata.write_h5ad(output_path)
        print("Conversion completed successfully!")

        return True

    except Exception as e:
        print(f"Error during conversion: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Convert RDS to H5AD')
    parser.add_argument('--input', required=True, help='Input RDS file path')
    parser.add_argument('--output', required=True, help='Output H5AD file path')

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file does not exist: {input_path}")
        return False

    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    success = convert_rds_to_h5ad(input_path, output_path)

    if success:
        print(f"Successfully converted {input_path} to {output_path}")
    else:
        print(f"Failed to convert {input_path}")

    return success

if __name__ == '__main__':
    success = main()
    exit(0 if success else 1)