#!/usr/bin/env Rscript

# Activate the renv environment if available
#if (file.exists("renv/activate.R")) {
#  source("renv/activate.R")
#  print("renv activated")
#} else {
#  warning("renv/activate.R not found. Proceeding without activating renv environment.")
#}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert_seurat.R <input_file.rds> <output_file.h5ad>")
}

input_file <- args[1]
output_file <- args[2]

# Load necessary libraries (adjust based on your package names and functions)
#library(Seurat)
#library(convert2anndata)
#library(anndata)
#library(reticulate)
# Assume that convert_seurat_to_sce and convert_to_anndata are available functions,
# possibly from a custom package or defined elsewhere.
# library(YourConversionPackage)  

# Load the Seurat object from the provided input file path
seurat_obj <- readRDS(input_file)

# Convert to SingleCellExperiment if necessary
sce <- convert_seurat_to_sce(seurat_obj)

# Convert to AnnData. Adjust parameters as needed.
ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)

# Save the AnnData object to the provided output file path
write_h5ad(ad, output_file)
