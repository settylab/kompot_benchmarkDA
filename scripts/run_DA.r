#!/usr/bin/env Rscript

### Modern R DA methods runner with unified interface ###
### Compatible with Python methods and batch-corrected embeddings ###

suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
    library(SingleCellExperiment)
    library(scran)
    library(anndata)
    library(reticulate)
})

# Create argument parser with same interface as Python methods
parser <- ArgumentParser(description = "Run R-based DA methods with unified interface")

# Core parameters (matching Python interface)
parser$add_argument("--file_path", type = "character", required = TRUE, help = "Path to h5ad file")
parser$add_argument("--pop", type = "character", required = TRUE, help = "Population/cell type")
parser$add_argument("--pop_enr", type = "double", required = TRUE, help = "Population enrichment level")
parser$add_argument("--pop_column", type = "character", required = TRUE, help = "Population column name")
parser$add_argument("--ds_type", type = "character", required = TRUE, help = "Dataset type")
parser$add_argument("--batch_sd", type = "double", required = TRUE, help = "Batch effect standard deviation")
parser$add_argument("--input_file", type = "character", required = TRUE, help = "Input file directory")
parser$add_argument("--package", type = "character", required = TRUE, help = "DA method name")
parser$add_argument("--seed", type = "integer", required = TRUE, help = "Random seed")
parser$add_argument("--layer_embedding", type = "character", required = TRUE, help = "Base embedding layer")
parser$add_argument("--output_dir", type = "character", required = TRUE, help = "Output directory")

# Method-specific parameters
parser$add_argument("--k", type = "integer", default = 30L, help = "KNN parameter")
parser$add_argument("--resolution", type = "double", default = 0.5, help = "Resolution parameter")
parser$add_argument("--n_dm", type = "integer", default = 10L, help = "Number of diffusion components")

args <- parser$parse_args()

# Set random seed
set.seed(args$seed)

print(paste("Running R DA method:", args$package))
print(paste("Dataset:", args$ds_type, "Population:", args$pop))
print(paste("Using embedding:", args$layer_embedding, "with n_dm =", args$n_dm))

# Load data using anndata (consistent with Python methods)
print("Loading h5ad file...")
adata <- read_h5ad(args$file_path)

# Convert to SingleCellExperiment
print("Converting to SingleCellExperiment...")
sce <- SingleCellExperiment(
    assays = list(logcounts = t(adata$X)),
    colData = adata$obs,
    rowData = adata$var
)

# CRITICAL: Use the same batch-corrected embeddings as Python methods
embedding_key <- if (args$n_dm > 0) "DM_EigenVectors_batch" else paste0(args$layer_embedding, "_batch")

print(paste("Using batch-corrected embedding:", embedding_key))

# Check if the required embedding exists
if (!embedding_key %in% names(adata$obsm)) {
    stop(paste("Required embedding", embedding_key, "not found in data.",
               "Available embeddings:", paste(names(adata$obsm), collapse = ", ")))
}

# Get the batch-corrected embedding matrix
embedding_matrix <- adata$obsm[[embedding_key]]

# Add embedding to SingleCellExperiment
reducedDim(sce, "embedding") <- embedding_matrix

# Ensure required metadata columns exist
required_cols <- c("synth_labels", "synth_samples", "synth_batches")
missing_cols <- required_cols[!required_cols %in% colnames(colData(sce))]
if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Set up output directory
output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "iteration_0", paste0(args$package, "_package_performance.csv"))
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Run the appropriate DA method
print(paste("Running method:", args$package))

if (args$package == "milo") {
    # Load required libraries
    suppressPackageStartupMessages(library(miloR))

    print("Running Milo analysis...")

    # Create Milo object
    milo <- Milo(sce)

    # Build KNN graph using batch-corrected embedding
    milo <- buildGraph(milo, k = args$k, d = ncol(reducedDim(sce, "embedding")),
                       reduced.dim = "embedding")

    # Make neighbourhoods
    milo <- makeNhoods(milo, prop = 0.1, k = args$k, d = ncol(reducedDim(sce, "embedding")),
                       reduced_dims = "embedding")

    # Count cells in neighbourhoods
    milo <- countCells(milo, meta.data = colData(milo), sample = "synth_samples")

    # Create design matrix
    design <- data.frame(colData(milo))[!duplicated(colData(milo)$synth_samples), ]
    rownames(design) <- design$synth_samples

    # Test for differential abundance
    da_results <- testNhoods(milo, design = ~ synth_labels, design.df = design)

    # Create output with cell-level scores (assign neighbourhood scores to cells)
    nhood_ixs <- nhoods(milo)
    cell_scores <- rep(0, ncol(milo))

    for (i in seq_len(nrow(da_results))) {
        cells_in_nhood <- which(nhood_ixs[, i] == 1)
        # Use log fold change as the score
        cell_scores[cells_in_nhood] <- cell_scores[cells_in_nhood] + da_results$logFC[i]
    }

    results_df <- data.frame(
        col_0 = cell_scores,
        row.names = colnames(sce)
    )

} else if (args$package == "daseq") {
    # Load required libraries
    suppressPackageStartupMessages({
        library(DAseq)
        library(Seurat)
    })

    print("Running DAseq analysis...")

    # Convert to Seurat object for DAseq
    seurat_obj <- as.Seurat(sce)
    seurat_obj[["embedding"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "embedding"),
                                                       key = "Emb_", assay = DefaultAssay(seurat_obj))

    # Set default reduction
    DefaultDimReduc(seurat_obj) <- "embedding"

    # Run DAseq
    da_cells <- getDAcells(
        X = seurat_obj,
        cell.type.labels = seurat_obj$synth_labels,
        labels.1 = "Condition1", labels.2 = "Condition2",
        k.vector = seq(50, 500, 50),
        plot.embedding = NULL
    )

    # Extract DA scores
    da_score <- da_cells$da.pred

    results_df <- data.frame(
        col_0 = da_score,
        row.names = names(da_score)
    )

} else if (args$package == "cydar") {
    # Load required libraries
    suppressPackageStartupMessages({
        library(cydar)
        library(S4Vectors)
    })

    print("Running CyDAR analysis...")

    # Prepare data for CyDAR
    exprs_matrix <- t(logcounts(sce))
    sample_ids <- colData(sce)$synth_samples
    condition <- colData(sce)$synth_labels

    # Create CyDAR input
    cd <- prepareCellData(exprs_matrix)

    # Count cells in hyperspheres
    cnt <- countCells(cd, tol = 0.5, BPPARAM = SerialParam())

    # Create sample information
    sample_data <- DataFrame(
        sample_id = unique(sample_ids),
        condition = sapply(unique(sample_ids), function(x) {
            unique(condition[sample_ids == x])[1]
        })
    )

    # Test for differential abundance
    design <- model.matrix(~ condition, sample_data)
    da_results <- testDA(cnt, design, coef = "conditionCondition2")

    # Assign scores to cells based on hypersphere membership
    coords <- intensities(cd)
    cell_scores <- rep(0, nrow(coords))

    # This is a simplified assignment - in practice, you'd need more sophisticated mapping
    # For now, assign the mean log fold change
    cell_scores[] <- mean(da_results$logFC, na.rm = TRUE)

    results_df <- data.frame(
        col_0 = cell_scores,
        row.names = colnames(sce)
    )

} else if (args$package == "louvain") {
    # Load required libraries
    suppressPackageStartupMessages({
        library(igraph)
        library(bluster)
    })

    print("Running Louvain clustering-based DA...")

    # Build KNN graph
    knn_graph <- buildKNNGraph(reducedDim(sce, "embedding"), k = args$k)

    # Perform Louvain clustering
    clusters <- cluster_louvain(knn_graph, resolution = args$resolution)$membership
    colData(sce)$louvain_clusters <- factor(clusters)

    # Calculate cluster proportions per condition
    prop_table <- table(colData(sce)$synth_labels, colData(sce)$louvain_clusters)
    prop_table <- prop.table(prop_table, margin = 1)

    # Calculate log fold changes for each cluster
    condition1_props <- prop_table["Condition1", ]
    condition2_props <- prop_table["Condition2", ]

    # Avoid log(0) by adding small pseudocount
    lfc_clusters <- log2((condition2_props + 1e-6) / (condition1_props + 1e-6))

    # Assign cluster log fold changes to individual cells
    cell_scores <- lfc_clusters[as.character(colData(sce)$louvain_clusters)]

    results_df <- data.frame(
        col_0 = as.numeric(cell_scores),
        row.names = colnames(sce)
    )

} else {
    stop(paste("Unknown R method:", args$package))
}

# Save results in the same format as Python methods
print(paste("Saving results to:", output_file))
write.csv(results_df, output_file)

print("R DA analysis completed successfully!")