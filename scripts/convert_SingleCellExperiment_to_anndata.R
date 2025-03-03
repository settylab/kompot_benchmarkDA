#!/usr/bin/env Rscript

# Function to check if required packages are installed, load them, and suppress startup messages
check_and_load_packages <- function(packages) {
  missing_packages <- character()
  for (package in packages) {
    suppressPackageStartupMessages({
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        missing_packages <- c(missing_packages, package)
      }
    })
  }
  if (length(missing_packages) > 0) {
    cat("The following packages are not installed: ", paste(missing_packages, collapse = ", "), "\n")
    cat("Install them by running the following command in R:\n")
    cat("install.packages(c(", paste(sprintf("'%s'", missing_packages), collapse = ", "), "))\n")
    quit(status = 1, save = "no") # Stop execution if any packages are missing
  }
}

# List of required packages
required_packages <- c("SingleCellExperiment", "anndata", "optparse")
check_and_load_packages(required_packages)

# Description and help
description <- paste(
  "This script converts a potentially old Seurat object or a SingleCellExperiment",
  "stored in an RDS file into an AnnData object stored as an H5AD file.",
  "The user can specify input and output file paths, with an option to change",
  "the output filename from .rds to .h5ad if no output is specified."
)

# Set up command-line options
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the input RDS file containing the SingleCellExperiment object. This option is required.",
    metavar = "file"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = paste(
      "Path to the output H5AD file. If not specified,",
      "the output path is derived by replacing the .rds",
      "extension of the input path with .h5ad."
    ),
    metavar = "file"
  ),
  make_option(c("-a", "--assay"),
    type = "character", default = "counts",
    help = "The assay to use as the main matrix (anndata.X). Defaults to 'counts'.",
    metavar = "assay_name"
  )
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list, description = description)
opt <- parse_args(opt_parser)

# Function to print messages with a timestamp
timestamped_cat <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ...)
}

# Check if input file is provided
if (is.null(opt$input)) {
  stop("No input file provided. Use --input to specify the RDS file.", call. = FALSE)
}

# Set output filename
if (is.null(opt$output)) {
  opt$output <- sub("\\.[rR][dD][sS]$", ".h5ad", opt$input, ignore.case = TRUE)
}

timestamped_cat("Loading data from:", opt$input, "\n")
data <- readRDS(opt$input)

# Determine the class of the loaded object and convert if necessary
if ("seurat" %in% tolower(class(data))) {
  timestamped_cat("Summary of input Seurat object:\n\n")
  suppressPackageStartupMessages(print(data))
  cat("\n")
  # Use tryCatch to safely check for Seurat v2 or v3 indicators
  raw_data <- tryCatch(
    {
      if (!is.null(data@raw.data)) data@raw.data else NULL
    },
    error = function(e) NULL
  )

  if (!is.null(raw_data)) {
    timestamped_cat("Old Seurat v2 object detected, attempting to update...\n")
    data <- Seurat::UpdateSeuratObject(data)
  }

  # Convert to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(data)
} else if ("SingleCellExperiment" %in% class(data)) {
  sce <- data
} else {
  suppressPackageStartupMessages({
    sce <- as(data, "SingleCellExperiment")
  })
}
timestamped_cat("Data loaded and converted successfully if needed.\n")

# Print a summary of the input SCE object
timestamped_cat("Summary of input SingleCellExperiment object:\n\n")
print(sce)
cat("\n")


timestamped_cat(sprintf("Processing assay '%s' for anndata.X...\n", opt$assay))
if (!(opt$assay %in% names(assays(sce)))) {
  timestamped_cat(
    sprintf(
      "Error: The specified assay '%s' is not available in the provided SingleCellExperiment object.",
      opt$assay
    ),
    "Use -a or --assay to specify an available assay.\n"
  )
  timestamped_cat("Available assays are: ", paste(names(assays(sce)), collapse = ", "), "\n")
  quit(status = 1, save = "no")
}
X <- t(assay(sce, opt$assay))
timestamped_cat(sprintf("Using '%s' assay as the main data matrix.\n", opt$assay))

timestamped_cat(sprintf("Processing assays other than '%s'...\n", opt$assay))
all_assays <- assays(sce)
all_assays <- all_assays[!names(all_assays) %in% opt$assay]
all_assays <- lapply(all_assays, t)
timestamped_cat("Assays processed.\n")

timestamped_cat("Processing dimensional reductions...\n")
available_reductions <- names(reducedDims(sce))
obsm <- lapply(available_reductions, function(rd_name) {
  reducedDim(sce, rd_name)
})
names(obsm) <- paste("X", tolower(available_reductions), sep = "_")
timestamped_cat("Dimensional reductions processed.\n")

timestamped_cat("Processing and filtering the obs/colData and var/rowData...\n")
obs_data <- as.data.frame(colData(sce, internal = TRUE))
var_data <- as.data.frame(rowData(sce, internal = TRUE))
reduction_prefixes <- paste0(tolower(available_reductions), "\\.")
reduction_columns <- grep("^reducedDims\\.", names(obs_data), value = TRUE)
if (length(reduction_columns) > 0) {
  obs_data <- obs_data[, !(names(obs_data) %in% reduction_columns), drop=FALSE]
}
timestamped_cat("obs/colData processed.\n")

timestamped_cat("Gathering metadata and pairwise matrices...\n")
uns <- c(metadata(sce), int_metadata(sce))
obsp <- list()
varp <- list()
varm <- list()
n_obs <- nrow(X)
n_var <- ncol(X)
for (name in names(uns)) {
  item <- uns[[name]]
  if (is.matrix(item) && all(dim(item) == c(n_obs, n_obs))) {
    obsp[[name]] <- item
    uns[[name]] <- NULL # Remove from uns if added to obsp
    timestamped_cat(sprintf("Transferred '%s' to pairwise observations (obsp).\n", name))
  } else if (is.matrix(item) && nrow(item) == n_obs) {
    obsm[[name]] <- item
    uns[[name]] <- NULL # Remove from uns if added to obsm
    timestamped_cat(sprintf("Transferred '%s' to observation matrices (obsm).\n", name))
  }
  if (is.matrix(item) && all(dim(item) == c(n_var, n_var))) {
    varp[[name]] <- item
    uns[[name]] <- NULL # Remove from uns if added to obsp
    timestamped_cat(sprintf("Transferred '%s' to pairwise variables (varp).\n", name))
  } else if (is.matrix(item) && nrow(item) == n_var) {
    varm[[name]] <- item
    uns[[name]] <- NULL # Remove from uns if added to varm
    timestamped_cat(sprintf("Transferred '%s' to variable matrices (varm).\n", name))
  }
}
timestamped_cat("Metadata and pairwise data organized.\n")


timestamped_cat("Making AnnData object...\n")
ad <- AnnData(
  X = X,
  layers = all_assays,
  obs = obs_data,
  var = var_data,
  obsm = obsm,
  varm = varm,
  obsp = obsp,
  varp = varp,
  uns = uns
)
timestamped_cat("AnnData object created.\n")

# Print a summary of the final AnnData object
timestamped_cat("Conversion complete. Summary of the output AnnData object:\n\n")
print(ad)
cat("\n")

timestamped_cat("Saving the AnnData object to:", opt$output, "\n")
write_h5ad(ad, opt$output)
timestamped_cat("Conversion complete: ", opt$output, "\n")
