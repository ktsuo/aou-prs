# Shared setup for manuscript figure scripts.

required_packages <- c(
  "tidyverse", "ggplot2", "ggrepel", "data.table", "stringr", "tidytext",
  "gridExtra", "patchwork", "readxl", "ggpubr", "ggh4x", "cowplot", "here"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install missing packages before running this notebook: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(stringr)
  library(tidytext)
  library(gridExtra)
  library(patchwork)
  library(readxl)
  library(ggpubr)
  library(ggh4x)
  library(cowplot)
  library(here)
})

# Project directories used throughout the scripts.
results_dir <- here("AoU_results")
figure_dir <- here("manuscript", "FIGURES")
revision_figure_dir <- here("manuscript", "FIGURES", "FIGURES_revision")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(revision_figure_dir, recursive = TRUE, showWarnings = FALSE)

require_existing_file <- function(path, label = path) {
  if (!file.exists(path)) {
    stop("Missing required file for ", label, ": ", path, call. = FALSE)
  }
  path
}

read_required_tsv <- function(path, ...) {
  data.table::fread(require_existing_file(path), ...)
}

write_public_tsv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", ...)
}
