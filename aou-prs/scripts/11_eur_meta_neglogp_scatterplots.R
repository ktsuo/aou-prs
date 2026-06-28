#!/usr/bin/env Rscript
# EUR versus AoU meta-analysis -log10(P) scatterplots with gene labels.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(here)
  library(patchwork)
  library(purrr)
})

traits <- c("Height", "BMI", "MCH", "MCV", "Neutrophil_count", "WBC_count")
trait_titles <- c(
  Height = "Height",
  BMI = "BMI",
  MCH = "MCH",
  MCV = "MCV",
  Neutrophil_count = "Neutrophil count",
  WBC_count = "White blood cell count"
)

genome_wide_p <- 5e-8
meta_label_threshold <- 20
summstats_dir <- here("data", "summstats")
annotation_file <- here("data", "annotations", "full_variant_qc_metrics_nearestgenes.txt.gz")
snp_dir <- here("results", "snplists")
figure_dir <- here("manuscript", "FIGURES")

dir.create(snp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

require_existing_file <- function(path, label = path) {
  if (!file.exists(path)) {
    stop("Missing required file for ", label, ": ", path, call. = FALSE)
  }
  path
}

find_one_file <- function(directory, pattern, label) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(files) != 1L) {
    stop(
      "Expected exactly one file for ", label, " matching '", pattern,
      "' in ", directory, "; found ", length(files), ".",
      call. = FALSE
    )
  }
  files[[1]]
}

load_regenie <- function(pop, traits, directory) {
  set_names(traits) %>%
    map(function(trait) {
      path <- find_one_file(directory, paste0("^aou_v7_", pop, "_", trait, "\\.regenie\\.gz$"), paste(pop, trait))
      fread(path, select = c("ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P")) %>%
        mutate(!!paste0(tolower(pop), "_neglog10p") := LOG10P)
    })
}

load_meta <- function(traits, directory) {
  set_names(traits) %>%
    map(function(trait) {
      path <- find_one_file(directory, paste0("^METAANALYSIS_", trait, "_ALLPOPS_AOU_RUN1_1\\.tbl\\.gz$"), paste("meta-analysis", trait))
      fread(path, select = c("MarkerName", "Allele1", "Allele2", "Zscore", "P-value")) %>%
        mutate(meta_neglog10p = -log10(`P-value`))
    })
}

read_annotations <- function(path, variant_ids) {
  annot <- fread(require_existing_file(path, "variant-to-nearest-gene annotation"), header = FALSE)
  if (ncol(annot) < 2L) {
    stop("Annotation file must contain at least two columns: variant ID and nearest gene.", call. = FALSE)
  }
  annot <- annot[, 1:2]
  setnames(annot, c("ID", "nearest_gene"))
  annot %>% filter(ID %in% variant_ids)
}

label_top_variant_per_gene <- function(df) {
  labeled <- df %>%
    filter(!is.na(nearest_gene)) %>%
    group_by(nearest_gene) %>%
    mutate(max_meta = max(meta_neglog10p, na.rm = TRUE)) %>%
    mutate(nearest_gene = if_else(meta_neglog10p == max_meta, nearest_gene, NA_character_)) %>%
    select(-max_meta) %>%
    ungroup()

  bind_rows(labeled, filter(df, is.na(nearest_gene)))
}

plot_neglogp <- function(df, title) {
  ggplot(df, aes(x = eur_neglog10p, y = meta_neglog10p)) +
    geom_point(size = 3, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.8, linetype = "dashed", color = "#6AA5CD") +
    geom_text_repel(
      aes(label = nearest_gene),
      size = 8,
      color = "#941494",
      min.segment.length = 0,
      max.overlaps = 20,
      segment.alpha = 0.7,
      segment.size = 0.2,
      na.rm = TRUE
    ) +
    ggtitle(title) +
    labs(x = expression(-log[10] * P ~ "(AoU EUR)"), y = expression(-log[10] * P ~ "(AoU meta-analysis)")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 32),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 26, color = "black")
    )
}

eur_gwas <- load_regenie("EUR", traits, summstats_dir)
afr_gwas <- load_regenie("AFR", traits, summstats_dir)
meta_gwas <- load_meta(traits, summstats_dir)
eur_meta <- map2(eur_gwas, meta_gwas, inner_join, by = join_by(ID == MarkerName))

meta_sig <- map(eur_meta, filter, meta_neglog10p > meta_label_threshold)
meta_sig_afr <- map2(meta_sig, afr_gwas, left_join, by = "ID")
meta_sig_afr_gws <- map(meta_sig_afr, filter, afr_neglog10p > -log10(genome_wide_p))
variant_ids <- unique(unlist(map(meta_sig_afr_gws, "ID")))
fwrite(data.table(ID = variant_ids), file.path(snp_dir, "meta_afr_gws_alltraits.txt"), col.names = FALSE)

annotations <- read_annotations(annotation_file, variant_ids)
annotated <- map2(eur_meta, meta_sig_afr_gws, function(df, sig_df) {
  sig_annot <- sig_df %>%
    select(ID) %>%
    left_join(annotations, by = "ID") %>%
    mutate(nearest_gene = if_else(is.na(nearest_gene), ID, nearest_gene)) %>%
    distinct(ID, nearest_gene)

  df %>%
    left_join(sig_annot, by = "ID") %>%
    label_top_variant_per_gene()
})

names(annotated) <- unname(trait_titles[names(annotated)])

plots <- imap(annotated, function(df, title) {
  plot <- plot_neglogp(df, title)
  safe_title <- gsub("[^A-Za-z0-9]+", "_", title)
  ggsave(file.path(figure_dir, paste0(safe_title, "_eur_vs_meta_neglogp.jpg")), plot, height = 7, width = 7, dpi = 300)
  plot
})

ggsave(
  file.path(figure_dir, "eur_vs_meta_neglogp_annotated.jpg"),
  wrap_plots(plots),
  height = 8,
  width = 14,
  dpi = 300
)

sessionInfo()
