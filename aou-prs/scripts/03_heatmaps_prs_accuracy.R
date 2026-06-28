#!/usr/bin/env Rscript

# Heatmaps of PRS accuracy across discovery GWAS and target populations.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

alpha <- 0.05
bench <- "UKB-EUR"
targets <- c("EUR", "AFR", "AMR")
disc_keep <- c(
  "AoU+UKB-ALL", "UKB-ALL", "AoU-ALL", "AoU-AFR",
  "AoU-AMR", "AoU-EUR", bench
)

multi_labels <- c(
  "AoU+UKB-ALL" = "AoU+UKB-MULTI",
  "AoU-ALL" = "AoU-MULTI",
  "UKB-ALL" = "UKB-MULTI"
)

rename_multi <- function(x) {
  dplyr::recode(x, !!!multi_labels)
}

# Quantitative traits ---------------------------------------------------------

q_base <- quant %>%
  filter(pop != "MID")

q_test <- compute_pval_diff_quant(q_base, bench) %>%
  mutate(
    pval_label = if_else(
      p_r2var <= alpha & r2 > .data[[paste0(bench, "_r2")]],
      "*",
      NA_character_
    )
  )

q_ref <- q_base %>%
  filter(Group1 == bench) %>%
  select(pop, pheno1, r2, r2_se, Group1) %>%
  mutate(
    `UKB-EUR_r2` = NA_real_,
    `UKB-EUR_se` = NA_real_,
    var_r2 = NA_real_,
    p_r2var = NA_real_,
    pval_label = NA_character_
  )

q_plot <- bind_rows(q_test, q_ref) %>%
  filter(Group1 %in% disc_keep) %>%
  mutate(Group1 = rename_multi(Group1))

plot_heatmap_quant(
  q_plot,
  targets,
  file.path("manuscript", "FIGURES", "heatmap_quant.jpg"),
  height = 6,
  width = 18
)

# Binary traits ---------------------------------------------------------------

b_plot <- binary_sel %>%
  filter(Pop != "MID", Group1 %in% disc_keep, method == "PRS-CS") %>%
  compute_auc_vs_null_pval() %>%
  mutate(
    pval_label = if_else(pval <= alpha, "*", NA_character_),
    Group1 = rename_multi(Group1)
  )

plot_heatmap_binary(
  b_plot,
  targets,
  file.path("manuscript", "FIGURES", "heatmap_binary.jpg"),
  height = 6,
  width = 12
)
