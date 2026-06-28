#!/usr/bin/env Rscript
# Export best-performing PRS metrics.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

quant_best <- quant %>%
  group_by(pheno, pop) %>%
  slice_max(r2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(pheno, pop, r2, r2_2.5, r2_97.5, discovery, dataset, Group1)

binary_best <- binary_sel %>%
  group_by(UKB_description, Pop) %>%
  slice_max(auc2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(UKB_description, Pop, auc2, auc2_2.5, auc2_97.5, discovery, dataset, Group1)

binary_best_aou <- binary_sel %>%
  filter(dataset == "AoU", discovery %in% c("EUR", "AFR", "AMR")) %>%
  group_by(UKB_description, Pop) %>%
  slice_max(auc2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(UKB_description, Pop, h2l_NKr2, auc2, auc2_2.5, auc2_97.5, discovery, dataset, Group1)

write.table(
  quant_best,
  file.path(results_dir, "best_PRS_quant_240122.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  binary_best,
  file.path(results_dir, "best_PRS_binary_240122.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)