#!/usr/bin/env Rscript
# Relative quantitative PRS accuracy for P+T versus PRS-CS.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

quant_phenos <- c(
  "Alanine_aminotransferase", "BMI", "DBP", "Eosinophil_count", "Height", "MCH", "MCV",
  "Neutrophil_count", "PP_adj", "RBC_count", "Reticulocyte_percentage", "SBP_adj", "Urea", "WBC_count"
)

pheno_labels <- c(
  Alanine_aminotransferase = "ALT",
  Eosinophil_count = "Eosinophil",
  Neutrophil_count = "Neutrophil",
  PP_adj = "Pulse pressure",
  RBC_count = "RBC",
  Reticulocyte_percentage = "Reticulocyte",
  SBP_adj = "SBP",
  WBC_count = "WBC"
)

load_prscs <- function() {
  read_required_tsv(file.path(results_dir, "1223_quant_accuracy_summary_all.tsv")) %>%
    filter(pheno %in% quant_phenos, method == "PRS-CS", cohort %in% c("AoU_afr", "AoU_eur", "AoU_amr")) %>%
    mutate(
      pheno1 = recode(pheno, !!!pheno_labels, .default = pheno),
      cohort = recode(cohort, AoU_afr = "AoU_AFR", AoU_eur = "AoU_EUR", AoU_amr = "AoU_AMR")
    ) %>%
    select(-any_of(c("pheno", "prefix", "discovery", "dataset", "Group")))
}

load_pt <- function() {
  read_required_tsv(file.path(results_dir, "1223_aouGWAS_pt_accuracy_summary_quant.tsv")) %>%
    filter(pheno %in% quant_phenos, grepl("S1$", prefix)) %>%
    mutate(
      method = "P+T",
      pheno1 = recode(pheno, !!!pheno_labels, .default = pheno)
    ) %>%
    select(-any_of(c("prefix", "pheno")))
}

plot_relative_accuracy <- function(dat) {
  ref <- dat %>%
    filter(cohort == "AoU_EUR", pop == "EUR") %>%
    select(pheno1, method, ref_r2 = r2)

  plot_dat <- dat %>%
    filter(cohort == "AoU_EUR") %>%
    left_join(ref, by = c("pheno1", "method")) %>%
    mutate(
      rel_acc = r2 / ref_r2,
      pop = factor(pop, levels = c("EUR", "AMR", "CSA", "EAS", "AFR")),
      method = factor(method, levels = c("P+T", "PRS-CS"))
    ) %>%
    filter(pop != "MID")

  decay <- plot_dat %>%
    group_by(pop, method) %>%
    summarize(
      y = max(rel_acc, na.rm = TRUE) + 0.1,
      label = if_else(pop == "EUR", "", sprintf("-%d%%", round(mean(1 - rel_acc, na.rm = TRUE) * 100))),
      .groups = "drop"
    )

  outliers <- plot_dat %>%
    filter(pop != "EUR") %>%
    group_by(pop, method) %>%
    mutate(
      q1 = quantile(rel_acc, 0.25, na.rm = TRUE),
      q3 = quantile(rel_acc, 0.75, na.rm = TRUE),
      is_outlier = rel_acc > q3 + 1.5 * (q3 - q1),
      x_dodge = as.numeric(pop) + if_else(method == "PRS-CS", 0.2, -0.2)
    ) %>%
    ungroup() %>%
    filter(is_outlier)

  ggplot(plot_dat, aes(x = pop, y = rel_acc, fill = method)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.2, alpha = 0.7) +
    geom_text(data = decay, aes(y = y, label = label, color = method), position = position_dodge(width = 0.8), size = 6) +
    geom_point(data = outliers, aes(x = x_dodge, color = method), size = 1) +
    geom_text_repel(
      data = outliers,
      aes(x = x_dodge, label = pheno1),
      size = 4,
      box.padding = 0.9,
      point.padding = 0.3,
      force = 10,
      direction = "both",
      segment.size = 0.1,
      min.segment.length = 0.05,
      max.overlaps = 10,
      nudge_x = 0.2,
      nudge_y = -0.03
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("P+T" = "#228833", "PRS-CS" = "#443399")) +
    scale_color_manual(values = c("P+T" = "#228833", "PRS-CS" = "#443399"), guide = "none") +
    labs(x = "Target Population", y = "Relative Accuracy", fill = "Method") +
    theme_classic2() +
    theme(
      text = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
}

quant_combined <- bind_rows(load_prscs(), load_pt())

rel_acc <- plot_relative_accuracy(quant_combined)
ggsave(file.path(figure_dir, "violin_relative_accuracy_AoU_EUR.jpg"), rel_acc, width = 13, height = 8, dpi = 300)

rel_acc_no_retic <- quant_combined %>%
  filter(pheno1 != "Reticulocyte") %>%
  plot_relative_accuracy()

ggsave(
  file.path(revision_figure_dir, "violin_relative_accuracy_AoU_EUR_no_retic.jpg"),
  rel_acc_no_retic,
  width = 13,
  height = 8,
  dpi = 300
)
