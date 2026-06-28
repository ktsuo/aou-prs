#!/usr/bin/env Rscript
# PRS accuracy across target cohorts stratified by EUR ancestry proportion.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

cohort_colors <- c(
  AoU_eur = color_eur,
  AoU_afr = color_afr,
  AoU_amr = color_amr,
  AoU_meta = color_meta
)

cohort_labels <- c(
  AoU_eur = "EUR (AoU)",
  AoU_afr = "AFR (AoU)",
  AoU_amr = "AMR (AoU)",
  AoU_meta = "ALL (AoU)"
)

quant_admix <- read_required_tsv(file.path(results_dir, "aou_gwas_admixed_quant_accuracy.tsv")) %>%
  mutate(
    pop = factor(
      pop,
      levels = paste0("Decile", 1:5),
      labels = paste("Quintile", 1:5)
    ),
    cohort = factor(cohort, levels = names(cohort_colors))
  )

boxplot_admix_quant <- ggplot(quant_admix, aes(x = pop, y = r2, fill = cohort, color = cohort)) +
  geom_boxplot(alpha = 0.3, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  scale_color_manual(values = cohort_colors, labels = cohort_labels, name = "Discovery GWAS") +
  scale_fill_manual(values = cohort_colors, labels = cohort_labels, name = "Discovery GWAS") +
  labs(
    y = bquote("Incremental" ~ italic(R^2)),
    x = "EUR Ancestry Proportion (Target Cohorts)",
    title = "Quantitative"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16, margin = unit(c(5, 0, 0, 0), "mm")),
    axis.title.y = element_text(size = 18, margin = unit(c(0, 3, 0, 0), "mm")),
    axis.text.x = element_text(size = 14),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold", margin = unit(c(0, 0, 5, 0), "mm")),
    panel.spacing = unit(0.1, "in"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(revision_figure_dir, "boxplot_admix_quant.jpg"),
  boxplot_admix_quant,
  width = 8,
  height = 6,
  dpi = 300
)
