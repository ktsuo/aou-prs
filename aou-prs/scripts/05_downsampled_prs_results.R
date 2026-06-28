#!/usr/bin/env Rscript
# Downsampled quantitative PRS accuracy plots.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

pheno_labels <- c(
  Alanine_aminotransferase = "ALT",
  Eosinophil_count = "Eosinophil",
  Neutrophil_count = "Neutrophil",
  PP_adj = "Pulse pressure",
  RBC_count = "RBC",
  WBC_count = "WBC"
)

target_pops <- c("EUR", "AFR", "AMR", "CSA", "EAS")
aou_disc <- c("EUR", "AFR", "AMR")
aou_groups <- paste("AoU", aou_disc, sep = "-")
group_levels <- unlist(lapply(target_pops, function(target) paste(aou_groups, target, sep = "-")))
group_sizes <- setNames(
  ifelse(
    sub("^AoU-([A-Z]+)-.*$", "\\1", group_levels) ==
      sub("^AoU-[A-Z]+-([A-Z]+)$", "\\1", group_levels),
    1.3,
    0.3
  ),
  group_levels
)

strip_target <- strip_themed(
  background_x = elem_list_rect(
    fill = alpha(pop_colors[target_pops], 0.1),
    color = pop_colors[target_pops],
    size = c(rep(2, 3), rep(0, 2))
  )
)

boxplot_theme <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 16, margin = unit(c(5, 0, 0, 0), "mm")),
      axis.title.y = element_text(size = 18, margin = unit(c(0, 3, 0, 0), "mm")),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold", margin = unit(c(0, 0, 5, 0), "mm")),
      panel.spacing = unit(0.1, "in"),
      strip.text.x = element_text(size = 14),
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.placement = "outside",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_downsampled_boxplot <- function(dat) {
  ggplot(dat, aes(x = cohort, y = r2)) +
    facet_wrap2(~pop, nrow = 1, strip.position = "bottom", strip = strip_target) +
    geom_boxplot(aes(fill = cohort, color = cohort, size = group_target)) +
    geom_point(aes(color = cohort), size = 1, position = position_jitter(width = 0.1, height = 0)) +
    scale_fill_manual(values = alpha(pop_colors, 0.1), guide = "none") +
    scale_color_manual(values = pop_colors, guide = guide_legend("Discovery GWAS (AoU)")) +
    scale_size_manual(values = group_sizes, guide = "none") +
    labs(y = bquote("Incremental" ~ italic(R^2)), x = "Target Population", title = "Quantitative") +
    boxplot_theme()
}

quant_down <- read_required_tsv(file.path(results_dir, "prscs_accuracy_quant_downsample_0125.tsv")) %>%
  mutate(
    Group1 = paste(ldref, cohort, sep = "-"),
    pheno1 = recode(pheno, !!!pheno_labels, .default = pheno)
  )

quant_aou <- quant_down %>%
  filter(pop != "MID", Group1 %in% aou_groups) %>%
  mutate(
    group_target = factor(paste(Group1, pop, sep = "-"), levels = group_levels),
    pop = factor(pop, levels = target_pops),
    cohort = factor(cohort, levels = aou_disc)
  )

boxplot_quant_aou_down <- make_downsampled_boxplot(quant_aou)
ggsave(
  file.path(revision_figure_dir, "boxplot_singleanc_QUANT_DOWNSAMPLED.jpg"),
  plot = boxplot_quant_aou_down,
  width = 10,
  height = 8,
  dpi = 300
)

quant_single <- quant_down %>%
  filter(pop != "MID") %>%
  mutate(
    Group2 = paste(cohort, " (", ldref, ")", sep = ""),
    pop = factor(pop, levels = target_pops)
  )

plot_r2_abs_allphenos(
  quant_single,
  disc_pops = c("EUR (UKB)", "EUR (AoU)", "AMR (AoU)", "AFR (AoU)"),
  pheno_colname = "pheno1",
  pheno_col = pheno1,
  discovery_col = Group2,
  discovery_colname = "Group2",
  metric = r2,
  target_pop_col = pop,
  metric_ci1 = r2_2.5,
  metric_ci2 = r2_97.5,
  metric_str = bquote(italic(R^2)),
  plotname = file.path(revision_figure_dir, "scatterplot_singleanc_QUANT_DOWNSAMPLED.jpg"),
  height = 12,
  width = 24,
  quant = TRUE
)
