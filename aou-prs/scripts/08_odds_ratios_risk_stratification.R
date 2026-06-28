#!/usr/bin/env Rscript
# Odds-ratio risk stratification plots for binary traits.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

phenotype_labels <- c(
  I25 = "Ischaemic heart disease",
  J44 = "COPD",
  J45 = "Asthma",
  X250.2 = "T2D",
  X272 = "Lipid metabolism disorders",
  X411.4 = "Coronary atherosclerosis",
  X530.1 = "GERD",
  X594.1 = "Calculus of kidney"
)

cohort_colors <- c(AoU_afr = color_afr, AoU_amr = color_amr, AoU_eur = color_eur)
cohort_labels <- c(AoU_afr = "AFR", AoU_amr = "AMR", AoU_eur = "EUR")

or_theme <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title.x = element_text(size = 16, margin = unit(c(5, 0, 0, 0), "mm")),
      axis.title.y = element_text(size = 18, margin = unit(c(0, 3, 0, 0), "mm")),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      panel.spacing = unit(0.1, "in"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

plot_or_quantiles <- function(dat, phenotype) {
  ggplot(filter(dat, Phenotype == phenotype), aes(x = Nt, y = OR, color = Cohort, group = Cohort)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = OR_2.5, ymax = OR_97.5), alpha = 0.5, width = 0.2) +
    facet_wrap(~Pop) +
    scale_color_manual(values = cohort_colors, labels = cohort_labels, name = "Discovery GWAS (AoU)") +
    scale_x_continuous(breaks = 1:5) +
    labs(x = "Quantile", y = "Odds Ratio", title = phenotype) +
    or_theme()
}

plot_or_slopes <- function(dat, phenotype) {
  ggplot(filter(dat, Phenotype == phenotype), aes(x = Pop, y = slope, color = Cohort)) +
    geom_errorbar(
      aes(ymin = slope - slope_se, ymax = slope + slope_se),
      position = position_dodge(width = 0.4),
      width = 0.1,
      linewidth = 0.5,
      alpha = 0.5
    ) +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    scale_color_manual(values = cohort_colors, labels = cohort_labels, name = "Discovery GWAS (AoU)") +
    labs(x = "Target Population", y = "OR Slope", title = phenotype) +
    or_theme() +
    theme(plot.margin = margin(30, 30, 30, 30))
}

or_quantiles <- read_required_tsv(file.path(results_dir, "AoU_single_PRS_ORs.02272025.tsv")) %>%
  filter(
    Pop != "MID",
    Phenotype %in% names(phenotype_labels),
    N_tiles == "quantiles",
    Referenced_group == "first_tile",
    !Nt %in% c("topVSmid", "topVSother")
  ) %>%
  mutate(
    Phenotype = recode(Phenotype, !!!phenotype_labels),
    across(c(OR, OR_2.5, OR_97.5, Nt), as.numeric),
    Pop = factor(Pop, levels = c("EUR", "AMR", "AFR", "CSA", "EAS"))
  )

phenotypes <- unique(or_quantiles$Phenotype)
lineplot_or <- wrap_plots(lapply(phenotypes, function(pheno) plot_or_quantiles(or_quantiles, pheno)), ncol = 2) +
  plot_layout(guides = "collect", axes = "collect")

ggsave(file.path(figure_dir, "OR_quantiles.jpg"), lineplot_or, width = 12, height = 15, dpi = 300)

or_slopes <- or_quantiles %>%
  group_by(Phenotype, Pop, Cohort) %>%
  summarize(
    slope = coef(lm(OR ~ Nt))[2],
    slope_se = summary(lm(OR ~ Nt))$coefficients[2, 2],
    .groups = "drop"
  ) %>%
  mutate(Pop = factor(Pop, levels = c("EUR", "AMR", "CSA", "AFR", "EAS")))

slope_plot <- wrap_plots(lapply(phenotypes, function(pheno) plot_or_slopes(or_slopes, pheno)), ncol = 2) +
  plot_layout(guides = "collect")

ggsave(file.path(figure_dir, "OR_slopes.jpg"), slope_plot, width = 12, height = 15, dpi = 300)
