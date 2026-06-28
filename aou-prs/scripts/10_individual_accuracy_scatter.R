#!/usr/bin/env Rscript
# Individual-level PRS accuracy scatterplot

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "02_plot_helpers.R"))

indivacc <- read_required_tsv(file.path(results_dir, "indivAcc_specs_revision.tsv")) %>%
  rename_with(~ "slope_se", .cols = any_of("slope stderr")) %>%
  mutate(
    slope = as.numeric(slope),
    slope_se = as.numeric(slope_se),
    trait = recode(
      trait,
      `White blood cell count` = "WBC",
      `Neutrophil count` = "Neutrophil",
      .default = trait
    )
  )

indivacc_plot <- ggplot(indivacc, aes(x = trait, y = slope, color = train_pop)) +
  geom_point(size = 0.3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = slope - 1.96 * slope_se, ymax = slope + 1.96 * slope_se),
    width = 0.1,
    linewidth = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c(EUR = color_eur, Multi = color_meta), name = "Discovery") +
  scale_y_continuous(
    breaks = c(-2, -1.5, -1, -0.5, 0),
    trans = scales::trans_new(
      name = "uneven",
      transform = function(x) ifelse(x < -1, x / 2 - 0.5, x),
      inverse = function(x) ifelse(x < -1, (x + 0.5) * 2, x)
    )
  ) +
  labs(x = "Trait", y = "Individual-level Accuracy: Slope") +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(revision_figure_dir, "individual_accuracy_scatter.jpg"),
  indivacc_plot,
  width = 6,
  height = 5,
  dpi = 300
)
