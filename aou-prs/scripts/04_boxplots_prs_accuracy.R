#!/usr/bin/env Rscript
# Boxplots of PRS accuracy across discovery GWAS groups.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

target_pops <- c("EUR", "AFR", "AMR", "CSA", "EAS")
aou_disc <- c("EUR", "AFR", "AMR")
aou_groups <- paste("AoU", aou_disc, sep = "-")
aou_group_target_levels <- unlist(lapply(target_pops, function(target) paste(aou_groups, target, sep = "-")))
aou_group_sizes <- setNames(
  ifelse(
    sub("^AoU-([A-Z]+)-.*$", "\\1", aou_group_target_levels) ==
      sub("^AoU-[A-Z]+-([A-Z]+)$", "\\1", aou_group_target_levels),
    1.3,
    0.3
  ),
  aou_group_target_levels
)

strip_target <- strip_themed(
  background_x = elem_list_rect(
    fill = alpha(pop_colors[target_pops], 0.1),
    color = pop_colors[target_pops],
    size = c(rep(2, 3), rep(0, 2))
  )
)

strip_target_multi <- strip_themed(
  background_x = elem_list_rect(
    fill = alpha(pop_colors[c("EUR", "AFR", "AMR")], 0.3),
    color = rep(0, 3)
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

prep_single_aou <- function(dat, target_col, metric_col) {
  dat %>%
    transmute(
      target = .data[[target_col]],
      discovery,
      Group1,
      group_target = paste(Group1, .data[[target_col]], sep = "-"),
      metric = .data[[metric_col]],
      pheno1
    ) %>%
    filter(target != "MID", Group1 %in% aou_groups) %>%
    mutate(
      target = factor(target, levels = target_pops),
      discovery = factor(discovery, levels = aou_disc),
      group_target = factor(group_target, levels = aou_group_target_levels)
    )
}

plot_single_aou_boxplot <- function(dat, y_lab, title, null_line = NULL) {
  p <- ggplot(dat, aes(x = discovery, y = metric)) +
    facet_wrap2(~target, nrow = 1, strip.position = "bottom", strip = strip_target) +
    geom_boxplot(aes(fill = discovery, color = discovery, size = group_target)) +
    geom_point(aes(color = discovery), size = 1, position = position_jitter(width = 0.1, height = 0)) +
    scale_fill_manual(values = alpha(pop_colors, 0.1), guide = "none") +
    scale_color_manual(values = pop_colors, guide = guide_legend("Discovery GWAS (AoU)")) +
    scale_size_manual(values = aou_group_sizes, guide = "none") +
    labs(y = y_lab, x = "Target Population", title = title) +
    boxplot_theme()

  if (!is.null(null_line)) {
    p <- p + geom_hline(yintercept = null_line, linetype = "dashed", linewidth = 0.7, color = "darkgrey")
  }
  p
}

plot_multi_boxplot <- function(dat) {
  dat <- dat %>%
    filter(!pop %in% c("MID", "CSA", "EAS"), Group1 %in% c("UKB-ALL", "AoU-ALL", "AoU+UKB-ALL")) %>%
    mutate(
      pop = factor(pop, levels = c("EUR", "AFR", "AMR")),
      dataset = factor(dataset, levels = c("UKB", "AoU", "AoU+UKB"))
    )

  ggplot(dat, aes(x = dataset, y = r2)) +
    facet_wrap2(~pop, nrow = 1, strip.position = "bottom", strip = strip_target_multi) +
    geom_boxplot(aes(fill = dataset, color = dataset)) +
    geom_point(aes(color = dataset), size = 1, position = position_jitter(width = 0.1, height = 0)) +
    scale_fill_manual(values = alpha(pop_colors, 0.1), guide = "none") +
    scale_color_manual(values = pop_colors, guide = guide_legend("Discovery GWAS (Multi-ancestry)")) +
    labs(y = bquote("Incremental" ~ italic(R^2)), x = "Target Population", title = "Quantitative") +
    boxplot_theme()
}

save_figure <- function(plot, filename, width = 8, height = 6) {
  ggsave(file.path(figure_dir, filename), plot = plot, width = width, height = height, dpi = 300)
}

quant_aou <- prep_single_aou(quant, target_col = "pop", metric_col = "r2")
binary_aou <- prep_single_aou(binary_sel, target_col = "Pop", metric_col = "auc2")
boxplot_quant_aou <- plot_single_aou_boxplot(quant_aou, bquote("Incremental" ~ italic(R^2)), "Quantitative")
boxplot_binary_aou <- plot_single_aou_boxplot(binary_aou, "AUC", "Binary", null_line = 0.5)
boxplot_quant_meta <- plot_multi_boxplot(quant)

save_figure(boxplot_quant_aou, "boxplot_singleanc_aou_QUANT.jpg")
save_figure(boxplot_binary_aou, "boxplot_singleanc_aou_BINARY.jpg")
save_figure(boxplot_quant_meta, "boxplot_multianc_ukb_aou_QUANT.jpg")
