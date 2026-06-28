# Shared plotting and statistical helper functions.
# Assumes setup/00_setup.R has already been sourced.

color_amr <- "#ED1E24"
color_eur <- "#6AA5CD"
color_eur1 <- "#6AC5CD"
color_afr <- "#941494"
color_csa <- "#FF9912"
color_eas <- "#108C44"
color_mid <- "#33CC33"
color_meta <- "#33CC33"

pop_colors <- c(
  "AoU-AFR" = color_afr,
  "AFR (AoU)" = color_afr,
  "AoU-AMR" = color_amr,
  "AMR (AoU)" = color_amr,
  "AoU-EUR" = color_eur,
  "EUR (AoU)" = color_eur,
  "UKB-EUR" = color_eur1,
  "AoU-ALL" = color_meta,
  "AoU+UKB-ALL" = color_csa,
  "EUR" = color_eur,
  "AFR" = color_afr,
  "AMR" = color_amr,
  "CSA" = color_csa,
  "EAS" = color_eas,
  "MID" = color_mid,
  "META" = color_csa,
  "EUR (UKB)" = color_eur1,
  "ALL (UKB)" = "#E4D00A",
  "ALL (AoU)" = "#33CC33",
  "ALL (AoU+UKB)" = "brown",
  "AoU+UKB" = "brown",
  "AoU" = "#33CC33",
  "UKB" = "#E4D00A",
  "A" = "brown",
  "B" = "coral"
)

single_ancestry_strip <- function() {
  pops <- c("EUR", "AFR", "AMR", "CSA", "EAS")
  strip_themed(
    background_x = elem_list_rect(
      fill = alpha(pop_colors[pops], 0.1),
      color = pop_colors[pops],
      size = c(rep(0.8, 3), 0, 0)
    )
  )
}

plot_accuracy_phenotype <- function(dat, discovery, metric, target_pop, ci_min, ci_max,
                                    title, y_label, null_line = NULL,
                                    point_size = 2, base_size = 18,
                                    show_legend = FALSE) {
  p <- ggplot(dat, aes(x = {{ discovery }}, y = {{ metric }}, color = {{ discovery }})) +
    geom_point(size = point_size, position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = {{ ci_min }}, ymax = {{ ci_max }}),
      width = 0,
      position = position_dodge(width = 0.5)
    ) +
    facet_grid2(cols = vars({{ target_pop }}), scales = "free", switch = "x", strip = single_ancestry_strip()) +
    scale_x_discrete(expand = c(0.9, 0)) +
    scale_color_manual(values = pop_colors) +
    labs(title = title, y = y_label, x = "Target Population") +
    guides(color = guide_legend("Discovery GWAS")) +
    theme_classic() +
    theme(
      text = element_text(size = base_size),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.line.x.bottom = element_blank(),
      strip.text.x = element_text(face = "bold"),
      strip.placement = "outside",
      plot.margin = margin(l = 20, t = 5, b = 20),
      legend.position = if (show_legend) "right" else "none"
    )

  if (!is.null(null_line)) {
    p <- p + geom_hline(yintercept = null_line, linetype = "dashed", linewidth = 0.5, color = "darkgrey")
  }

  p
}

plot_r2_abs_allphenos <- function(input_data, disc_pops, pheno_colname, pheno_col,
                                  discovery_col, discovery_colname, metric,
                                  target_pop_col, metric_ci1, metric_ci2,
                                  metric_str, plotname, height, width, quant) {
  dat <- input_data %>% filter({{ discovery_col }} %in% disc_pops)
  dat[[discovery_colname]] <- factor(dat[[discovery_colname]], levels = disc_pops)

  phenos <- if (quant) {
    c(
      "BMI", "Height", "ALT", "Urea", "MCH", "MCV", "RBC", "WBC",
      "Neutrophil", "Eosinophil", "DBP", "SBP", "Pulse pressure"
    )
  } else {
    unique(dat[[pheno_colname]])
  }

  plot_one <- function(phenotype_name) {
    pheno_dat <- dat %>% filter({{ pheno_col }} == phenotype_name)
    plot_accuracy_phenotype(
      dat = pheno_dat,
      discovery = {{ discovery_col }},
      metric = {{ metric }},
      target_pop = {{ target_pop_col }},
      ci_min = {{ metric_ci1 }},
      ci_max = {{ metric_ci2 }},
      title = phenotype_name,
      y_label = metric_str,
      null_line = if (quant) NULL else 0.5,
      point_size = if (quant) 2 else 3,
      base_size = if (quant) 18 else 16,
      show_legend = !quant
    )
  }

  plots <- lapply(phenos, plot_one)

  if (quant) {
    combined_plot <- assemble_quant_accuracy_grid(plots)
  } else {
    combined_plot <- wrap_plots(plots) + plot_layout(guides = "collect", axis_titles = "collect")
  }

  ggsave(plotname, combined_plot, height = height, width = width)
  invisible(combined_plot)
}

assemble_quant_accuracy_grid <- function(plots) {
  group_label <- function(label) {
    ggdraw() + draw_label(label, fontface = "bold", angle = 90, x = 0.2, y = 0.5, size = 20)
  }

  anthro <- cowplot::plot_grid(group_label("Anthro"), plots[[1]], plots[[2]], nrow = 1, rel_widths = c(0.05, 1, 1))
  biomarkers <- cowplot::plot_grid(group_label("Biomarkers"), plots[[3]], plots[[4]], nrow = 1, rel_widths = c(0.05, 1, 1))
  blood_panel <- cowplot::plot_grid(
    group_label("Blood Panel"), plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]],
    nrow = 1,
    rel_widths = c(0.05, rep(1, 7))
  )
  blood_pressure <- cowplot::plot_grid(
    group_label("Blood Pressure"), plots[[12]], plots[[13]],
    nrow = 1,
    rel_widths = c(0.05, 1, 1)
  )

  anthro + biomarkers + blood_panel + blood_pressure +
    plot_layout(design = "
      AA#####
      BB#####
      CCCCCCC
      DD#####
    ")
}

plot_heatmap_quant <- function(input_data, target_pops, plotname, height, width) {
  plots <- lapply(target_pops, function(target_pop) {
    dat <- input_data %>%
      filter(.data$pop == .env$target_pop) %>%
      select(pheno1, value = r2, Group1, pval_label) %>%
      add_best_metric_border("value") %>%
      mutate(
        Group1 = factor(Group1, levels = discovery_order()),
        pheno_group = phenotype_group(pheno1),
        pheno1 = factor(pheno1, levels = quant_pheno_order())
      )

    heatmap_plot(dat, target_pop, fill_col = "value", fill_label = bquote("Incremental" ~ italic(R^2))) +
      facet_grid(~pheno_group, scales = "free", space = "free_x", switch = "x") +
      theme(
        strip.background.y = element_rect(fill = "white", color = "gray80"),
        strip.text = element_text(size = 7.5, face = "bold"),
        axis.title.y = element_text(size = 18)
      )
  })

  combined_plot <- combine_three_target_heatmaps(plots)
  ggsave(plotname, combined_plot, height = height, width = width)
  invisible(combined_plot)
}

plot_heatmap_binary <- function(input_data, target_pops, plotname, height, width) {
  plots <- lapply(target_pops, function(target_pop) {
    dat <- input_data %>%
      filter(.data$Pop == .env$target_pop) %>%
      select(pheno1, value = auc2, Group1, pval_label) %>%
      add_best_metric_border("value") %>%
      mutate(Group1 = factor(Group1, levels = discovery_order()))

    heatmap_plot(dat, target_pop, fill_col = "value", fill_label = "AUC")
  })

  combined_plot <- combine_three_target_heatmaps(plots)
  ggsave(plotname, combined_plot, height = height, width = width)
  invisible(combined_plot)
}

discovery_order <- function() {
  c("AoU+UKB-MULTI", "UKB-MULTI", "UKB-EUR", "AoU-MULTI", "AoU-EUR", "AoU-AFR", "AoU-AMR")
}

quant_pheno_order <- function() {
  c(
    "BMI", "Height", "ALT", "Urea", "MCH", "MCV", "RBC", "WBC",
    "Neutrophil", "Eosinophil", "Reticulocyte", "DBP", "SBP", "Pulse pressure"
  )
}

phenotype_group <- function(pheno) {
  case_when(
    pheno %in% c("BMI", "Height") ~ "Anthro",
    pheno %in% c("ALT", "Urea") ~ "Biomarkers",
    pheno %in% c("DBP", "SBP", "Pulse pressure") ~ "Blood Pressure",
    pheno %in% c("Eosinophil", "MCH", "MCV", "Neutrophil", "RBC", "Reticulocyte", "WBC") ~ "Blood Panel",
    TRUE ~ "Other"
  )
}

add_best_metric_border <- function(dat, metric_col) {
  dat %>%
    group_by(pheno1) %>%
    mutate(border = .data[[metric_col]] == max(.data[[metric_col]])) %>%
    ungroup()
}

heatmap_plot <- function(dat, pop, fill_col, fill_label) {
  ggplot(dat, aes(x = pheno1, y = Group1, fill = .data[[fill_col]])) +
    geom_tile() +
    geom_tile(data = filter(dat, border), aes(x = pheno1, y = Group1), fill = NA, color = "black", linewidth = 0.5) +
    geom_text(aes(label = pval_label), size = 7) +
    scale_fill_gradient(
      limits = c(min(dat[[fill_col]]), max(dat[[fill_col]])),
      high = pop_colors[[pop]],
      low = "white",
      name = fill_label
    ) +
    labs(y = "Discovery GWAS", title = paste("Target Population:", pop)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1.5, "line"),
      legend.text = element_text(size = 9),
      panel.spacing = unit(0.1, "in"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

combine_three_target_heatmaps <- function(plots) {
  cowplot::plot_grid(
    plots[[1]],
    NULL,
    plots[[2]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()),
    NULL,
    plots[[3]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()),
    nrow = 1,
    rel_widths = c(1, 0.1, 0.8, 0.1, 0.8)
  )
}

compute_pval_diff_quant <- function(input_data, benchmark) {
  benchmark_r2_col <- paste0(benchmark, "_r2")
  benchmark_se_col <- paste0(benchmark, "_se")

  benchmark_data <- input_data %>%
    filter(Group1 == benchmark) %>%
    select(pheno1, pop, !!benchmark_r2_col := r2, !!benchmark_se_col := r2_se)

  input_data %>%
    filter(Group1 != benchmark) %>%
    select(pop, pheno1, r2, r2_se, Group1) %>%
    left_join(benchmark_data, by = c("pop", "pheno1")) %>%
    mutate(
      var_r2 = ((r2 - .data[[benchmark_r2_col]])^2) / (r2_se^2 + .data[[benchmark_se_col]]^2),
      p_r2var = pchisq(var_r2, df = 1, lower.tail = FALSE)
    )
}

compute_auc_vs_null_pval <- function(input_data) {
  input_data %>%
    select(Pop, pheno1, auc2, auc2_97.5, auc2_2.5, Group1) %>%
    mutate(
      SE = (auc2_97.5 - auc2_2.5) / (2 * 1.96),
      t_stat = (auc2 - 0.5) / SE,
      pval = 2 * pt(-abs(t_stat), df = 1)
    )
}
