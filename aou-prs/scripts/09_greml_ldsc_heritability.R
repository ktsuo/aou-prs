#!/usr/bin/env Rscript
# GREML and LDSC heritability comparison plots.

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

if (!requireNamespace("deming", quietly = TRUE)) {
  stop("Install the 'deming' package before running this script.", call. = FALSE)
}

clean_pheno <- function(x) {
  recode(
    x,
    PP_adj = "PP",
    RBC_count = "RBC",
    WBC_count = "WBC",
    Neutrophil_count = "Neutrophil",
    Eosinophil_count = "Eosinophil",
    Alanine_aminotransferase = "ALT",
    .default = gsub("_", " ", x)
  )
}

read_ukb_hsq <- function(path) {
  fread(path, fill = TRUE)[1:4, ] %>%
    mutate(
      Pheno = sub("\\.hsq$", "", basename(path)),
      Pop = "EUR (UKB)",
      Method = "REML"
    )
}

load_reml_data <- function() {
  reml_aou <- read_required_tsv(file.path(results_dir, "reml_downsample_summary_03202025.tsv")) %>%
    mutate(Pop = paste0(Pop, " (AoU)"))

  ukb_dir <- file.path(results_dir, "UKB_downsampled_GREML", "updated_runs_raw")
  ukb_files <- list.files(ukb_dir, pattern = "\\.hsq$", full.names = TRUE)
  if (length(ukb_files) == 0) {
    stop("No UKB .hsq files found in ", ukb_dir, call. = FALSE)
  }

  pheno_map <- fread(require_existing_file(here("updated_pheno_list_250304.tsv"))) %>%
    select(phenocode, AoU_phenotype)

  reml_ukb_raw <- map_dfr(ukb_files, read_ukb_hsq)
  reml_ukb_mapped <- reml_ukb_raw %>%
    left_join(pheno_map, by = c("Pheno" = "phenocode")) %>%
    mutate(Pheno = AoU_phenotype) %>%
    select(-AoU_phenotype)

  pp_sbp <- reml_ukb_raw %>%
    filter(Pheno %in% c("PP", "SBP")) %>%
    mutate(Pheno = if_else(Pheno == "PP", "PP_adj", Pheno))

  bind_rows(reml_aou, reml_ukb_mapped, pp_sbp) %>%
    filter(Source == "V(G)/Vp", !is.na(Pheno)) %>%
    mutate(Pheno = clean_pheno(Pheno))
}

fit_deming_slope <- function(x, y, se_x, se_y) {
  fit <- deming::deming(y ~ x + 0, data.frame(x = x, y = y), xstd = se_x, ystd = se_y)
  tibble(
    slope = unname(tail(fit$coefficients, 1)),
    ci_lower = unname(tail(fit$ci[, 1], 1)),
    ci_upper = unname(tail(fit$ci[, 2], 1))
  )
}

compare_pops <- function(dat, pop_x, pop_y) {
  wide <- dat %>%
    filter(Pop %in% c(pop_x, pop_y)) %>%
    select(Pop, Pheno, Variance, SE) %>%
    pivot_wider(names_from = Pop, values_from = c(Variance, SE))

  fit_deming_slope(
    wide[[paste0("Variance_", pop_x)]],
    wide[[paste0("Variance_", pop_y)]],
    wide[[paste0("SE_", pop_x)]],
    wide[[paste0("SE_", pop_y)]]
  )
}

plot_h2_scatter <- function(dat, pop_x, pop_y, slope) {
  plot_dat <- dat %>%
    filter(Pop %in% c(pop_x, pop_y)) %>%
    select(Pop, Pheno, Variance, SE) %>%
    pivot_wider(names_from = Pop, values_from = c(Variance, SE)) %>%
    transmute(
      Pheno,
      x = .data[[paste0("Variance_", pop_x)]],
      y = .data[[paste0("Variance_", pop_y)]],
      se_x = .data[[paste0("SE_", pop_x)]],
      se_y = .data[[paste0("SE_", pop_y)]]
    )

  lim <- range(c(plot_dat$x - plot_dat$se_x, plot_dat$x + plot_dat$se_x,
                 plot_dat$y - plot_dat$se_y, plot_dat$y + plot_dat$se_y), na.rm = TRUE)
  txt_x <- lim[2] - diff(lim) * 0.05
  txt_y <- lim[1] + diff(lim) * 0.01

  ggplot(plot_dat, aes(x = x, y = y)) +
    geom_errorbar(aes(ymin = y - se_y, ymax = y + se_y, color = Pheno), alpha = 0.3, linewidth = 0.8) +
    geom_errorbarh(aes(xmin = x - se_x, xmax = x + se_x, color = Pheno), alpha = 0.3) +
    geom_point(aes(color = Pheno), size = 3.5) +
    geom_text_repel(
      aes(label = Pheno),
      size = 7,
      max.overlaps = 20,
      box.padding = 1.7,
      point.padding = 1.3,
      segment.size = 0.5,
      segment.alpha = 0.3,
      color = "grey40",
      force = 10
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_abline(intercept = 0, slope = slope, color = "#8B0000") +
    annotate("text", x = txt_x, y = txt_y, label = sprintf("Slope = %.2f", slope), size = 6, color = "#8B0000") +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(
      x = paste0("GREML heritability estimates from ", pop_x),
      y = paste0("GREML heritability estimates from ", pop_y)
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title = element_text(size = 24),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 20),
      plot.margin = margin(30, 30, 30, 30)
    )
}

plot_ldsc_scatter <- function(dat, slope) {
  lim <- range(c(dat$AoU_h2 - dat$AoU_h2_se, dat$AoU_h2 + dat$AoU_h2_se,
                 dat$UKB_h2 - dat$UKB_h2_se, dat$UKB_h2 + dat$UKB_h2_se), na.rm = TRUE)

  ggplot(dat, aes(x = AoU_h2, y = UKB_h2)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    geom_errorbar(aes(ymin = UKB_h2 - UKB_h2_se, ymax = UKB_h2 + UKB_h2_se), color = "darkgray") +
    geom_errorbarh(aes(xmin = AoU_h2 - AoU_h2_se, xmax = AoU_h2 + AoU_h2_se), color = "darkgray") +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = slope, color = "#8B0000") +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(
      x = "LDSC heritability estimates from EUR (AoU)",
      y = "LDSC heritability estimates from EUR (UKB)"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 24),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 20),
      plot.margin = margin(30, 30, 30, 30)
    )
}

reml_all <- load_reml_data()
write.table(reml_all, file.path(results_dir, "GREML_allres_240409.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

pop_pairs <- list(
  eur_afr = c("EUR (AoU)", "AFR (AoU)"),
  eur_amr = c("EUR (AoU)", "AMR (AoU)"),
  afr_amr = c("AFR (AoU)", "AMR (AoU)"),
  eur_ukb = c("EUR (AoU)", "EUR (UKB)")
)

slopes <- imap_dfr(pop_pairs, function(pair, label) {
  compare_pops(reml_all, pair[1], pair[2]) %>% mutate(comparison = label)
})

ldsc_raw <- read_excel(
  require_existing_file(here("manuscript", "Supplementary_Tables.xlsx")),
  sheet = 6,
  skip = 3,
  col_names = FALSE,
  col_types = "numeric"
)

quant_ldsc <- ldsc_raw[1:14, -c(1, 5, 6, 9, 10)]
names(quant_ldsc) <- c("Pheno", "AoU_h2", "AoU_h2_se", "UKB_h2", "UKB_h2_se")
quant_ldsc$Pheno <- c("BMI", "Height", "ALT", "Urea", "MCV", "MCH", "WBC", "RBC", "Neutrophil", "Eosinophil", "Reticulocyte", "DBP", "SBP", "PP")

ldsc_slope <- fit_deming_slope(quant_ldsc$AoU_h2, quant_ldsc$UKB_h2, quant_ldsc$AoU_h2_se, quant_ldsc$UKB_h2_se)

p_eur_afr <- plot_h2_scatter(reml_all, "EUR (AoU)", "AFR (AoU)", slopes$slope[slopes$comparison == "eur_afr"])
p_eur_amr <- plot_h2_scatter(reml_all, "EUR (AoU)", "AMR (AoU)", slopes$slope[slopes$comparison == "eur_amr"])
p_afr_amr <- plot_h2_scatter(reml_all, "AFR (AoU)", "AMR (AoU)", slopes$slope[slopes$comparison == "afr_amr"])
p_eur_ukb <- plot_h2_scatter(reml_all, "EUR (AoU)", "EUR (UKB)", slopes$slope[slopes$comparison == "eur_ukb"])
p_ldsc <- plot_ldsc_scatter(quant_ldsc, ldsc_slope$slope)

ggsave(file.path(revision_figure_dir, "LDSC_heritability_comparison_EURvsEUR_DEMING.jpg"), p_ldsc, width = 8, height = 8, dpi = 300)
ggsave(file.path(revision_figure_dir, "REML_heritability_comparison_EURvsEUR_DEMING.jpg"), p_eur_ukb, width = 8, height = 8, dpi = 300)
ggsave(file.path(revision_figure_dir, "REML_heritability_comparison_EURvsAFR_DEMING.jpg"), p_eur_afr, width = 8, height = 8, dpi = 300)

combined_reml <- grid.arrange(p_eur_amr, p_afr_amr, ncol = 2, widths = c(1.2, 1.2), heights = 1.2, padding = unit(2, "cm"))
ggsave(file.path(revision_figure_dir, "REML_heritability_comparison_EURvsAFR_AMR_DEMING.jpg"), combined_reml, width = 16, height = 8, dpi = 300)
