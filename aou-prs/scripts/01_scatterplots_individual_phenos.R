#!/usr/bin/env Rscript

# Scatterplots (plot per pheno)

source(file.path("setup", "00_setup.R"))
source(file.path("setup", "01_load_main_data.R"))
source(file.path("setup", "02_plot_helpers.R"))

#############
### QUANT ###
#############

## single ancestry (UKB & AoU) ##
pops_single <- c('EUR (UKB)', 'EUR (AoU)', 'AMR (AoU)', 'AFR (AoU)')
quant_single <- quant %>% filter(pop != 'MID') %>% mutate(Group2 = paste(discovery, " (", dataset, ")", sep=""))
quant_single$pop <- factor(quant_single$pop, levels=c('EUR', 'AFR', 'AMR', 'CSA', 'EAS'))
plot_r2_abs_allphenos(quant_single, pops_single, 'pheno1', pheno1, Group2, 'Group2', r2, pop, r2_2.5, r2_97.5, bquote(italic(R^2)), './manuscript/FIGURES/scatterplot_singleanc_QUANT.jpg', 12, 24, TRUE)

## multi-ancestry ##
pops_multi <- c('ALL (UKB)', 'ALL (AoU)', 'ALL (AoU+UKB)')
quant_multi <- quant %>% filter(pop != 'MID') %>% mutate(Group2 = paste(discovery, " (", dataset, ")", sep=""))
quant_multi$pop <- factor(quant_multi$pop, levels=c('EUR', 'AFR', 'AMR', 'CSA', 'EAS'))
plot_r2_abs_allphenos(quant_multi, pops_multi, 'pheno1', pheno1, Group2, 'Group2', r2, pop, r2_2.5, r2_97.5, bquote(italic(R^2)), './manuscript/FIGURES/scatterplot_multianc_QUANT.jpg', 12, 24, TRUE)

#############
### BINARY ##
#############

## single ancestry (UKB & AoU) ##
binary_single <- binary_sel %>% filter(Pop != 'MID', method == 'PRS-CS') %>% mutate(Group2 = paste(discovery, " (", dataset, ")", sep=""))
binary_single$Pop <- factor(binary_single$Pop, levels=c('EUR', 'AFR', 'AMR', 'CSA', 'EAS'))

plot_r2_abs_allphenos(binary_single, pops_single, 'pheno1', pheno1, Group2, 'Group2', auc2, Pop, auc2_2.5, auc2_97.5, 'AUC', './manuscript/FIGURES/scatterplot_singleanc_BINARY.jpg', 12, 18, FALSE)


## multi-ancestry ##
binary_multi <- binary_sel %>% filter(Pop != 'MID', method == 'PRS-CS') %>% mutate(Group2 = paste(discovery, " (", dataset, ")", sep=""))
binary_multi$Pop <- factor(binary_multi$Pop, levels=c('EUR', 'AFR', 'AMR', 'CSA', 'EAS'))
plot_r2_abs_allphenos(binary_multi, pops_multi, 'pheno1', pheno1, Group2, 'Group2', auc2, Pop, auc2_2.5, auc2_97.5, 'AUC', './manuscript/FIGURES/scatterplot_multianc_BINARY.jpg', 12, 18, FALSE)
