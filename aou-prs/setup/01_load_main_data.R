# Load and harmonize main PRS accuracy results used by multiple figures.

binary <- read_required_tsv(file.path(results_dir, '1223_binary_accuracy_summary_all.tsv'))
quant <- read_required_tsv(file.path(results_dir, '1223_quant_accuracy_summary_all.tsv'))

# subset to PRS-CS only
binary <- binary %>% filter(method == 'PRS-CS')
quant <- quant %>% filter(method == 'PRS-CS')

# read in pheno names mapping file for binary phenos
pheno_names <- read_required_tsv(file.path(results_dir, 'AoU_UKBB_phenotypes_for_r2_plots.txt'))
pheno_names <- pheno_names %>% select(UKB_description, AoU_phenotype)

# include only following quant phenos
quant_in <- c('Alanine_aminotransferase', 'BMI', 'DBP', 'Eosinophil_count', 'Height', 'MCH', 'MCV', 'Neutrophil_count', 'PP_adj',
'RBC_count', 'Reticulocyte_percentage', 'SBP_adj', 'Urea', 'WBC_count')
quant <- quant %>% filter(pheno %in% quant_in)

# map binary names
binary <- left_join(binary, pheno_names, by=c('Phenotype' = 'AoU_phenotype'))

binary <- binary %>% mutate(Group1 = paste(dataset, discovery, sep='-'))
quant <- quant %>% mutate(Group1 = paste(dataset, discovery, sep='-'))

# filter to select phenos for binary
binary_phenos <- c("I25 Chronic ischaemic heart disease", "J44 Other chronic obstructive pulmonary disease", "J45 Asthma",
                "Type 2 diabetes","Disorders of lipoid metabolism", "Coronary atherosclerosis", "Esophagitis, GERD and related diseases",
                "Calculus of kidney")
binary_sel <- binary %>% filter(UKB_description %in% binary_phenos)

quant <- quant %>% mutate(pheno1 = case_when(
                          pheno == 'Alanine_aminotransferase' ~ 'ALT',
                          pheno == 'Eosinophil_count' ~ 'Eosinophil',
                          pheno == 'Neutrophil_count' ~ 'Neutrophil',
                          pheno == 'PP_adj' ~ 'Pulse pressure',
                          pheno == 'RBC_count' ~ 'RBC',
                          pheno == 'Reticulocyte_percentage' ~ 'Reticulocyte',
                          pheno == 'SBP_adj' ~ 'SBP',
                          pheno == 'WBC_count' ~ 'WBC',
                          .default = pheno))

binary_sel <- binary_sel %>% mutate(pheno1 = case_when(
  UKB_description == 'I25 Chronic ischaemic heart disease' ~ 'Ischaemic heart disease',
  UKB_description == 'J44 Other chronic obstructive pulmonary disease' ~ 'COPD',
  UKB_description == 'J45 Asthma' ~ 'Asthma',
  UKB_description == 'Type 2 diabetes' ~ 'T2D',
  UKB_description == "Disorders of lipoid metabolism" ~ 'Lipid metabolism disorders',
  UKB_description == "Esophagitis, GERD and related diseases" ~ 'GERD',
  UKB_description == "Calculus of kidney" ~ 'Kidney stones',
  .default = UKB_description
))