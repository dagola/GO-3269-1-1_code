
seed <- 20190611
set.seed(seed)

##============================================================================+
## Project setup ----
##============================================================================+
project_name <- "GO-3269-1-1"
work_dirs <- list(
  n01 = list(
    gola = "Documents/"
  ),
  n04 = list(
    gola = "Documents/"
  ),
  gola.local = list(
    gola = "Documents/Gastwissenschaftler"
  ),
  gola.fritz.box = list(
    gola = "Documents/Gastwissenschaftler"
  ),
  damian.novalocal = list(
    damian = "documents"
  )
)

home <- Sys.getenv("HOME")
user <- Sys.info()["user"]
system <- Sys.info()["nodename"]
main_dir <- file.path(home,
                      work_dirs[[system]][[user]],
                      project_name)

##============================================================================+
## Directories ----
##============================================================================+
data_dirs <- list(
  n01 = list(
    gola = file.path("/imbs/projects/", project_name)
  ),
  n04 = list(
    gola = file.path("/imbs/projects", project_name)
  ),
  gola.local = list(
    gola = file.path("/imbs/projects", project_name)
  ),
  gola.fritz.box = list(
    gola = file.path("/imbs/projects", project_name)
  ),
  damian.novalocal = list(
    damian = file.path("/data", project_name)
  )
)
data_dir <- file.path(data_dirs[[system]][[user]])

dir.create(path = input_dir <- file.path(data_dir, "input"), showWarnings = FALSE, recursive = TRUE)
dir.create(path = proc_dir <- file.path(data_dir, "proc"), showWarnings = FALSE, recursive = TRUE)
dir.create(path = output_dir <- file.path(data_dir, "output"), showWarnings = FALSE, recursive = TRUE)
dir.create(path = registries_dir <- file.path(data_dir, "registries"), showWarnings = FALSE, recursive = FALSE)
dir.create(path = scripts_dir <- file.path(main_dir, "scripts"), showWarnings = FALSE, recursive = TRUE)
dir.create(path = functions_dir <- file.path(main_dir, "functions"), showWarnings = FALSE, recursive = TRUE)

##============================================================================+
## Packages ----
##============================================================================+
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(R.utils)
pacman::p_load(mlrMBO)
pacman::p_load(data.table)

##============================================================================+
## Custom functions ----
##============================================================================+
R.utils::sourceDirectory(functions_dir, modifiedOnly = FALSE)

##============================================================================+
## Parameter ----
##============================================================================+
autosomes <- 1:22
nikpay_sample_size <- 187599
additional_pca_samples <- 200
ld_window_size <- 10000 
ld_step_size <- 0.1*ld_window_size
ld_threshold <- 0.2
pca_min_maf <- 0.01
num_pcs <- 20
training_set_size <- 10000

phasing_burn_iterations <- 14
phasing_prune_iterations <- 16
phasing_main_iterations <- 40
phasing_states <- 500
phasing_window_size <- 2
phasing_effective_size <- 11418

impute2_buffer <- 1000
impute2_hap_imputing <- 2500
impute2_width <- 2.5e6
impute2_filter_rules <- c("'TYPE!=Biallelic_SNP'")

##============================================================================+
## Files ----
##============================================================================+
ukb_file_prefix <- switch(
  system,
  damian.novalocal = file.path(input_dir, "UKBB", "ukb_imp"),
  file.path(input_dir, "UKBB", "imputation", "ukb_imp")
)
ukb_sample_file <- switch(
  system,
  damian.novalocal = file.path(input_dir, "UKBB", "ukb48012_imp_chr1_v3_s487317.sample"),
  file.path(input_dir, "UKBB", "imputation", "ukb48012_imp_chr1_v3_s487317.sample")
)
ukb_cal_bed_file_prefix <- switch(
  system,
  damian.novalocal = file.path(input_dir, "UKBB", "ukb_cal"),
  file.path(input_dir, "UKBB", "genotypes", "ukb_cal")
)
ukb_cal_bim_file_prefix <- switch(
  system,
  damian.novalocal = file.path(input_dir, "UKBB", "ukb_snp"),
  file.path(input_dir, "UKBB", "genotypes", "ukb_snp")
)
ukb_cal_fam_file <- switch(
  system,
  damian.novalocal = file.path(input_dir, "UKBB", "ukb48012_cal_chr1_v2_s488285.fam"),
  file.path(input_dir, "UKBB", "genotypes", "ukb48012_cal_chr1_v2_s488285.fam")
)
eb_file_prefix <- file.path(input_dir, "EB", "GSA2")
eb_kinship_file <- file.path(input_dir, "EB", "plink.genome.gz")
de_file_prefix <- file.path(input_dir, "GerMIFS", "german")
nikpay_summary_statistics_file <- file.path(input_dir, "NikpayM_26343387_GCST003116", "harmonised", "26343387-GCST003116-EFO_0000378-build37.f.tsv.gz")
g1k_file_prefix <- switch(
  system,
  damian.novalocal = "/data/1000G/release/20130502/converted/proc/1kg.phase3.v5a.seximputed.nodups.nolongindels",
  file.path("/imbs/external_data/annotation_and_references/hapmap/1000G/1000GP_Phase3/1kg.phase3.v5a.seximputed.nodups.nolongindels")
)

ukb_phenotypes_file <- file.path(input_dir, "phenotypes", "UKB", "ukb35188.tab")
eb_phenotypes_file <- file.path(input_dir, "phenotypes", "EB", "phenot_damian.RData")
g6_file_prefix <- file.path(input_dir, "G6", "G6_1000Gimputed_SnpQcedInfo08")
g6_sample_file <- file.path(input_dir, "G6", "G6_1000Gimputed_SnpQcedInfo08_chr1.sample")
g7_file_prefix <- file.path(input_dir, "G7", "Genotype_QCed_before_imputing", "g7_genotype_qced")

khera_weights_file <- file.path(input_dir, "CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt")

##============================================================================+
## tuning parameters ----
##============================================================================+
tuning_stratify <- TRUE
tuning_cv_iters <- 10
tuning_iters <- 1000
tuning_time_budget <- 2*7*24*60*60 # 4 weeks

##============================================================================+
## Parameter ----
##============================================================================+
autosomes <- 1:22
nikpay_sample_size <- 187599
additional_pca_samples <- 200
ld_window_size <- 10000 
ld_step_size <- 0.1*ld_window_size
ld_threshold <- 0.2
pca_min_maf <- 0.01
num_pcs <- 20
training_set_size <- 10000

phasing_burn_iterations <- 14
phasing_prune_iterations <- 16
phasing_main_iterations <- 40
phasing_states <- 500
phasing_window_size <- 2
phasing_effective_size <- 11418

impute2_buffer <- 1000
impute2_hap_imputing <- 2500
impute2_width <- 2.5e6
impute2_filter_rules <- c("'TYPE!=Biallelic_SNP'")

##============================================================================+
## Files ----
##============================================================================+
# 1000 Genomes Project files for imputation of G7
g1k_imputation_reference_file_prefix <- "/imbs/external_data/annotation_and_references/hapmap/1000G/1000GP_Phase3/1000GP_Phase3"
genetic_map_file_template <- "/imbs/external_data/annotation_and_references/hapmap/1000G/1000GP_Phase3/genetic_map_chr%d_combined_b37.txt"

# Long range LD
long_range_ld_ranges_file <- file.path(input_dir, "long_range_ld.range")
long_range_ld_ranges <- fread(long_range_ld_ranges_file, col.names = c("CHR", "START", "END"), select = 1:3)

# UKB called genotypes file
ukb_called_genotypes_file_prefix <- file.path(proc_dir, "ukb.called")

# Converted files
ukb_converted_file_prefix <- file.path(proc_dir, "ukb")
eb_converted_file_prefix <- file.path(proc_dir, "eb")
g1k_converted_file_prefix <- file.path(proc_dir, "1kg")
de_converted_file_prefix <- file.path(proc_dir, "de")
g6_converted_file_prefix <- file.path(proc_dir, "g6")
g7_converted_file_prefix <- file.path(proc_dir, "g7")

# Common SNP set files
ukb_eb_snps_pos_file <- file.path(proc_dir, "ukb_eb.common.snps_pos")
ukb_eb_de_snps_pos_file <- file.path(proc_dir, "ukb_eb_de.common.snps_pos")
ukb_eb_de_nikpay_snps_pos_file <- file.path(proc_dir, "ukb_eb_de_nikpay.common.snps_pos")
ukb_eb_de_nikpay_g1k_snps_pos_file <- file.path(proc_dir, "ukb_eb_de_nikpay_1kg.common.snps_pos")
ukb_eb_snps_file <- file.path(proc_dir, "ukb_eb.common.snps")
ukb_eb_de_snps_file <- file.path(proc_dir, "ukb_eb_de.common.snps")
ukb_eb_de_nikpay_snps_file <- file.path(proc_dir, "ukb_eb_de_nikpay.common.snps")
ukb_eb_de_nikpay_g1k_snps_file <- file.path(proc_dir, "ukb_eb_de_nikpay_1kg.common.snps")
ukb_eb_plink_range_file <- file.path(proc_dir, "ukb_eb.common.plink_range")
ukb_eb_de_plink_range_file <- file.path(proc_dir, "ukb_eb_de.common.plink_range")
ukb_eb_de_nikpay_plink_range_file <- file.path(proc_dir, "ukb_eb_de_nikpay.common.plink_range")
ukb_eb_de_nikpay_g1k_plink_range_file <- file.path(proc_dir, "ukb_eb_de_nikpay_1kg.common.plink_range")
ukb_eb_bgenix_range_file <- file.path(proc_dir, "ukb_eb.common.bgenix_range")
ukb_eb_de_bgenix_range_file <- file.path(proc_dir, "ukb_eb_de.common.bgenix_range")
ukb_eb_de_nikpay_bgenix_range_file <- file.path(proc_dir, "ukb_eb_de_nikpay.common.bgenix_range")
ukb_eb_de_nikpay_g1k_bgenix_range_file <- file.path(proc_dir, "ukb_eb_de_nikpay_1kg.common.bgenix_range")

# Close relatives files
ukb_remove_closely_relateds_file <- file.path(proc_dir, "ukb_close_relatives_to_remove.samples")
eb_remove_closely_relateds_file <- file.path(proc_dir, "eb_close_relatives_to_remove.samples")
de_remove_closely_relateds_file <- file.path(proc_dir, "de_close_relatives_to_remove.samples")
g6_remove_closely_relateds_file <- file.path(proc_dir, "g6_close_relatives_to_remove.samples")
g7_remove_closely_relateds_file <- file.path(proc_dir, "g7_close_relatives_to_remove.samples")

# Flip files
eb_flip_file <- file.path(proc_dir, "eb_flip.snps")
de_flip_file <- file.path(proc_dir, "de_flip.snps")
g1k_flip_file <- file.path(proc_dir, "1kg_flip.snps")

# Harmonized files
ukb_harmonized_file_prefix <- file.path(proc_dir, "ukb.ukb_harmonized")
eb_harmonized_file_prefix <- file.path(proc_dir, "eb.ukb_harmonized")
g1k_harmonized_file_prefix <- file.path(proc_dir, "1kg.ukb_harmonized")
nikpay_summary_statistics_harmonized_file <- file.path(proc_dir, "nikpay.ukb_harmonized.summary_statistics")
ukb_merge_list_file <- file.path(proc_dir, "ukb_merge.list")
eb_merge_list_file <- file.path(proc_dir, "eb_merge.list")
g1k_merge_list_file <- file.path(proc_dir, "1kg_merge.list")

# 1KG population information
g1k_pop_info_file <- switch(
  system,
  damian.novalocal = "/data/1000G/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  file.path("/imbs/external_data/annotation_and_references/hapmap/1000G/1000GP_Phase3/1000GP_Phase3.sample")
)
pca_EUR_samples_file <- file.path(proc_dir, "pca.EUR.samples")

# Additional PCA sample files
ukb_pca_samples_file <- file.path(proc_dir, "ukb_pca.samples")
eb_pca_samples_file <- file.path(proc_dir, "eb_pca.samples")
de_pca_samples_file <- file.path(proc_dir, "de_pca.samples")

# PCA output files
ukb_pca_file_prefix <- file.path(proc_dir, "ukb_pca")
eb_pca_file_prefix <- file.path(proc_dir, "eb_pca")
de_pca_file_prefix <- file.path(proc_dir, "de_pca")
g6_pca_file_prefix <- file.path(proc_dir, "g6_pca")
g7_pca_file_prefix <- file.path(proc_dir, "g7_pca")

# PCA merge list file
pca_merge_list_file <- file.path(proc_dir, "pca.merge.list")

# PCA Datasets
pca_data_file_prefix <- file.path(proc_dir, "pca_data")

# GRMs
grms_list_file <- file.path(proc_dir, "pca_grms.list")
grms_list_EUR_file <- file.path(proc_dir, "pca_grms.EUR.list")

# PCA plot
pca_plot_file <- file.path(output_dir, "pca.pdf")
pca_EUR_plot_file <- file.path(output_dir, "pca_EUR.pdf")

# Training samples
ukb_training_samples_file <- file.path(proc_dir, "ukb_training.samples")
eb_training_samples_file <- file.path(proc_dir, "eb_training.samples")
de_training_samples_file <- file.path(proc_dir, "de_training.samples")
de_downsampled_training_samples_file <- file.path(proc_dir, "de_downsampled_training.samples")
combined_training_samples_file <- file.path(proc_dir, "combined_training.samples")

# Training datasets
ukb_training_file_prefix <- file.path(proc_dir, "ukb_training")
eb_training_file_prefix <- file.path(proc_dir, "eb_training")
de_training_file_prefix <- file.path(proc_dir, "de_training")
de_downsampled_training_file_prefix <- file.path(proc_dir, "de_downsampled_training")
combined_training_merge_list_file <- file.path(proc_dir, "combined_training.merge_list")
combined_training_file_prefix <- file.path(proc_dir, "combined_training")

# Tasks
ukb_task_file <- file.path(proc_dir, "ukb_task.rds")
ukb_nikpay_task_file <- file.path(proc_dir, "ukb_nikpay_task.rds")
eb_task_file <- file.path(proc_dir, "eb_task.rds")
eb_nikpay_task_file <- file.path(proc_dir, "eb_nikpay_task.rds")
de_task_file <- file.path(proc_dir, "de_task.rds")
de_nikpay_task_file <- file.path(proc_dir, "de_nikpay_task.rds")
de_downsampled_task_file <- file.path(proc_dir, "de_downsampled_task.rds")
de_downsampled_nikpay_task_file <- file.path(proc_dir, "de_downsampled_nikpay_task.rds")
combined_task_file <- file.path(proc_dir, "combined_task.rds")
combined_nikpay_task_file <- file.path(proc_dir, "combined_nikpay_task.rds")

# Models ----
ukb_model_file <- file.path(proc_dir, "ukb_model.rds")
ukb_nikpay_model_file <- file.path(proc_dir, "ukb_nikpay_model.rds")
eb_model_file <- file.path(proc_dir, "eb_model.rds")
eb_nikpay_model_file <- file.path(proc_dir, "eb_nikpay_model.rds")
de_model_file <- file.path(proc_dir, "de_model.rds")
de_nikpay_model_file <- file.path(proc_dir, "de_nikpay_model.rds")
de_downsampled_model_file <- file.path(proc_dir, "de_downsampled_model.rds")
de_downsampled_nikpay_model_file <- file.path(proc_dir, "de_downsampled_nikpay_model.rds")
combined_model_file <- file.path(proc_dir, "combined_model.rds")
combined_nikpay_model_file <- file.path(proc_dir, "combined_nikpay_model.rds")
## GWS SNPs
ukb_gws_model_file <- file.path(proc_dir, "ukb_gws_model.rds")
ukb_gws_nikpay_model_file <- file.path(proc_dir, "ukb_gws_nikpay_model.rds")
eb_gws_model_file <- file.path(proc_dir, "eb_gws_model.rds")
eb_gws_nikpay_model_file <- file.path(proc_dir, "eb_gws_nikpay_model.rds")
de_gws_model_file <- file.path(proc_dir, "de_gws_model.rds")
de_gws_nikpay_model_file <- file.path(proc_dir, "de_gws_nikpay_model.rds")
de_gws_downsampled_model_file <- file.path(proc_dir, "de_gws_downsampled_model.rds")
de_gws_downsampled_nikpay_model_file <- file.path(proc_dir, "de_gws_downsampled_nikpay_model.rds")
combined_gws_model_file <- file.path(proc_dir, "combined_gws_model.rds")
combined_gws_nikpay_model_file <- file.path(proc_dir, "combined_gws_nikpay_model.rds")
## Top 3000 SNPs
ukb_t3000_model_file <- file.path(proc_dir, "ukb_t3000_model.rds")
ukb_t3000_nikpay_model_file <- file.path(proc_dir, "ukb_t3000_nikpay_model.rds")
eb_t3000_model_file <- file.path(proc_dir, "eb_t3000_model.rds")
eb_t3000_nikpay_model_file <- file.path(proc_dir, "eb_t3000_nikpay_model.rds")
de_t3000_model_file <- file.path(proc_dir, "de_t3000_model.rds")
de_t3000_nikpay_model_file <- file.path(proc_dir, "de_t3000_nikpay_model.rds")
de_t3000_downsampled_model_file <- file.path(proc_dir, "de_t3000_downsampled_model.rds")
de_t3000_downsampled_nikpay_model_file <- file.path(proc_dir, "de_t3000_downsampled_nikpay_model.rds")
combined_t3000_model_file <- file.path(proc_dir, "combined_t3000_model.rds")
combined_t3000_nikpay_model_file <- file.path(proc_dir, "combined_t3000_nikpay_model.rds")
## Top 300000 SNPs
ukb_t300000_model_file <- file.path(proc_dir, "ukb_t300000_model.rds")
ukb_t300000_nikpay_model_file <- file.path(proc_dir, "ukb_t300000_nikpay_model.rds")
eb_t300000_model_file <- file.path(proc_dir, "eb_t300000_model.rds")
eb_t300000_nikpay_model_file <- file.path(proc_dir, "eb_t300000_nikpay_model.rds")
de_t300000_model_file <- file.path(proc_dir, "de_t300000_model.rds")
de_t300000_nikpay_model_file <- file.path(proc_dir, "de_t300000_nikpay_model.rds")
de_t300000_downsampled_model_file <- file.path(proc_dir, "de_t300000_downsampled_model.rds")
de_t300000_downsampled_nikpay_model_file <- file.path(proc_dir, "de_t300000_downsampled_nikpay_model.rds")
combined_t300000_model_file <- file.path(proc_dir, "combined_t300000_model.rds")
combined_t300000_nikpay_model_file <- file.path(proc_dir, "combined_t300000_nikpay_model.rds")
## Top 3000000 SNPs
ukb_t3000000_model_file <- file.path(proc_dir, "ukb_t3000000_model.rds")
ukb_t3000000_nikpay_model_file <- file.path(proc_dir, "ukb_t3000000_nikpay_model.rds")
eb_t3000000_model_file <- file.path(proc_dir, "eb_t3000000_model.rds")
eb_t3000000_nikpay_model_file <- file.path(proc_dir, "eb_t3000000_nikpay_model.rds")
de_t3000000_model_file <- file.path(proc_dir, "de_t3000000_model.rds")
de_t3000000_nikpay_model_file <- file.path(proc_dir, "de_t3000000_nikpay_model.rds")
de_t3000000_downsampled_model_file <- file.path(proc_dir, "de_t3000000_downsampled_model.rds")
de_t3000000_downsampled_nikpay_model_file <- file.path(proc_dir, "de_t3000000_downsampled_nikpay_model.rds")
combined_t3000000_model_file <- file.path(proc_dir, "combined_t3000000_model.rds")
combined_t3000000_nikpay_model_file <- file.path(proc_dir, "combined_t3000000_nikpay_model.rds")
## Khera Model
khera_model_file <- file.path(proc_dir, "khera_model.rds")

# Training GWAS
ukb_training_file_prefix <- file.path(proc_dir, "ukb_training")
eb_training_file_prefix <- file.path(proc_dir, "eb_training")
de_training_file_prefix <- file.path(proc_dir, "de_training")

# Testing samples
ukb_testing_samples_file <- file.path(proc_dir, "ukb_testing.samples")
eb_testing_samples_file <- file.path(proc_dir, "eb_testing.samples")
de_testing_samples_file <- file.path(proc_dir, "de_testing.samples")
de_downsampled_testing_samples_file <- file.path(proc_dir, "de_downsampled_testing.samples")
g6_testing_samples_file <- file.path(proc_dir, "g6_testing.samples")
g7_testing_samples_file <- file.path(proc_dir, "g7_testing.samples")

# Testing datasets
ukb_testing_file_prefix <- file.path(proc_dir, "ukb_testing")
eb_testing_file_prefix <- file.path(proc_dir, "eb_testing")
de_testing_file_prefix <- file.path(proc_dir, "de_testing")
de_downsampled_testing_file_prefix <- file.path(proc_dir, "de_downsampled_testing")

# Tasks
ukb_testing_task_file <- file.path(proc_dir, "ukb_testing_task.rds")
eb_testing_task_file <- file.path(proc_dir, "eb_testing_task.rds")
de_testing_task_file <- file.path(proc_dir, "de_testing_task.rds")
de_downsampled_testing_task_file <- file.path(proc_dir, "de_downsampled_testing_task.rds")
g6_testing_task_file <- file.path(proc_dir, "g6_testing_task.rds")
g7_testing_task_file <- file.path(proc_dir, "g7_testing_task.rds")

testing_performance_file <- file.path(proc_dir, "testing_performance.rds")
testing_results_file <- file.path(proc_dir, "testing_results.rds")
