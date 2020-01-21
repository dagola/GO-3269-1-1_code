
source("init.R")

# Get sample list ----
ukb_samples <- data.table::fread(
  file = sprintf("%s.fam", ukb_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS"),
  colClasses = c("character", "character", "character", "character", "integer", "integer")
)

de_samples <- data.table::fread(
  file = sprintf("%s.fam", de_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS")
)

# Get close relatives sample list ----
ukb_close_samples <- data.table::fread(
  file = ukb_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID"),
  colClasses = c("character", "character")
)

eb_close_samples <- data.table::fread(
  file = eb_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID"),
  colClasses = c("character", "character")
)

de_close_samples <- data.table::fread(
  file = de_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID"),
  colClasses = c("character", "character")
)

# Get non EUR by PCA ----
pca_eigenvec <- data.table::fread(sprintf("%s.eigenvec", pca_data_file_prefix))
pca_EUR_eigenvec <- data.table::fread(sprintf("%s.EUR.eigenvec", pca_data_file_prefix))
ukb_proj <- data.table::fread(sprintf("%s.proj.eigenvec", ukb_converted_file_prefix))
ukb_EUR_proj <- data.table::fread(sprintf("%s.EUR.proj.eigenvec", ukb_converted_file_prefix))
de_proj <- data.table::fread(sprintf("%s.proj.eigenvec", de_converted_file_prefix))
de_EUR_proj <- data.table::fread(sprintf("%s.EUR.proj.eigenvec", de_converted_file_prefix))

g1k_pop_info <- data.table::fread(g1k_pop_info_file, header = FALSE, col.names = c("IID", "POP", "GROUP", "SEX"))
ukb_pop_info <- data.table::fread(ukb_pca_samples_file, header = FALSE, col.names = c("FID", "IID"))
de_pop_info <- data.table::fread(de_pca_samples_file, header = FALSE, col.names = c("FID", "IID"))
pop_info <- rbind(
  g1k_pop_info[, .(ID = IID, POP = POP, GROUP = GROUP)], 
  ukb_pop_info[, .(ID = IID, POP = "UKB", GROUP = "EUR")], 
  de_pop_info[, .(ID = IID, POP = "DE", GROUP = "EUR")], 
  fill = TRUE
)

pca_EUR_data <- melt(
  data = pop_info[GROUP == "EUR" & !(POP %in% c("UKB", "DE"))][
    pca_EUR_eigenvec, on = c("ID" = "V2"), nomatch = 0L],
  id.vars = c("V1", "ID", "POP", "GROUP"), 
  variable.name = "PC"
)
pca_EUR_data[, c("MIN", "MAX", "MEAN", "SD") := list(min(value), max(value), mean(value), sd(value)), by = PC]

ukb_EUR_data <- melt(ukb_EUR_proj, id.vars = c("V1", "V2"), variable.name = "PC")
de_EUR_data <- melt(de_EUR_proj, id.vars = c("V1", "V2"), variable.name = "PC")

pca_EUR_boxed <- rbind(de_EUR_data, ukb_EUR_data)[
  unique(pca_EUR_data[, .(PC, MIN, MAX, MEAN, SD)]), 
  nomatch = 0L, 
  on = "PC"]
pca_EUR_boxed[, NONEUR := value < MEAN - 6*SD | value > MEAN + 6*SD]

# Define UKB training samples and tasks ----
set.seed(seed)
ukb_phenotypes_data <- get_UKB_phenotype_data()
ukb_phenotypes_long <- melt(
  ukb_phenotypes_data, 
  id.vars = c("ID", "SEX", "BIRTH_YEAR", "BIRTH_MONTH", "ETHNICITY_0", "ETHNICITY_1", "ETHNICITY_2"),
  measure = patterns(ICD9 = "ICD9", ICD10 = "ICD10", OPCS4 = "OPCS4")
)
rm(ukb_phenotypes_data); gc()
ukb_phenotypes_long[, CAD := 0]
ukb_phenotypes_long[
  grepl("^410|^411|^412|^42979$", ICD9) |
    grepl("^I21|^I22|^I23|^I241$|^I252$", ICD10) | 
    grepl("^K40[1-4]$|^K41[1-4]$|^K45[1–5]$|K49[1-2]$|^K49[8–9]$|K502$|^K75[1–4]$|^K75[8–9]$", OPCS4), 
  CAD := 1
  ]
ukb_phenotypes_long[unique(ukb_phenotypes_long[CAD == 1, .(ID)]), CAD := 1, on = "ID"]
ukb_phenotypes <- unique(
  ukb_phenotypes_long[, .(ID, SEX, BIRTH_YEAR, BIRTH_MONTH, ETHNICITY_0, ETHNICITY_1, ETHNICITY_2, CAD)]
)
rm(ukb_phenotypes_long); gc()
ukb_training_samples <- ajoin(
  x = ajoin(
    x = ukb_phenotypes[!grepl("^\\-", ID)][ukb_samples[!grepl("^\\-", IID)], on = c("ID" = "IID"), nomatch = 0L],
    y = ukb_close_samples, 
    by = c("ID" = "IID")
  ),
  y = pca_EUR_boxed[which(NONEUR)],
  by = c("ID" = "V2")
)[, {
  perc <- training_set_size/.N
  .SD[sample(.N, perc*.N), .(FID, IID = ID), by = CAD]
}]
data.table::fwrite(
  x = ukb_training_samples[, .(FID, IID, CAD = CAD)], 
  file = ukb_training_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
ukb_training_samples <- data.table::fread(
  file = ukb_training_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
ukb_nikpay_task <- mlr::makeClassifTask(
  id = "UKB_NIKPAY", 
  data = as.data.frame(
    ukb_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(ukb_training_file_prefix),
      SUMMARY_STATISTICS_FILE = factor(nikpay_summary_statistics_harmonized_file)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(ukb_nikpay_task, file = ukb_nikpay_task_file)
ukb_task <- mlr::makeClassifTask(
  id = "UKB", 
  data = as.data.frame(
    ukb_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(ukb_training_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(ukb_task, file = ukb_task_file)

# Define DE training samples and tasks ----
set.seed(seed)
de_training_samples <- ajoin(
  x = de_samples,
  y = de_close_samples,
  by = c("IID")
)[, {
  perc <- training_set_size/.N
  .SD[sample(.N, perc*.N), .(FID, IID), by = STATUS]
}]
data.table::fwrite(
  x = de_training_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = de_training_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
de_training_samples <- data.table::fread(
  file = de_training_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_nikpay_task <- mlr::makeClassifTask(
  id = "DE_NIKPAY", 
  data = as.data.frame(
    de_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_training_file_prefix),
      SUMMARY_STATISTICS_FILE = factor(nikpay_summary_statistics_harmonized_file)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_nikpay_task, file = de_nikpay_task_file)
de_task <- mlr::makeClassifTask(
  id = "DE", 
  data = as.data.frame(
    de_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_training_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_task, file = de_task_file)
set.seed(seed)
de_downsampled_training_samples <- ajoin(
  x = de_samples, 
  y = de_close_samples, 
  by = c("IID")
)[, {
  controls <- .SD[STATUS==1]
  cases <- .SD[STATUS==2][sample(.N, 234)] # ~3% cases as in other datasets
  rbind(controls, cases)[, .(FID, IID, STATUS)]
}]
data.table::fwrite(
  x = de_downsampled_training_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = de_downsampled_training_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
de_downsampled_training_samples <- data.table::fread(
  file = de_downsampled_training_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_downsampled_nikpay_task <- mlr::makeClassifTask(
  id = "DE_DOWNSAMPLED_NIKPAY",
  data = as.data.frame(
    de_downsampled_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_downsampled_training_file_prefix),
      SUMMARY_STATISTICS_FILE = factor(nikpay_summary_statistics_harmonized_file)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_downsampled_nikpay_task, file = de_downsampled_nikpay_task_file)
de_downsampled_task <- mlr::makeClassifTask(
  id = "DE_DOWNSAMPLED",
  data = as.data.frame(
    de_downsampled_training_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_downsampled_training_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_downsampled_task, file = de_downsampled_task_file)

# Create training data sets ----
reg_training_data <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "training_data"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "data.table", "checkmate")
)
ids <- batchtools::batchMap(
  fun = plink2_subset,
  input.prefix = c(ukb_converted_file_prefix, eb_converted_file_prefix, de_converted_file_prefix, de_converted_file_prefix),
  output.prefix = c(ukb_training_file_prefix, eb_training_file_prefix, de_training_file_prefix, de_downsampled_training_file_prefix),
  keep = c(ukb_training_samples_file, eb_training_samples_file, de_training_samples_file, de_downsampled_training_samples_file),
  more.args = list(
    plink2.exec = plink2_exec,
    extract = ukb_eb_de_snps_file
  ),
  reg = reg_training_data
)
ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_training_data, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 40000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    account = "GO-3269-1-1",
    name = sprintf("%s_Training_Data", project_name)
  )
)

batchtools::waitForJobs(reg = reg_training_data)
