
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

g6_samples <- data.table::fread(
  file = sprintf("%s.fam", g6_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS"),
  colClasses = c("character", "character", "character", "character", "integer", "integer")
)

g7_samples <- data.table::fread(
  file = sprintf("%s.fam", g7_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS"),
  colClasses = c("character", "character", "character", "character", "integer", "integer")
)

# Get close relatives sample list ----
ukb_close_samples <- data.table::fread(
  file = ukb_remove_closely_relateds_file,
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

g6_close_samples <- data.table::fread(
  file = g6_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID"),
  colClasses = c("character", "character")
)

g7_close_samples <- data.table::fread(
  file = g6_remove_closely_relateds_file,
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

# Define UKB testing samples and tasks ----
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
ukb_training_samples <- data.table::fread(
  file = ukb_training_samples_file, 
  col.names = c("FID", "IID", "CAD"),
  colClasses = c("character", "character", "integer")
)
ukb_testing_samples <- ajoin(
  x = ajoin(
    x = ajoin(
      x = ukb_phenotypes[ukb_samples, on = c("ID" = "IID"), nomatch = 0L],
      y = ukb_close_samples, 
      by = c("ID" = "IID")
    ),
    y = pca_EUR_boxed[which(NONEUR)],
    by = c("ID" = "V2")
  ), 
  ukb_training_samples, 
  by = c("ID" = "IID")
)
data.table::fwrite(
  x = ukb_testing_samples[, .(FID, IID = ID, CAD = CAD)], 
  file = ukb_testing_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
ukb_testing_samples <- data.table::fread(
  file = ukb_testing_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
ukb_testing_task <- mlr::makeClassifTask(
  id = "UKB_TESTING", 
  data = as.data.frame(
    ukb_testing_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(ukb_converted_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(ukb_testing_task, file = ukb_testing_task_file)

# Define DE training samples and tasks ----
set.seed(seed)
de_training_samples <- data.table::fread(
  file = de_training_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_testing_samples <- ajoin(
  x = ajoin(
    x = ajoin(
      x = de_samples,
      y = de_close_samples,
      by = c("IID")
    ),
    y = pca_EUR_boxed[which(NONEUR)],
    by = c("IID" = "V2")
  ),
  y = de_training_samples,
  by = c("IID")
)
data.table::fwrite(
  x = de_testing_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = de_testing_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
de_testing_samples <- data.table::fread(
  file = de_testing_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_testing_task <- mlr::makeClassifTask(
  id = "DE_TESTING", 
  data = as.data.frame(
    de_testing_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_converted_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_testing_task, file = de_testing_task_file)

# G6 Testing task ----
data.table::fwrite(
  x = g6_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = g6_testing_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
g6_testing_samples <- data.table::fread(
  file = g6_testing_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
g6_testing_task <- mlr::makeClassifTask(
  id = "DE_TESTING", 
  data = as.data.frame(
    g6_testing_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(g6_converted_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(g6_testing_task, file = g6_testing_task_file)

# G7 Testing task ----
data.table::fwrite(
  x = g7_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = g7_testing_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
g7_testing_samples <- data.table::fread(
  file = g7_testing_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
g7_testing_task <- mlr::makeClassifTask(
  id = "G7_TESTING", 
  data = as.data.frame(
    g7_testing_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(g7_converted_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(g7_testing_task, file = g7_testing_task_file)

# Define DE_DOWNSAMPLED testing samples and tasks ----
set.seed(seed)
de_downsampled_training_samples <- data.table::fread(
  file = de_downsampled_training_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_downsampled_testing_samples <- ajoin(
  x = ajoin(
    x = ajoin(
      x = de_samples,
      y = de_close_samples,
      by = c("IID")
    ),
    y = pca_EUR_boxed[which(NONEUR)],
    by = c("IID" = "V2")
  ),
  y = de_downsampled_training_samples,
  by = c("IID")
)
data.table::fwrite(
  x = de_downsampled_testing_samples[, .(FID, IID, CAD = STATUS - 1)], 
  file = de_downsampled_testing_samples_file, 
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
de_downsampled_testing_samples <- data.table::fread(
  file = de_downsampled_testing_samples_file, 
  col.names = c("FID", "IID", "CAD")
)
de_downsampled_testing_task <- mlr::makeClassifTask(
  id = "DE_DOWNSAMPLED_TESTING", 
  data = as.data.frame(
    de_downsampled_testing_samples[, .(
      FID = factor(FID), 
      IID = factor(IID), 
      STATUS = factor(CAD, levels = c(0, 1), labels = c("Healthy", "CAD")),
      TARGET_FILE_PREFIX = factor(de_converted_file_prefix)
    )]), 
  target = "STATUS",
  positive = "CAD"
)
saveRDS(de_downsampled_testing_task, file = de_downsampled_testing_task_file)
