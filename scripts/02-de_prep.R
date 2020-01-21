
# Preparation of DE data

source("init.R")

# Convert DE to PLINK and reduce to common SNPs ----
reg_convert_de <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_de"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = imbs::plink_subset, 
  bfile = de_file_prefix,
  output.prefix = de_converted_file_prefix, 
  more.args = list( 
    extract.custom = sprintf("--extract range %s", ukb_eb_de_plink_range_file),
    exec = plink_exec
  ),
  reg = reg_convert_de
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_de, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-DE", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_de)

# Set variant IDs to CHR:POS format ----
reg_set_ids_de <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_de"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file))
  }, 
  file = sprintf("%s.bim", de_converted_file_prefix),
  reg = reg_set_ids_de
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_de, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-DE", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_de)

# Set family ID to DE ----
reg_set_fid_de <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_de"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$1=\"DE\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file))
  }, 
  file = sprintf("%s.fam", de_converted_file_prefix),
  reg = reg_set_fid_de
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_de, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-DE", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_de)

# Find close relatives ----
reg_kinship_de <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "kinship_de"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_fast_related, 
  input.prefix = de_converted_file_prefix,
  output.prefix = sprintf("%s.kinship", de_converted_file_prefix), 
  more.args = list( 
    plink2.exec = plink2_exec,
    degree = 0.0884, # 2nd degree
    screen.degree = 0.0442, # 3rd degree
    screen.thin = 0.1,
    screen.maf = 0.2
  ),
  reg = reg_kinship_de
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_kinship_de, 
  resources = list(
    ntasks = 1, ncpus = 60, memory = 90000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Kinship-DE", project_name)
  )
)

batchtools::waitForJobs(reg = reg_kinship_de)

# Find maximum independent set of close relatives
de_high_kinship <- data.table::fread(
  file = sprintf("%s.kinship.kin0", de_converted_file_prefix), 
  colClasses = c("character", "character", "character", "character", "integer", "numeric", "numeric", "numeric")
)

reg_mis_de <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "mis_de"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table", "igraph")
)

ids <- batchtools::batchMap(
  fun = find_maximum_independet_set, 
  args = list(kinships = list(de_high_kinship)),
  reg = reg_mis_de
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_mis_de, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_MIS-DE", project_name)
  )
)

batchtools::waitForJobs(reg = reg_mis_de)

samples_to_keep <- data.table::data.table(
  ID = reduceResultsList(reg = reg_mis_de)[[1]]
)
samples_to_remove <- batchtools::ajoin(
  x = de_high_kinship[, .(FID = c(`#FID1`, FID2), ID = c(ID1, ID2))],
  y = samples_to_keep, by = c("ID" = "ID")
)
data.table::fwrite(
  x = unique(samples_to_remove), 
  file = de_remove_closely_relateds_file, 
  quote = FALSE, 
  sep = " ", 
  col.names = FALSE
)
