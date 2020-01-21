
# Preparation of DE data

source("init.R")

# Convert 1KG to PLINK and reduce to common SNPs ----
reg_convert_g1k <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_1kg"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = imbs::plink_subset, 
  bfile = g1k_file_prefix,
  output.prefix = g1k_converted_file_prefix, 
  more.args = list( 
    extract.custom = sprintf("--extract range %s", ukb_eb_de_nikpay_g1k_plink_range_file),
    exec = plink_exec
  ),
  reg = reg_convert_g1k
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_g1k, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 30000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-1KG", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_g1k)

# Set variant IDs to CHR:POS format ----
reg_set_ids_g1k <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_1kg"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file))
  }, 
  file = sprintf("%s.bim", g1k_converted_file_prefix),
  reg = reg_set_ids_g1k
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_g1k, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-1KG", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_g1k)
