
# Preparation of UKB data

source("init.R")

# Convert UKB to PLINK ----
reg_convert_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_bgen2bed_conversion, 
  bgen.file = sprintf("%s_chr%d_v3.bgen", ukb_file_prefix, autosomes),
  sample.file = ukb_sample_file,
  output.prefix = sprintf("%s.chr%02d", ukb_converted_file_prefix, autosomes), 
  more.args = list( 
    extract = sprintf("--extract range %s", ukb_eb_plink_range_file),
    hardcall.threshold = 0.49999,
    plink2.exec = plink2_exec
  ),
  reg = reg_convert_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_ukb, 
  resources = list(
    ntasks = 1, ncpus = 60, memory = 90000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_ukb)

# Merge chromosomes to one single file ----
reg_concat_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "concat_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = function(bed.files, bim.files, fam.files, output.prefix) {
    
    # Copy fam file
    cat("Copying fam file...")
    imbs::system_call(bin = "cat", c(fam.files[1], ">", sprintf("%s.fam", output.prefix)))
    
    # Combine bim files
    cat("Merging bim files...")
    imbs::system_call(bin = "cat", c(paste(bim.files, collapse = " "), ">", sprintf("%s.bim", output.prefix)))
    
    # Combine bed files
    cat("Merging bed files...")
    imbs::system_call(bin = "bash", sprintf("-c '(echo -en \"\\x6C\\x1B\\x01\"; for g in %s; do tail -qc +4 $g; done) > %s.bed'", paste(bed.files, collapse = " "), output.prefix))
    
    cat("Done.")
    
  }, 
  output.prefix = sprintf("%s", ukb_converted_file_prefix), 
  more.args = list( 
    bed.files = sprintf("%s.chr%02d.bed", ukb_converted_file_prefix, autosomes),
    bim.files = sprintf("%s.chr%02d.bim", ukb_converted_file_prefix, autosomes),
    fam.files = sprintf("%s.chr%02d.fam", ukb_converted_file_prefix, autosomes)
  ),
  reg = reg_concat_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_concat_ukb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Concat-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_concat_ukb)

# Set variant IDs to CHR:POS format ----
reg_set_ids_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file))
  }, 
  file = sprintf("%s.bim", ukb_converted_file_prefix),
  reg = reg_set_ids_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_ukb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_ukb)

# Set family ID to UKB in merged imputed file ----
reg_set_fid_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$1=\"UKB\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file))
  }, 
  file = sprintf("%s.fam", ukb_converted_file_prefix),
  reg = reg_set_fid_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_ukb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_ukb)

# Combine UKB genotype calls ----
reg_merge_cal_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "merge_cal_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = function(bed.files, bim.files, fam.files, output.prefix) {
    
    # Copy fam file
    cat("Copying fam file...")
    imbs::system_call(bin = "cat", c(fam.files[1], ">", sprintf("%s.fam", output.prefix)))
    
    # Combine bim files
    cat("Merging bim files...")
    imbs::system_call(bin = "zcat", c(paste(bim.files, collapse = " "), ">", sprintf("%s.bim", output.prefix)))
    
    # Combine bed files
    cat("Merging bed files...")
    imbs::system_call(bin = "bash", sprintf("-c '(echo -en \"\\x6C\\x1B\\x01\"; for g in %s; do gzip -cd $g | tail -qc +4; done) > %s.bed'", paste(bed.files, collapse = " "), output.prefix))
    
    cat("Done.")
    
  }, 
  output.prefix = ukb_called_genotypes_file_prefix, 
  more.args = list(
    bed.files = sprintf("%s_chr%d_v2.bed.gz", ukb_cal_bed_file_prefix, autosomes),
    bim.files = sprintf("%s_chr%d_v2.bim.gz", ukb_cal_bim_file_prefix, autosomes),
    fam.files = ukb_cal_fam_file
  ),
  reg = reg_merge_cal_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_merge_cal_ukb, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Merge-Called-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_merge_cal_ukb)

# Set family ID to UKB in called files ----
reg_set_fid_ukb_cal <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_ukb_cal"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$1=\"UKB\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file))
  }, 
  file = sprintf("%s.fam", ukb_called_genotypes_file_prefix),
  reg = reg_set_fid_ukb_cal
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_ukb_cal, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-UKB-Called", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_ukb_cal)

# Find close relatives ----
reg_kinship_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "kinship_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_fast_related, 
  input.prefix = sprintf("%s", ukb_called_genotypes_file_prefix),
  output.prefix = sprintf("%s.kinship", ukb_called_genotypes_file_prefix), 
  more.args = list( 
    plink2.exec = plink2_exec,
    degree = 0.0884, # 2nd degree
    screen.degree = 0.0442, # 3rd degree
    screen.thin = 0.1,
    screen.maf = 0.2,
    keep = sprintf("--keep %s.fam", ukb_converted_file_prefix)
  ),
  reg = reg_kinship_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_kinship_ukb, 
  resources = list(
    ntasks = 1, ncpus = 60, memory = 90000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Kinship-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_kinship_ukb)

# Find maximum independent set of close relatives
ukb_high_kinship <- data.table::fread(
  file = sprintf("%s.kinship.kin0", ukb_called_genotypes_file_prefix), 
  colClasses = c("character", "character", "character", "character", "integer", "numeric", "numeric", "numeric")
)

reg_mis_ukb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "mis_ukb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table", "igraph")
)

ids <- batchtools::batchMap(
  fun = find_maximum_independet_set, 
  args = list(kinships = list(ukb_high_kinship)),
  reg = reg_mis_ukb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_mis_ukb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_MIS-UKB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_mis_ukb)

samples_to_keep <- data.table::data.table(
  ID = reduceResultsList(reg = reg_mis_ukb)[[1]]
)
samples_to_remove <- batchtools::ajoin(
  x = ukb_high_kinship[, .(FID = c(`#FID1`, FID2), ID = c(ID1, ID2))],
  y = samples_to_keep, by = c("ID" = "ID")
)
data.table::fwrite(
  x = unique(samples_to_remove),
  file = ukb_remove_closely_relateds_file, 
  quote = FALSE, 
  sep = " ", 
  col.names = FALSE
)
