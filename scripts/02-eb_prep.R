
# Preparation of EB data

source("init.R")

# Convert EB to PLINK ----
reg_convert_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_vcf2bed_conversion, 
  vcf.file = sprintf("%s_chr%d.vcf.gz", eb_file_prefix, autosomes),
  output.prefix = sprintf("%s.chr%02d", eb_converted_file_prefix, autosomes), 
  more.args = list( 
    extract = sprintf("--extract range %s", ukb_eb_plink_range_file),
    hardcall.threshold = 0.49999,
    plink2.exec = plink2_exec
  ),
  reg = reg_convert_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_eb, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 60000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_eb)

# Merge chromosomes to one single file ----
reg_concat_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "concat_eb"), 
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
    imbs::system_call(bin = "bash", sprintf("-c '(echo -en \"\\x6C\\x1B\\x01\"; for g in %s; do tail -qc +4 $g; done) > %s.bed'", paste(bed.files, collapse = " "),output.prefix))
    
    cat("Done.")
    
  }, 
  output.prefix = sprintf("%s", eb_converted_file_prefix), 
  more.args = list( 
    bed.files = sprintf("%s.chr%02d.bed", eb_converted_file_prefix, autosomes),
    bim.files = sprintf("%s.chr%02d.bim", eb_converted_file_prefix, autosomes),
    fam.files = sprintf("%s.chr%02d.fam", eb_converted_file_prefix, autosomes)
  ),
  reg = reg_concat_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_concat_eb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Concat-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_concat_eb)

# Set variant IDs to CHR:POS format ----
reg_set_ids_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = TRUE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file))
  }, 
  file = sprintf("%s.bim", eb_converted_file_prefix),
  reg = reg_set_ids_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_eb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_eb)

# Set family ID to EB ----
reg_set_fid_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = TRUE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$1=\"EB\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file))
  }, 
  file = sprintf("%s.fam", eb_converted_file_prefix),
  reg = reg_set_fid_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_eb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_eb)

# Find close relatives ----
reg_kinship_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "kinship_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_fast_related, 
  input.prefix = sprintf("%s", eb_converted_file_prefix),
  output.prefix = sprintf("%s.kinship", eb_converted_file_prefix), 
  more.args = list( 
    plink2.exec = plink2_exec,
    degree = 0.0884, # 2nd degree
    screen.degree = 0.0442, # 3th degree
    screen.thin = 0.1,
    screen.maf = 0.2
  ),
  reg = reg_kinship_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_kinship_eb, 
  resources = list(
    ntasks = 1, ncpus = 60, memory = 90000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Kinship-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_kinship_eb)

# Find maximum independent set of close relatives
eb_high_kinship <- data.table::fread(
  file = sprintf("%s.kinship.kin0", eb_converted_file_prefix), 
  colClasses = c("character", "character", "character", "character", "integer", "numeric", "numeric", "numeric")
)

reg_mis_eb <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "mis_eb"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table", "igraph")
)

ids <- batchtools::batchMap(
  fun = find_maximum_independet_set, 
  args = list(kinships = list(eb_high_kinship)),
  reg = reg_mis_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_mis_eb, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_MIS-EB", project_name)
  )
)

batchtools::waitForJobs(reg = reg_mis_eb)

samples_to_keep <- data.table::data.table(
  ID = reduceResultsList(reg = reg_mis_eb)[[1]]
)
samples_to_remove <- batchtools::ajoin(
  x = eb_high_kinship[, .(FID = c(`#FID1`, FID2), ID = c(ID1, ID2))],
  y = samples_to_keep, by = c("ID" = "ID")
)
data.table::fwrite(
  x = unique(samples_to_remove),
  file = eb_remove_closely_relateds_file, 
  quote = FALSE, 
  sep = " ", 
  col.names = FALSE
)
