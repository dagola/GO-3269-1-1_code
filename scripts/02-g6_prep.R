
# Preparation of GerMIFS6 data

source("init.R")

# Fix GerMIFS6 sample file ----
imbs::system_call(
  bin = "cut",
  args = c(
    "-f1,2,3,7,8",
    "-d' '",
    g6_sample_file,
    ">", file.path(proc_dir, "g6.sample")
  )
)

# Convert GerMIFS6 to PLINK ----
reg_convert_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_g6"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_gen2bed_conversion, 
  gen.file = sprintf("%s_chr%d.gen", g6_file_prefix, autosomes),
  sample.file = file.path(proc_dir, "g6.sample"),
  output.prefix = sprintf("%s.chr%02d", g6_converted_file_prefix, autosomes), 
  chr = sprintf("--oxford-single-chr %d", autosomes),
  more.args = list( 
    extract = sprintf("--extract range %s", ukb_eb_plink_range_file),
    snps.only = "--snps-only",
    hardcall.threshold = 0.49999,
    import.dosage.certainty = 0,
    plink2.exec = plink2_exec
  ),
  reg = reg_convert_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_g6, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 20000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_g6)

# Merge chromosomes to one single file ----
reg_concat_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "concat_g6"), 
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
  output.prefix = sprintf("%s", g6_converted_file_prefix), 
  more.args = list( 
    bed.files = sprintf("%s.chr%02d.bed", g6_converted_file_prefix, autosomes),
    bim.files = sprintf("%s.chr%02d.bim", g6_converted_file_prefix, autosomes),
    fam.files = sprintf("%s.chr%02d.fam", g6_converted_file_prefix, autosomes)
  ),
  reg = reg_concat_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_concat_g6, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Concat-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_concat_g6)

# Set variant IDs to CHR:POS format ----
reg_set_ids_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_g6"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file))
  }, 
  file = sprintf("%s.bim", g6_converted_file_prefix),
  reg = reg_set_ids_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_g6, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_g6)

# Set family ID to G6 in merged imputed file ----
reg_set_fid_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_g6"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(command = "awk", args = c("'{$1=\"G6\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file))
  }, 
  file = sprintf("%s.fam", g6_converted_file_prefix),
  reg = reg_set_fid_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_g6, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_g6)

# Find close relatives ----
reg_kinship_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "kinship_g6"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = TRUE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_fast_related, 
  input.prefix = sprintf("%s", g6_converted_file_prefix),
  output.prefix = sprintf("%s.kinship", g6_converted_file_prefix), 
  more.args = list( 
    plink2.exec = plink2_exec,
    degree = 0.0884, # 2nd degree
    screen.degree = 0.0442, # 3rd degree
    screen.thin = 0.1,
    screen.maf = 0.2
  ),
  reg = reg_kinship_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_kinship_g6, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 20000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Kinship-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_kinship_g6)

# Find maximum independent set of close relatives
g6_high_kinship <- data.table::fread(
  file = sprintf("%s.kinship.kin0", g6_converted_file_prefix), 
  colClasses = c("character", "character", "character", "character", "integer", "numeric", "numeric", "numeric")
)

reg_mis_g6 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "mis_g6"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table", "igraph")
)

ids <- batchtools::batchMap(
  fun = find_maximum_independet_set, 
  args = list(kinships = list(g6_high_kinship)),
  reg = reg_mis_g6
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_mis_g6, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_MIS-G6", project_name)
  )
)

batchtools::waitForJobs(reg = reg_mis_g6)

samples_to_keep <- data.table::data.table(
  ID = reduceResultsList(reg = reg_mis_g6)[[1]]
)
samples_to_remove <- batchtools::ajoin(
  x = g6_high_kinship[, .(FID = c(`#FID1`, FID2), ID = c(ID1, ID2))],
  y = samples_to_keep, by = c("ID" = "ID")
)
data.table::fwrite(
  x = unique(samples_to_remove),
  file = g6_remove_closely_relateds_file, 
  quote = FALSE, 
  sep = " ", 
  col.names = FALSE
)
