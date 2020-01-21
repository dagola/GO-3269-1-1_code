
# Preparation of GerMIFS7 data
# Includes imputation against 1000G reference

source("init.R")

# Harmonize with reference ----
reg_harmonize_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "harmonize_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "checkmate")
)

ids <- batchtools::batchMap(
  fun = imbs::harmonize_genotypes, 
  chr.filter = autosomes,
  output = sprintf("%s_chr%02d", g7_converted_file_prefix, autosomes),
  force.chr = autosomes,
  more.args = list( 
    input = g7_file_prefix,
    ref = g1k_imputation_reference_file_prefix,
    input.type = "PLINK_BED",
    ref.type = "VCF",
    output.type = "PLINK_BED",
    update.id = TRUE,
    min.ld = 0.6,
    min.variants = 5,
    variants = 500,
    check.ld = FALSE,
    maf.align = 0.1,
    update.reference.allele = TRUE,
    ambiguous.snp.filter = FALSE,
    keep = TRUE,
    exec = genotype_harmonizer_exec
  ),
  reg = reg_harmonize_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_harmonize_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 8000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Harmonize-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_harmonize_g7)

# Phase chromosomes ----
reg_phase_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "phase_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "checkmate")
)

ids <- batchtools::batchMap(
  fun = shapeit,
  input.file.prefix = sprintf("%s_chr%02d", g7_converted_file_prefix, autosomes), 
  input.type = "PLINK_BED", 
  map.file = sprintf(genetic_map_file_template, autosomes), 
  output.file.prefix = sprintf("%s_chr%02d", g7_converted_file_prefix, autosomes),
  more.args = list( 
    noped = TRUE,
    burn = phasing_burn_iterations,
    prune = phasing_prune_iterations,
    main = phasing_main_iterations, 
    states = phasing_states,
    window = phasing_window_size,
    effective.size = phasing_effective_size,
    force = TRUE,
    exec = shapeit2_exec
  ),
  reg = reg_phase_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_phase_g7, 
  resources = list(
    ntasks = 1, ncpus = 8, memory = 60000,
    partition = "batch,prio", walltime = 0,
    account = "GO-3269-1-1",
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Phase-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_phase_g7)

# Impute chunks ----
chromosome_lenghts <- "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes"

reg_impute_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "impute_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

parameters <- CJ(CHR = autosomes, interval = lapply(BBmisc::chunk(1:max(data.table::fread(chromosome_lenghts)$V2), chunk.size = impute2_width), function(chunk) c(min(chunk), max(chunk))))
parameters[, study.haplotypes.file := sprintf("%s_chr%02d.haps", g7_converted_file_prefix, CHR)]
parameters[, genetic.map.file := sprintf(genetic_map_file_template, CHR)]
parameters[, output.file := sprintf("%s_chr%02d_%s.impute2", g7_converted_file_prefix, CHR, lapply(interval, paste, collapse = "_"))]
parameters[, reference.haplotypes.file := sprintf("%s_chr%d.hap.gz", g1k_imputation_reference_file_prefix, CHR)]
parameters[, reference.legend.file := sprintf("%s_chr%d.legend.gz", g1k_imputation_reference_file_prefix, CHR)]
parameters[, CHR := NULL]

ids <- batchtools::batchMap(
  fun = impute, 
  args = parameters,
  more.args = list(
    gz = TRUE,
    buffer = impute2_buffer,
    filter.rules = impute2_filter_rules,
    ref.template.haplotypes = impute2_hap_imputing,
    verbose = TRUE,
    exec = impute2_exec
  ),
  reg = reg_impute_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_impute_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 32000,
    partition = "batch,prio", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    account = "GO-3269-1-1",
    name = sprintf("%s_Impute-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_impute_g7)

# Combine imputed chromosome chunks ----
reg_combine_chunks_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "combine_chunks_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

imputed_chunks <- list.files(
  path = proc_dir, 
  pattern = "g7.*impute2.gz", 
  full.names = TRUE
)
imputed_chunks <- split(
  x = imputed_chunks, 
  f = gsub(".*/g7_chr([0-9]+).*impute2\\.gz", "\\1", imputed_chunks)
)

ids <- batchtools::batchMap(
  fun = function(input.files, output.file) {
    system2(command = "cat", args = c(input.files, ">", output.file))
  },
  input.files = lapply(imputed_chunks[sprintf("%02d", autosomes)], gtools::mixedsort),
  output.file = sprintf("%s_chr%02d.impute2.gz", g7_converted_file_prefix, autosomes),
  reg = reg_combine_chunks_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_combine_chunks_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch,prio", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    account = "GO-3269-1-1",
    name = sprintf("%s_Impute-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_combine_chunks_g7)

# Convert GerMIFS7 to PLINK ----
reg_convert_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "convert_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

system2(
  command = "awk", 
  args = c(
    "'NR<=2 { ; print $1,$2,$3,$6,$7; next} { ; $7=$7-1; print $1,$2,$3,$6,$7}'",
    file.path(proc_dir, "g7_chr01.sample"), ">", file.path(proc_dir, "g7_fixed.sample")
  )
)

ids <- batchtools::batchMap(
  fun = plink2_gen2bed_conversion, 
  gen.file = sprintf("%s_chr%02d.impute2.gz", g7_converted_file_prefix, autosomes),
  sample.file = file.path(proc_dir, "g7_fixed.sample"),
  output.prefix = sprintf("%s.chr%02d", g7_converted_file_prefix, autosomes), 
  chr = sprintf("--oxford-single-chr %d", autosomes),
  more.args = list( 
    extract = sprintf("--extract range %s", ukb_eb_plink_range_file),
    snps.only = "--snps-only",
    hardcall.threshold = 0.49999,
    import.dosage.certainty = 0,
    plink2.exec = plink2_exec
  ),
  reg = reg_convert_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_convert_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 20000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Convert-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_convert_g7)

# Merge chromosomes to one single file ----
reg_concat_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "concat_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = function(bed.files, bim.files, fam.files, output.prefix) {
    
    # Copy fam file
    cat("Copying fam file...")
    imbs::system_call(
      bin = "cat", 
      args = c(fam.files[1], ">", sprintf("%s.fam", output.prefix))
    )
    
    # Combine bim files
    cat("Merging bim files...")
    imbs::system_call(
      bin = "cat", 
      args = c(paste(bim.files, collapse = " "), ">", sprintf("%s.bim", output.prefix))
    )
    
    # Combine bed files
    cat("Merging bed files...")
    imbs::system_call(
      bin = "bash", 
      args = sprintf("-c '(echo -en \"\\x6C\\x1B\\x01\"; for g in %s; do tail -qc +4 $g; done) > %s.bed'", paste(bed.files, collapse = " "), output.prefix)
    )
    
    cat("Done.")
    
  }, 
  output.prefix = sprintf("%s", g7_converted_file_prefix), 
  more.args = list( 
    bed.files = sprintf("%s.chr%02d.bed", g7_converted_file_prefix, autosomes),
    bim.files = sprintf("%s.chr%02d.bim", g7_converted_file_prefix, autosomes),
    fam.files = sprintf("%s.chr%02d.fam", g7_converted_file_prefix, autosomes)
  ),
  reg = reg_concat_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_concat_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 6000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Concat-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_concat_g7)

# Set variant IDs to CHR:POS format ----
reg_set_ids_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_ids_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(
      command = "awk", 
      args = c(
        "'{$2=$1\":\"$4; print $0}'", file, ">", tmp, "&&", "mv", tmp, file
      )
    )
  }, 
  file = sprintf("%s.bim", g7_converted_file_prefix),
  reg = reg_set_ids_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_ids_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-IDs-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_ids_g7)

# Set family ID to G7 in merged imputed file ----
reg_set_fid_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "set_fid_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE
)

ids <- batchtools::batchMap(
  fun = function(file) {
    tmp <- tempfile()
    system2(
      command = "awk", 
      args = c(
        "'{$1=\"G7\"; print $0}'", file, ">", tmp, "&&", "cp", tmp, file
      )
    )
  }, 
  file = sprintf("%s.fam", g7_converted_file_prefix),
  reg = reg_set_fid_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_set_fid_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 1000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Set-FID-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_set_fid_g7)

# Find close relatives ----
reg_kinship_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "kinship_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = plink2_fast_related, 
  input.prefix = sprintf("%s", g7_converted_file_prefix),
  output.prefix = sprintf("%s.kinship", g7_converted_file_prefix), 
  more.args = list( 
    plink2.exec = plink2_exec,
    degree = 0.0884, # 2nd degree
    screen.degree = 0.0442, # 3rd degree
    screen.thin = 0.1,
    screen.maf = 0.2
  ),
  reg = reg_kinship_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_kinship_g7, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 20000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Kinship-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_kinship_g7)

# Find maximum independent set of close relatives
g7_high_kinship <- data.table::fread(
  file = sprintf("%s.kinship.kin0", g7_converted_file_prefix), 
  colClasses = c("character", "character", "character", "character", "integer", "numeric", "numeric", "numeric")
)

reg_mis_g7 <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "mis_g7"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table", "igraph")
)

ids <- batchtools::batchMap(
  fun = find_maximum_independet_set, 
  args = list(kinships = list(g7_high_kinship)),
  reg = reg_mis_g7
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_mis_g7, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 0,
    atonce = 1,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_MIS-G7", project_name)
  )
)

batchtools::waitForJobs(reg = reg_mis_g7)

samples_to_keep <- data.table::data.table(
  ID = reduceResultsList(reg = reg_mis_g7)[[1]]
)
samples_to_remove <- batchtools::ajoin(
  x = g7_high_kinship[, .(FID = c(`#FID1`, FID2), ID = c(ID1, ID2))],
  y = samples_to_keep, by = c("ID" = "ID")
)
data.table::fwrite(
  x = unique(samples_to_remove),
  file = g7_remove_closely_relateds_file, 
  quote = FALSE, 
  sep = " ", 
  col.names = FALSE
)
