
source("init.R")

# Select samples from UKB and EB for PCA enrichment ----
# Get sample list
ukb_samples <- data.table::fread(
  file = sprintf("%s.fam", ukb_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS")
)

eb_samples <- data.table::fread(
  file = sprintf("%s.fam", eb_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS")
)

de_samples <- data.table::fread(
  file = sprintf("%s.fam", de_converted_file_prefix),
  header = FALSE,
  col.names = c("FID", "IID", "PID", "MID", "SEX", "STATUS")
)
# Get close relatives sample list
ukb_close_samples <- data.table::fread(
  file = ukb_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID")
)

eb_close_samples <- data.table::fread(
  file = eb_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID")
)

de_close_samples <- data.table::fread(
  file = de_remove_closely_relateds_file,
  header = FALSE,
  col.names = c("FID", "IID")
)

# Choose random samples for PCA enrichment
ukb_pca_samples <- ajoin(ukb_samples, ukb_close_samples, by = c("FID", "IID"))[sample(1:.N, additional_pca_samples)]
data.table::fwrite(
  x = ukb_pca_samples[, .(FID, IID)],
  file = ukb_pca_samples_file,
  sep = " ", quote = FALSE,
  col.names = FALSE
)
eb_pca_samples <- ajoin(eb_samples, eb_close_samples, by = c("FID", "IID"))[sample(1:.N, additional_pca_samples)]
data.table::fwrite(
  x = eb_pca_samples[, .(FID, IID)],
  file = eb_pca_samples_file,
  sep = " ", quote = FALSE,
  col.names = FALSE
)
de_pca_samples <- ajoin(de_samples, de_close_samples, by = c("FID", "IID"))[sample(1:.N, additional_pca_samples)]
data.table::fwrite(
  x = de_pca_samples[, .(FID, IID)],
  file = de_pca_samples_file,
  sep = " ", quote = FALSE,
  col.names = FALSE
)

# Extract additional PCA samples
reg_pca_subset <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_subset"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = imbs::plink_subset, 
  bfile = c(ukb_converted_file_prefix, eb_converted_file_prefix, de_converted_file_prefix), 
  keep = c(ukb_pca_samples_file, eb_pca_samples_file, de_pca_samples_file), 
  output.prefix = c(ukb_pca_file_prefix, eb_pca_file_prefix, de_pca_file_prefix), 
  more.args = list(
    extract = ukb_eb_de_nikpay_g1k_snps_file, 
    exec = plink_exec
  ),
  reg = reg_pca_subset
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_subset, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 40000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Subset", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_subset)

# Merge 1kG, UKB and EB samples ----
# Define list of files to merge into one dataset
merge_list <- c(
  "1kG" = g1k_converted_file_prefix, 
  "UKB" = ukb_pca_file_prefix, 
  "EB" = eb_pca_file_prefix,
  "DE" = de_pca_file_prefix
)
writeLines(text = merge_list, con = pca_merge_list_file)

# Define list of EUR samples to keep
g1k_pop_info <- data.table::fread(g1k_pop_info_file)
pca_EUR_samples <- rbind(g1k_pop_info[super_pop == "EUR", .(FID = sample, IID = sample)], ukb_pca_samples[, .(FID = FID, IID = IID)], eb_pca_samples[, .(FID = FID, IID = IID)], de_pca_samples[, .(FID = FID, IID = IID)])
data.table::fwrite(
  x = pca_EUR_samples,
  file = pca_EUR_samples_file, 
  quote = FALSE,
  sep = " ", 
  col.names = FALSE
)

# Perform merge of 1kG and additional UKB and EB samples for PCA
reg_pca_merge <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_merge"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = imbs::plink_merge_list, 
  keep = c("", sprintf("--keep %s", pca_EUR_samples_file)),
  output.prefix = c(pca_data_file_prefix, sprintf("%s.EUR", pca_data_file_prefix)),
  more.args = list(
    merge.list = pca_merge_list_file,
    extract = sprintf("--extract %s", ukb_eb_nikpay_g1k_snps_file),
    exec = plink_exec
  ),
  reg = reg_pca_merge
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_merge, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 10000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Merge", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_merge)

# LD pruning ----
reg_pca_pruning <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_pruning"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = imbs::plink_ld_pruning, 
  args = data.table(
    bfile = c(rep(pca_data_file_prefix, times = length(autosomes)), rep(sprintf("%s.EUR", pca_data_file_prefix), times = length(autosomes))),
    output.prefix = c(sprintf("%s.chr%02d", pca_data_file_prefix, autosomes), sprintf("%s.EUR.chr%02d", pca_data_file_prefix, autosomes)),
    chr = autosomes
  ),
  more.args = list(
    exec = plink_exec,
    window.size = ld_window_size,
    step.size = ld_step_size,
    threshold = ld_threshold,
    maf.filter = sprintf("--maf %f", pca_min_maf),
    exclude.long.range.ld = sprintf("--exclude range %s", long_range_ld_ranges_file)
  ),
  reg = reg_pca_pruning
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_pruning, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 8000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Pruning", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_pruning)

# Combine per-chromosome LD pruning results
system2(command = "cat", args = c(sprintf("%s.chr%02d.prune.in", pca_data_file_prefix, autosomes), ">", sprintf("%s.prune.in", pca_data_file_prefix)))
system2(command = "cat", args = c(sprintf("%s.EUR.chr%02d.prune.in", pca_data_file_prefix, autosomes), ">", sprintf("%s.EUR.prune.in", pca_data_file_prefix)))

# Create GCTA GRM ----
reg_pca_make_grms <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_make_grms"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_make_grm,
  args = data.table(
    bfile = c(rep(pca_data_file_prefix, times = length(autosomes)), rep(sprintf("%s.EUR", pca_data_file_prefix), times = length(autosomes))),
    output.prefix = c(sprintf("%s.chr%02d", pca_data_file_prefix, autosomes), sprintf("%s.EUR.chr%02d", pca_data_file_prefix, autosomes)),
    extract = c(sprintf("%s.chr%02d.prune.in", pca_data_file_prefix, autosomes), sprintf("%s.EUR.chr%02d.prune.in", pca_data_file_prefix, autosomes)),
    chr = autosomes
  ),
  more.args = list(
    exec = gcta_exec
  ),
  reg = reg_pca_make_grms
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_make_grms, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 8000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Make-GRMs", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_make_grms)

writeLines(text = sprintf("%s.chr%02d", pca_data_file_prefix, autosomes), con = grms_list_file)
writeLines(text = sprintf("%s.EUR.chr%02d", pca_data_file_prefix, autosomes), con = grms_list_EUR_file)

reg_pca_merge_grms <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_merge_grms"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_merge_grms,
  grms.list = c(grms_list_file, grms_list_EUR_file),
  output.prefix = c(pca_data_file_prefix, sprintf("%s.EUR", pca_data_file_prefix)),
  more.args = list(
    exec = gcta_exec
  ),
  reg = reg_pca_merge_grms
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_merge_grms, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 8000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Merge-GRMs", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_merge_grms)

# Calculate GCTA PCs and PC loadings ----
reg_pca <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_pca,
  grm.prefix = c(pca_data_file_prefix, sprintf("%s.EUR", pca_data_file_prefix)),
  bfile = c(pca_data_file_prefix, sprintf("%s.EUR", pca_data_file_prefix)),
  output.prefix = c(pca_data_file_prefix, sprintf("%s.EUR", pca_data_file_prefix)),
  extract = c(sprintf("--extract %s.prune.in", pca_data_file_prefix), sprintf("--extract %s.EUR.prune.in", pca_data_file_prefix)),
  num.pcs = num_pcs,
  more.args = list(
    exec = gcta_exec
  ),
  reg = reg_pca
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca, 
  resources = list(
    ntasks = 1, ncpus = 60, memory = 90000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca)

# PCA projection
reg_pca_projection <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_projection"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_pca_projection,
  args = data.table(
    bfile = rep(c(ukb_converted_file_prefix, eb_converted_file_prefix, de_converted_file_prefix), each = 2),
    output.prefix = c(sprintf(c("%s", "%s.EUR"), ukb_converted_file_prefix), sprintf(c("%s", "%s.EUR"), eb_converted_file_prefix), sprintf(c("%s", "%s.EUR"), de_converted_file_prefix)),
    pca.loadings.prefix = rep(c(sprintf(c("%s", "%s.EUR"), pca_data_file_prefix)), times = 3),
    extract = rep(c(sprintf("--extract %s.prune.in", pca_data_file_prefix), sprintf("--extract %s.EUR.prune.in", pca_data_file_prefix)), times = 3)
  ),
  more.args = list(
    num.pcs = num_pcs,
    exec = gcta_exec
  ),
  reg = reg_pca_projection
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_projection, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 31000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Projection", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_projection)

# Make PCA plots ----
pca_eigenvec <- data.table::fread(sprintf("%s.eigenvec", pca_data_file_prefix))
pca_EUR_eigenvec <- data.table::fread(sprintf("%s.EUR.eigenvec", pca_data_file_prefix))
ukb_proj <- data.table::fread(sprintf("%s.proj.eigenvec", ukb_converted_file_prefix))
ukb_EUR_proj <- data.table::fread(sprintf("%s.EUR.proj.eigenvec", ukb_converted_file_prefix))
eb_proj <- data.table::fread(sprintf("%s.proj.eigenvec", eb_converted_file_prefix))
eb_EUR_proj <- data.table::fread(sprintf("%s.EUR.proj.eigenvec", eb_converted_file_prefix))
de_proj <- data.table::fread(sprintf("%s.proj.eigenvec", de_converted_file_prefix))
de_EUR_proj <- data.table::fread(sprintf("%s.EUR.proj.eigenvec", de_converted_file_prefix))

g1k_pop_info <- data.table::fread(g1k_pop_info_file, header = FALSE, col.names = c("IID", "POP", "GROUP", "SEX"))
ukb_pop_info <- data.table::fread(ukb_pca_samples_file, header = FALSE, col.names = c("FID", "IID"))
eb_pop_info <- data.table::fread(eb_pca_samples_file, header = FALSE, col.names = c("FID", "IID"))
de_pop_info <- data.table::fread(de_pca_samples_file, header = FALSE, col.names = c("FID", "IID"))
pop_info <- rbind(
  g1k_pop_info[, .(ID = IID, POP = POP, GROUP = GROUP)], 
  ukb_pop_info[, .(ID = IID, POP = "UKB", GROUP = "EUR")], 
  eb_pop_info[, .(ID = IID, POP = "EB", GROUP = "EUR")], 
  de_pop_info[, .(ID = IID, POP = "DE", GROUP = "EUR")], 
  fill = TRUE
)

pca_plot_data <- pop_info[pca_eigenvec, on = c("ID" = "V2")]
pca_EUR_plot_data <- pop_info[pca_EUR_eigenvec, on = c("ID" = "V2")]

pca_plot <- ggplot(data = pca_plot_data) + 
  geom_point(aes(x = V3, y = V4, color = POP), alpha = 0.25) + # Base PCA plot
  geom_point(aes(x = V3, y = V4, color = "UKB_PROJ"), data = ukb_proj[sample(.N, 1000)]) + # UKB projections
  geom_point(aes(x = V3, y = V4, color = "EB_PROJ"), data = eb_proj[sample(.N, 1000)]) + # EB projections
  geom_point(aes(x = V3, y = V4, color = "DE_PROJ"), data = de_proj[sample(.N, 1000)]) # DE projections

cowplot::ggsave(plot = pca_plot, filename = pca_plot_file, width = 12, height = 12)

pca_EUR_plot <- ggplot(data = pca_EUR_plot_data) + 
  geom_point(aes(x = V3, y = V4, color = POP), alpha = 0.25) + # Base PCA plot
  geom_point(aes(x = V3, y = V4, color = "UKB_PROJ"), data = ukb_EUR_proj[sample(.N, 1000)]) + # UKB projections
  geom_point(aes(x = V3, y = V4, color = "EB_PROJ"), data = eb_EUR_proj[sample(.N, 1000)]) + # EB projections
  geom_point(aes(x = V3, y = V4, color = "DE_PROJ"), data = de_EUR_proj[sample(.N, 1000)]) # DE projections

cowplot::ggsave(plot = pca_EUR_plot, filename = pca_EUR_plot_file, width = 12, height = 12)

# Per dataset PCA ----
reg_pca_per_dataset <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "pca_per_dataset"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "checkmate")
)

ids <- batchtools::batchMap(
  fun = plink2_pca, 
  input.prefix = c(ukb_converted_file_prefix, eb_converted_file_prefix, de_converted_file_prefix), 
  remove = sprintf("--remove %s", c(ukb_remove_closely_relateds_file, eb_remove_closely_relateds_file, de_remove_closely_relateds_file)), 
  output.prefix = c(ukb_pca_file_prefix, eb_pca_file_prefix, de_pca_file_prefix), 
  more.args = list(
    extract = sprintf("--extract %s.prune.in", pca_data_file_prefix), 
    approx = TRUE,
    num.pcs = num_pcs,
    plink2.exec = plink2_exec
  ),
  reg = reg_pca_per_dataset
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_pca_per_dataset, 
  resources = list(
    ntasks = 1, ncpus = 20, memory = 40000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_PCA-Datasets", project_name)
  )
)

batchtools::waitForJobs(reg = reg_pca_per_dataset)

# G6 ----

# PCA projection
reg_g6_pca_projection <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "g6_pca_projection"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_pca_projection,
  args = data.table(
    bfile = rep(g6_converted_file_prefix, each = 2),
    output.prefix = sprintf(c("%s", "%s.EUR"), g6_converted_file_prefix),
    pca.loadings.prefix = sprintf(c("%s", "%s.EUR"), pca_data_file_prefix),
    extract = c(sprintf("--extract %s.prune.in", pca_data_file_prefix), sprintf("--extract %s.EUR.prune.in", pca_data_file_prefix))
  ),
  more.args = list(
    num.pcs = num_pcs,
    exec = gcta_exec
  ),
  reg = reg_g6_pca_projection
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_g6_pca_projection, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 91000,
    partition = "prio", walltime = 0,
    account = project_name,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_G6-PCA-Projection", project_name)
  )
)

batchtools::waitForJobs(reg = reg_g6_pca_projection)

# PCA on single dataset
reg_g6_pca <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "g6_pca"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "checkmate")
)

ids <- batchtools::batchMap(
  fun = plink2_pca, 
  input.prefix = g6_converted_file_prefix, 
  remove = sprintf("--remove %s", g6_remove_closely_relateds_file), 
  output.prefix = g6_pca_file_prefix, 
  more.args = list(
    extract = sprintf("--extract %s.prune.in", pca_data_file_prefix), 
    approx = FALSE,
    num.pcs = num_pcs,
    plink2.exec = plink2_exec
  ),
  reg = reg_g6_pca
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_g6_pca, 
  resources = list(
    ntasks = 1, ncpus = 20, memory = 40000,
    partition = "prio", walltime = 0,
    account = project_name,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_G6-PCA", project_name)
  )
)

batchtools::waitForJobs(reg = reg_g6_pca)

# G7 ----

# PCA projection
reg_g7_pca_projection <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "g7_pca_projection"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs")
)

ids <- batchtools::batchMap(
  fun = gcta_pca_projection,
  args = data.table(
    bfile = rep(g7_converted_file_prefix, each = 2),
    output.prefix = sprintf(c("%s", "%s.EUR"), g7_converted_file_prefix),
    pca.loadings.prefix = sprintf(c("%s", "%s.EUR"), pca_data_file_prefix),
    extract = c(sprintf("--extract %s.prune.in", pca_data_file_prefix), sprintf("--extract %s.EUR.prune.in", pca_data_file_prefix))
  ),
  more.args = list(
    num.pcs = num_pcs,
    exec = gcta_exec
  ),
  reg = reg_g7_pca_projection
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_g7_pca_projection, 
  resources = list(
    ntasks = 1, ncpus = 10, memory = 91000,
    partition = "prio", walltime = 0,
    account = project_name,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_G7-PCA-Projection", project_name)
  )
)

batchtools::waitForJobs(reg = reg_g7_pca_projection)

# PCA on single dataset
reg_g7_pca <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "g7_pca"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("imbs", "checkmate")
)

ids <- batchtools::batchMap(
  fun = plink2_pca, 
  input.prefix = g7_converted_file_prefix, 
  remove = sprintf("--remove %s", g7_remove_closely_relateds_file), 
  output.prefix = g7_pca_file_prefix, 
  more.args = list(
    extract = sprintf("--extract %s.prune.in", pca_data_file_prefix), 
    approx = FALSE,
    num.pcs = num_pcs,
    plink2.exec = plink2_exec
  ),
  reg = reg_g7_pca
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_g7_pca, 
  resources = list(
    ntasks = 1, ncpus = 20, memory = 40000,
    partition = "prio", walltime = 0,
    account = project_name,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_G7-PCA", project_name)
  )
)

batchtools::waitForJobs(reg = reg_g7_pca)


