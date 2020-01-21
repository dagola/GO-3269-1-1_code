
source("init.R")

# Define measures ----
measures <- list(aucpr, # use Area under Precision-Recall Curve as performance measure during tuning
                 mlr::mmce, mlr::npv, mlr::fpr, mlr::f1, mlr::fnr, mlr::ssr,
                 mlr::tp, mlr::tn, mlr::gpr, mlr::lsr, mlr::acc,
                 mlr::wkappa, mlr::ppv, mlr::logloss, mlr::ber, mlr::tpr,
                 mlr::brier, mlr::gmean, mlr::fdr, mlr::tnr, mlr::qsr, mlr::bac,
                 mlr::brier.scaled, mlr::fp, mlr::fn, mlr::kappa, mlr::auc)

# Load tasks ----
ukb_nikpay_task <- readRDS(ukb_nikpay_task_file)
ukb_task <- readRDS(ukb_task_file)
de_nikpay_task <- readRDS(de_nikpay_task_file)
de_task <- readRDS(de_task_file)
de_downsampled_nikpay_task <- readRDS(de_downsampled_nikpay_task_file)
de_downsampled_task <- readRDS(de_downsampled_task_file)
eb_nikpay_task <- readRDS(eb_nikpay_task_file)
eb_task <- readRDS(eb_task_file)
combined_nikpay_task <- readRDS(combined_nikpay_task_file)
combined_task <- readRDS(combined_task_file)

# Define MBO control ----
mbo_ctrl <- mlrMBO::makeMBOControl(
  impute.y.fun = function(x, y, opt.path, ...) 0, # This is the worst AUCPR,
  final.method = "best.predicted" # better if target function is noisy
)
mbo_ctrl <- mlrMBO::setMBOControlTermination(
  control = mbo_ctrl,
  iters = tuning_iters, #tuning_iters
  time.budget = tuning_time_budget #tuning_time_budget
)
mbo_ctrl <- mlrMBO::setMBOControlInfill(
  control = mbo_ctrl, 
  crit = mlrMBO::makeMBOInfillCritAEI()
)
surrogate_lrn <- mlr::makeImputeWrapper(
  learner = mlr::makeLearner("regr.ranger", predict.type = "se", replace = FALSE, keep.inbag = TRUE),
  classes = list(
    numeric = mlr::imputeMax(2),
    factor = mlr::imputeConstant("__miss__")
  )
)
ctrl <- mlr:::makeTuneControlMBO(
  learner = surrogate_lrn,
  mbo.control = mbo_ctrl
)

# Define resampling strategies ----
inner <- mlr::makeResampleDesc("CV", stratify = tuning_stratify, iters = tuning_cv_iters)
outer <- list(
  "UKB" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(ukb_task),
    test.inds = c(sample(which(ukb_task$env$data$STATUS==mlr::getTaskDesc(ukb_task)$positive), 5), sample(which(ukb_task$env$data$STATUS==mlr::getTaskDesc(ukb_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(ukb_task)
  ),
  "UKB_NIKPAY" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(ukb_nikpay_task), 
    test.inds = c(sample(which(ukb_nikpay_task$env$data$STATUS==mlr::getTaskDesc(ukb_nikpay_task)$positive), 5), sample(which(ukb_nikpay_task$env$data$STATUS==mlr::getTaskDesc(ukb_nikpay_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(ukb_nikpay_task)
  ),
  "DE" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(de_task),
    test.inds = c(sample(which(de_task$env$data$STATUS==mlr::getTaskDesc(de_task)$positive), 5), sample(which(de_task$env$data$STATUS==mlr::getTaskDesc(de_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(de_task)
  ),
  "DE_NIKPAY" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(de_nikpay_task), 
    test.inds = c(sample(which(de_nikpay_task$env$data$STATUS==mlr::getTaskDesc(de_nikpay_task)$positive), 5), sample(which(de_nikpay_task$env$data$STATUS==mlr::getTaskDesc(de_nikpay_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(de_nikpay_task)
  ),
  "DE_DOWNSAMPLED" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(de_downsampled_task), 
    test.inds = c(sample(which(de_downsampled_task$env$data$STATUS==mlr::getTaskDesc(de_downsampled_task)$positive), 5), sample(which(de_downsampled_task$env$data$STATUS==mlr::getTaskDesc(de_downsampled_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(de_downsampled_task)
  ),
  "DE_DOWNSAMPLED_NIKPAY" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(de_downsampled_nikpay_task), 
    test.inds = c(sample(which(de_downsampled_nikpay_task$env$data$STATUS==mlr::getTaskDesc(de_downsampled_nikpay_task)$positive), 5), sample(which(de_downsampled_nikpay_task$env$data$STATUS==mlr::getTaskDesc(de_downsampled_nikpay_task)$positive), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(de_downsampled_nikpay_task)
  ),
  "EB" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(eb_task), 
    test.inds = c(sample(which(eb_task$env$data$STATUS==mlr::getTaskDesc(eb_task)$positive), 5), sample(which(eb_task$env$data$STATUS==mlr::getTaskDesc(eb_task)$negative), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(eb_task)
  ),
  "EB_NIKPAY" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(eb_nikpay_task), 
    test.inds = c(sample(which(eb_nikpay_task$env$data$STATUS==mlr::getTaskDesc(eb_nikpay_task)$positive), 5), sample(which(eb_nikpay_task$env$data$STATUS==mlr::getTaskDesc(eb_nikpay_task)$negative), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(eb_nikpay_task)
  ),
  "COMBINED" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(combined_task), 
    test.inds = c(sample(which(combined_task$env$data$STATUS==mlr::getTaskDesc(combined_task)$positive), 5), sample(which(combined_task$env$data$STATUS==mlr::getTaskDesc(combined_task)$negative), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(combined_task)
  ),
  "COMBINED_NIKPAY" = mlr::makeFixedHoldoutInstance(
    train.inds = 1:getTaskSize(combined_nikpay_task), 
    test.inds = c(sample(which(combined_nikpay_task$env$data$STATUS==mlr::getTaskDesc(combined_nikpay_task)$positive), 5), sample(which(combined_nikpay_task$env$data$STATUS==mlr::getTaskDesc(combined_nikpay_task)$negative), 5)), # just a dummy to be able to use batchmark
    size = getTaskSize(combined_nikpay_task)
  )
)

# Define learner ----
prsice_learner <- mlr::makeLearner(
  cl = "classif.prsice",
  predict.type = "prob", 
  par.vals = list(
    summary.statistics.chr.col = '"#CHROM"',
    summary.statistics.pos.col = "POS",
    summary.statistics.snp.id.col = "ID",
    summary.statistics.effect.allele.col = "A1",
    summary.statistics.non.effect.allele.col = "AX",
    summary.statistics.stat.col = "BETA",
    summary.statistics.pvalue.col = "P",
    summary.statistics.maf.cols = "A1_FREQ",
    summary.statistics.beta = TRUE,
    summary.statistics.index = FALSE,
    target.type = "bed",
    target.nonfounders = TRUE,
    ld.type = "bed",
    ld.file.prefix = g1k_converted_file_prefix,
    clumping.pearson = FALSE,
    model = "add",
    score = "sum",
    exclude.range = long_range_ld_ranges,
    keep.ambig = FALSE,
    scratch.dir = ramdisk,
    scratch.dir.files = c("base.file", "ld.file", "target.file"),
    memory = 1500,
    prsice.exec = prsice_exec,
    plink2.exec = plink2_exec
  ),
  config = list(on.learner.error = "warn")
)

prsice_tuning_par_set <- ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam(id = "summary.statistics.maf.thresholds", lower = 0, upper = 0.1),
  ParamHelpers::makeNumericParam(id = "target.geno", lower = 0.9, upper = 1),
  ParamHelpers::makeNumericParam(id = "target.maf", lower = 0, upper = 0.1),
  ParamHelpers::makeLogicalParam(id = "ld.external", requires = quote(clumping == TRUE)),
  ParamHelpers::makeNumericParam(id = "ld.geno", lower = 0.9, upper = 1, requires = quote(!is.na(ld.external) & clumping == TRUE & ld.external == TRUE)),
  ParamHelpers::makeNumericParam(id = "ld.maf", lower = 0, upper = 0.1, requires = quote(!is.na(ld.external) & clumping == TRUE & ld.external == TRUE)),
  ParamHelpers::makeLogicalParam(id = "clumping"),
  ParamHelpers::makeIntegerParam(id = "clumping.kb", lower = 125, upper = 5000, requires = quote(clumping == TRUE)),
  ParamHelpers::makeNumericParam(id = "clumping.r2", lower = 0.1, upper = 0.8, requires = quote(clumping == TRUE)),
  ParamHelpers::makeNumericParam(id = "pval.level", lower = 5e-8, upper = 1),
  ParamHelpers::makeDiscreteParam(id = "missing.handling", values = c("IMPUTE", "SET_ZERO", "CENTER"))
)

prsice_tuning_learner <- mlr::makeTuneWrapper(
  learner = prsice_learner, 
  resampling = inner, 
  measures = measures,
  par.set = prsice_tuning_par_set, 
  show.info = TRUE,
  control = ctrl # add mlrMBO tuning control
)

# Start batchmark ----
reg_training <- batchtools::makeExperimentRegistry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)
reg_training <- batchtools::loadRegistry(
  file.dir = file.path(registries_dir, "training"), 
  writeable = TRUE
)

mlr::batchmark(
  learners = list(
    "prsice" = prsice_tuning_learner
  ),
  tasks = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task,
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task,
    "eb_task" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  measures = measures, 
  models = TRUE,
  resamplings = outer,
  reg = reg_training
)

batch_ids <- findNotDone(reg = reg_training)
batch_ids <- ajoin(batch_ids, findRunning(reg = reg_training))
batch_ids[, chunk := 1]
batchtools::submitJobs(
  ids = batch_ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training", project_name),
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training
)

batchtools::waitForJobs(reg = reg_training)

# Reduce results ----
bmr <- reduceBatchmarkResults(reg = reg_training, keep.pred = FALSE, show.info = TRUE)

ukb_model <- bmr[["results"]][["UKB"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(ukb_model, ukb_model_file)
ukb_nikpay_model <- bmr[["results"]][["UKB_NIKPAY"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(ukb_nikpay_model, ukb_nikpay_model_file)
de_model <- bmr[["results"]][["DE"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(de_model, de_model_file)
de_nikpay_model <- bmr[["results"]][["DE_NIKPAY"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(de_nikpay_model, de_nikpay_model_file)
de_downsampled_model <- bmr[["results"]][["DE_DOWNSAMPLED"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(de_downsampled_model, de_downsampled_model_file)
de_downsampled_nikpay_model <- bmr[["results"]][["DE_DOWNSAMPLED_NIKPAY"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(de_downsampled_nikpay_model, de_downsampled_nikpay_model_file)
eb_model <- bmr[["results"]][["EB"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(eb_model, eb_model_file)
eb_nikpay_model <- bmr[["results"]][["EB_NIKPAY"]][["classif.prsice.tuned"]][["models"]][[1]]
saveRDS(eb_nikpay_model, eb_nikpay_model_file)

# DE ----
reg_training2 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training2"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"

prsice_tuning_learner_de <- prsice_tuning_learner
prsice_tuning_learner_de$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de.rds")
prsice_tuning_learner_de_nikpay <- prsice_tuning_learner
prsice_tuning_learner_de_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_nikpay.rds")
prsice_tuning_learner_de_downsampled <- prsice_tuning_learner
prsice_tuning_learner_de_downsampled$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled.rds")
prsice_tuning_learner_de_downsampled_nikpay <- prsice_tuning_learner
prsice_tuning_learner_de_downsampled_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "de" = prsice_tuning_learner_de, 
    "de_nikpay" = prsice_tuning_learner_de_nikpay, 
    "de_downsampled" = prsice_tuning_learner_de_downsampled, 
    "de_downsampled_nikpay" = prsice_tuning_learner_de_downsampled_nikpay),
  task = list(
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task
  ),
  reg = reg_training2
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 3,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training2
)

batchtools::waitForJobs(reg = reg_training2)

models <- batchtools::reduceResultsDataTable(reg = reg_training2, fun = function(res) list(res))

de_model <- models[job.id == 1, result[[1]]]
saveRDS(de_model, de_model_file)
de_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_nikpay_model, de_nikpay_model_file)
de_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_downsampled_model, de_downsampled_model_file)
de_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_downsampled_nikpay_model, de_downsampled_nikpay_model_file)

models <- batchtools::reduceResultsDataTable(reg = reg_training2, fun = function(res) list(res))

de_model <- models[job.id == 1, result[[1]]]
saveRDS(de_model, de_model_file)
de_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_nikpay_model, de_nikpay_model_file)
de_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_downsampled_model, de_downsampled_model_file)
de_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_downsampled_nikpay_model, de_downsampled_nikpay_model_file)

## DE - Only genome-wide significant ----
reg_training_de_gws <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_de_gws"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_de_gws <- prsice_tuning_learner
prsice_tuning_learner_de_gws$opt.pars$pars$pval.level$lower <- 5e-8
prsice_tuning_learner_de_gws$opt.pars$pars$pval.level$upper <- 5e-8
prsice_tuning_learner_de_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_gws.rds")
prsice_tuning_learner_de_nikpay_gws <- prsice_tuning_learner_de_gws
prsice_tuning_learner_de_nikpay_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_nikpay_gws.rds")
prsice_tuning_learner_de_downsampled_gws <- prsice_tuning_learner_de_gws
prsice_tuning_learner_de_downsampled_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_gws.rds")
prsice_tuning_learner_de_downsampled_nikpay_gws <- prsice_tuning_learner_de_gws
prsice_tuning_learner_de_downsampled_nikpay_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_nikpay_gws.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "de_gws" = prsice_tuning_learner_de_gws, 
    "de_nikpay_gws" = prsice_tuning_learner_de_nikpay_gws, 
    "de_downsampled_gws" = prsice_tuning_learner_de_downsampled_gws, 
    "de_downsampled_nikpay_gws" = prsice_tuning_learner_de_downsampled_nikpay_gws),
  task = list(
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task
  ),
  reg = reg_training_de_gws
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_de_gws
)

batchtools::waitForJobs(reg = reg_training_de_gws)

models <- batchtools::reduceResultsDataTable(reg = reg_training_de_gws, fun = function(res) res)

de_gws_model <- models[job.id == 1, result[[1]]]
saveRDS(de_gws_model, de_gws_model_file)
de_gws_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_gws_nikpay_model, de_gws_nikpay_model_file)
de_gws_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_gws_downsampled_model, de_gws_downsampled_model_file)
de_gws_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_gws_downsampled_nikpay_model, de_gws_downsampled_nikpay_model_file)

## DE - Only top 3000 ----
reg_training_de_t3000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_de_t3000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_de_t3000 <- prsice_tuning_learner
prsice_tuning_learner_de_t3000$opt.pars$pars$pval.level$lower <- 2.54e-06
prsice_tuning_learner_de_t3000$opt.pars$pars$pval.level$upper <- 2.54e-06
prsice_tuning_learner_de_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_t3000.rds")
prsice_tuning_learner_de_nikpay_t3000 <- prsice_tuning_learner_de_t3000
prsice_tuning_learner_de_nikpay_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_nikpay_t3000.rds")
prsice_tuning_learner_de_downsampled_t3000 <- prsice_tuning_learner_de_t3000
prsice_tuning_learner_de_downsampled_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_t3000.rds")
prsice_tuning_learner_de_downsampled_nikpay_t3000 <- prsice_tuning_learner_de_t3000
prsice_tuning_learner_de_downsampled_nikpay_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_nikpay_t3000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "de_t3000" = prsice_tuning_learner_de_t3000, 
    "de_nikpay_t3000" = prsice_tuning_learner_de_nikpay_t3000, 
    "de_downsampled_t3000" = prsice_tuning_learner_de_downsampled_t3000, 
    "de_downsampled_nikpay_t3000" = prsice_tuning_learner_de_downsampled_nikpay_t3000),
  task = list(
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task
  ),
  reg = reg_training_de_t3000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_de_t3000
)

batchtools::waitForJobs(reg = reg_training_de_t3000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_de_t3000, fun = function(res) res)

de_t3000_model <- models[job.id == 1, result[[1]]]
saveRDS(de_t3000_model, de_t3000_model_file)
de_t3000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_t3000_nikpay_model, de_t3000_nikpay_model_file)
de_t3000_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_t3000_downsampled_model, de_t3000_downsampled_model_file)
de_t3000_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_t3000_downsampled_nikpay_model, de_t3000_downsampled_nikpay_model_file)

## DE - Only top 300000 ----
reg_training_de_t300000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_de_t300000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_de_t300000 <- prsice_tuning_learner
prsice_tuning_learner_de_t300000$opt.pars$pars$pval.level$lower <- 4.12907e-02
prsice_tuning_learner_de_t300000$opt.pars$pars$pval.level$upper <- 4.12907e-02
prsice_tuning_learner_de_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_t300000.rds")
prsice_tuning_learner_de_nikpay_t300000 <- prsice_tuning_learner_de_t300000
prsice_tuning_learner_de_nikpay_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_nikpay_t300000.rds")
prsice_tuning_learner_de_downsampled_t300000 <- prsice_tuning_learner_de_t300000
prsice_tuning_learner_de_downsampled_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_t300000.rds")
prsice_tuning_learner_de_downsampled_nikpay_t300000 <- prsice_tuning_learner_de_t300000
prsice_tuning_learner_de_downsampled_nikpay_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_nikpay_t300000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "de_t300000" = prsice_tuning_learner_de_t300000, 
    "de_nikpay_t300000" = prsice_tuning_learner_de_nikpay_t300000, 
    "de_downsampled_t300000" = prsice_tuning_learner_de_downsampled_t300000, 
    "de_downsampled_nikpay_t300000" = prsice_tuning_learner_de_downsampled_nikpay_t300000),
  task = list(
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task
  ),
  reg = reg_training_de_t300000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_de_t300000
)

batchtools::waitForJobs(reg = reg_training_de_t300000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_de_t300000, fun = function(res) res)

de_t300000_model <- models[job.id == 1, result[[1]]]
saveRDS(de_t300000_model, de_t300000_model_file)
de_t300000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_t300000_nikpay_model, de_t300000_nikpay_model_file)
de_t300000_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_t300000_downsampled_model, de_t300000_downsampled_model_file)
de_t300000_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_t300000_downsampled_nikpay_model, de_t300000_downsampled_nikpay_model_file)

## DE - Only top 3000000 ----
reg_training_de_t3000000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_de_t3000000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_de_t3000000 <- prsice_tuning_learner
prsice_tuning_learner_de_t3000000$opt.pars$pars$pval.level$lower <- 4.904758e-01
prsice_tuning_learner_de_t3000000$opt.pars$pars$pval.level$upper <- 4.904758e-01
prsice_tuning_learner_de_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_t3000000.rds")
prsice_tuning_learner_de_nikpay_t3000000 <- prsice_tuning_learner_de_t3000000
prsice_tuning_learner_de_nikpay_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_nikpay_t3000000.rds")
prsice_tuning_learner_de_downsampled_t3000000 <- prsice_tuning_learner_de_t3000000
prsice_tuning_learner_de_downsampled_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_t3000000.rds")
prsice_tuning_learner_de_downsampled_nikpay_t3000000 <- prsice_tuning_learner_de_t3000000
prsice_tuning_learner_de_downsampled_nikpay_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_de_downsampled_nikpay_t3000000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "de_t3000000" = prsice_tuning_learner_de_t3000000, 
    "de_nikpay_t3000000" = prsice_tuning_learner_de_nikpay_t3000000, 
    "de_downsampled_t3000000" = prsice_tuning_learner_de_downsampled_t3000000, 
    "de_downsampled_nikpay_t3000000" = prsice_tuning_learner_de_downsampled_nikpay_t3000000),
  task = list(
    "de" = de_task,
    "de_nikpay" = de_nikpay_task,
    "de_downsampled" = de_downsampled_task,
    "de_downsampled_nikpay" = de_downsampled_nikpay_task
  ),
  reg = reg_training_de_t3000000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_de_t3000000
)

batchtools::waitForJobs(reg = reg_training_de_t3000000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_de_t3000000, fun = function(res) res)

de_t3000000_model <- models[job.id == 1, result[[1]]]
saveRDS(de_t3000000_model, de_t3000000_model_file)
de_t3000000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(de_t3000000_nikpay_model, de_t3000000_nikpay_model_file)
de_t3000000_downsampled_model <- models[job.id == 3, result[[1]]]
saveRDS(de_t3000000_downsampled_model, de_t3000000_downsampled_model_file)
de_t3000000_downsampled_nikpay_model <- models[job.id == 4, result[[1]]]
saveRDS(de_t3000000_downsampled_nikpay_model, de_t3000000_downsampled_nikpay_model_file)

# EB training ----
reg_training_eb <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_eb"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"

prsice_tuning_learner_eb <- prsice_tuning_learner
prsice_tuning_learner_eb$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb.rds")
prsice_tuning_learner_eb_nikpay <- prsice_tuning_learner
prsice_tuning_learner_eb_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "eb" = prsice_tuning_learner_eb, 
    "eb_nikpay" = prsice_tuning_learner_eb_nikpay
  ),
  task = list(
    "eb" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  reg = reg_training_eb
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_EB", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_eb
)

batchtools::waitForJobs(reg = reg_training_eb)

models <- batchtools::reduceResultsDataTable(reg = reg_training_eb, fun = function(res) list(res))

eb_model <- models[job.id == 1, result[[1]]]
saveRDS(eb_model, eb_model_file)
eb_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(eb_nikpay_model, eb_nikpay_model_file)

## EB - Only genome-wide significant -----
reg_training_eb_gws <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_eb_gws"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_eb_gws <- prsice_tuning_learner
prsice_tuning_learner_eb_gws$opt.pars$pars$pval.level$lower <- 5e-8
prsice_tuning_learner_eb_gws$opt.pars$pars$pval.level$upper <- 5e-8
prsice_tuning_learner_eb_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_gws.rds")
prsice_tuning_learner_eb_gws_nikpay <- prsice_tuning_learner_eb_gws
prsice_tuning_learner_eb_gws_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_gws_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "eb_gws" = prsice_tuning_learner_eb_gws, 
    "eb_gws_nikpay" = prsice_tuning_learner_eb_gws_nikpay
  ),
  task = list(
    "eb" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  reg = reg_training_eb_gws
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_EB", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_eb_gws
)

batchtools::waitForJobs(reg = reg_training_eb_gws)

models <- batchtools::reduceResultsDataTable(reg = reg_training_eb_gws, fun = function(res) list(res))

eb_gws_model <- models[job.id == 1, result[[1]]]
saveRDS(eb_gws_model, eb_gws_model_file)
eb_gws_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(eb_gws_nikpay_model, eb_gws_nikpay_model_file)

## EB - Only top 3000 -----
reg_training_eb_t3000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_eb_t3000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_eb_t3000 <- prsice_tuning_learner
prsice_tuning_learner_eb_t3000$opt.pars$pars$pval.level$lower <- 2.54e-06
prsice_tuning_learner_eb_t3000$opt.pars$pars$pval.level$upper <- 2.54e-06
prsice_tuning_learner_eb_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t3000.rds")
prsice_tuning_learner_eb_t3000_nikpay <- prsice_tuning_learner_eb_t3000
prsice_tuning_learner_eb_t3000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t3000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "eb_t3000" = prsice_tuning_learner_eb_t3000, 
    "eb_t3000_nikpay" = prsice_tuning_learner_eb_t3000_nikpay
  ),
  task = list(
    "eb" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  reg = reg_training_eb_t3000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_EB", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_eb_t3000
)

batchtools::waitForJobs(reg = reg_training_eb_t3000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_eb_t3000, fun = function(res) list(res))

eb_t3000_model <- models[job.id == 1, result[[1]]]
saveRDS(eb_t3000_model, eb_t3000_model_file)
eb_t3000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(eb_t3000_nikpay_model, eb_t3000_nikpay_model_file)

## EB - Only top 300000 -----
reg_training_eb_t300000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_eb_t300000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_eb_t300000 <- prsice_tuning_learner
prsice_tuning_learner_eb_t300000$opt.pars$pars$pval.level$lower <- 4.12907e-02
prsice_tuning_learner_eb_t300000$opt.pars$pars$pval.level$upper <- 4.12907e-02
prsice_tuning_learner_eb_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t300000.rds")
prsice_tuning_learner_eb_t300000_nikpay <- prsice_tuning_learner_eb_t300000
prsice_tuning_learner_eb_t300000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t300000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "eb_t300000" = prsice_tuning_learner_eb_t300000, 
    "eb_t300000_nikpay" = prsice_tuning_learner_eb_t300000_nikpay
  ),
  task = list(
    "eb" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  reg = reg_training_eb_t300000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_EB", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_eb_t300000
)

batchtools::waitForJobs(reg = reg_training_eb_t300000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_eb_t300000, fun = function(res) list(res))

eb_t300000_model <- models[job.id == 1, result[[1]]]
saveRDS(eb_t300000_model, eb_t300000_model_file)
eb_t300000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(eb_t300000_nikpay_model, eb_t300000_nikpay_model_file)

## EB - Only top 3000000 -----
reg_training_eb_t3000000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_eb_t3000000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_eb_t3000000 <- prsice_tuning_learner
prsice_tuning_learner_eb_t3000000$opt.pars$pars$pval.level$lower <- 4.904758e-01
prsice_tuning_learner_eb_t3000000$opt.pars$pars$pval.level$upper <- 4.904758e-01
prsice_tuning_learner_eb_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t3000000.rds")
prsice_tuning_learner_eb_t3000000_nikpay <- prsice_tuning_learner_eb_t3000000
prsice_tuning_learner_eb_t3000000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_eb_t3000000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "eb_t3000000" = prsice_tuning_learner_eb_t3000000, 
    "eb_t3000000_nikpay" = prsice_tuning_learner_eb_t3000000_nikpay
  ),
  task = list(
    "eb" = eb_task,
    "eb_nikpay" = eb_nikpay_task
  ),
  reg = reg_training_eb_t3000000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_EB", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_eb_t3000000
)

batchtools::waitForJobs(reg = reg_training_eb_t3000000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_eb_t3000000, fun = function(res) list(res))

eb_t3000000_model <- models[job.id == 1, result[[1]]]
saveRDS(eb_t3000000_model, eb_t3000000_model_file)
eb_t3000000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(eb_t3000000_nikpay_model, eb_t3000000_nikpay_model_file)

## UKB ----
reg_training_ukb <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_ukb"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"

prsice_tuning_learner_ukb <- prsice_tuning_learner
prsice_tuning_learner_ukb$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb.rds")
prsice_tuning_learner_ukb_nikpay <- prsice_tuning_learner
prsice_tuning_learner_ukb_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "ukb" = prsice_tuning_learner_ukb, 
    "ukb_nikpay" = prsice_tuning_learner_ukb_nikpay
  ),
  task = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task
  ),
  reg = reg_training_ukb
)
ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training_UKB", project_name),
    ntasks = 2,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_ukb
)

batchtools::waitForJobs(reg = reg_training_ukb)

models <- batchtools::reduceResultsDataTable(reg = reg_training_ukb, fun = function(res) res)

ukb_model <- models[job.id == 1, result[[1]]]
saveRDS(ukb_model, ukb_model_file)
ukb_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(ukb_nikpay_model, ukb_nikpay_model_file)

## UKB - Only genome-wide significant ----
reg_training_ukb_gws <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_ukb_gws"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_ukb_gws <- prsice_tuning_learner
prsice_tuning_learner_ukb_gws$opt.pars$pars$pval.level$lower <- 5e-8
prsice_tuning_learner_ukb_gws$opt.pars$pars$pval.level$upper <- 5e-8
prsice_tuning_learner_ukb_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_gws.rds")
prsice_tuning_learner_ukb_nikpay_gws <- prsice_tuning_learner_ukb_gws
prsice_tuning_learner_ukb_nikpay_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_nikpay_gws.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "ukb_gws" = prsice_tuning_learner_ukb_gws, 
    "ukb_nikpay_gws" = prsice_tuning_learner_ukb_nikpay_gws),
  task = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task
  ),
  reg = reg_training_ukb_gws
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_ukb_gws
)

batchtools::waitForJobs(reg = reg_training_ukb_gws)

models <- batchtools::reduceResultsDataTable(reg = reg_training_ukb_gws, fun = function(res) list(res))

ukb_gws_model <- models[job.id == 1, result[[1]]]
saveRDS(ukb_gws_model, ukb_gws_model_file)
ukb_gws_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(ukb_gws_nikpay_model, ukb_gws_nikpay_model_file)

## UKB - Only top 3000 ----
reg_training_ukb_t3000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_ukb_t3000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_ukb_t3000 <- prsice_tuning_learner
prsice_tuning_learner_ukb_t3000$opt.pars$pars$pval.level$lower <- 2.54e-06
prsice_tuning_learner_ukb_t3000$opt.pars$pars$pval.level$upper <- 2.54e-06
prsice_tuning_learner_ukb_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_t3000.rds")
prsice_tuning_learner_ukb_nikpay_t3000 <- prsice_tuning_learner_ukb_t3000
prsice_tuning_learner_ukb_nikpay_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_nikpay_t3000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "ukb_t3000" = prsice_tuning_learner_ukb_t3000, 
    "ukb_nikpay_t3000" = prsice_tuning_learner_ukb_nikpay_t3000),
  task = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task
  ),
  reg = reg_training_ukb_t3000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_ukb_t3000
)

batchtools::waitForJobs(reg = reg_training_ukb_t3000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_ukb_t3000, fun = function(res) list(res))

ukb_t3000_model <- models[job.id == 1, result[[1]]]
saveRDS(ukb_t3000_model, ukb_t3000_model_file)
ukb_t3000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(ukb_t3000_nikpay_model, ukb_t3000_nikpay_model_file)

## UKB - Only top 300000 ----
reg_training_ukb_t300000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_ukb_t300000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_ukb_t300000 <- prsice_tuning_learner
prsice_tuning_learner_ukb_t300000$opt.pars$pars$pval.level$lower <- 4.12907e-02
prsice_tuning_learner_ukb_t300000$opt.pars$pars$pval.level$upper <- 4.12907e-02
prsice_tuning_learner_ukb_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_t300000.rds")
prsice_tuning_learner_ukb_nikpay_t300000 <- prsice_tuning_learner_ukb_t300000
prsice_tuning_learner_ukb_nikpay_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_nikpay_t300000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "ukb_t300000" = prsice_tuning_learner_ukb_t300000, 
    "ukb_nikpay_t300000" = prsice_tuning_learner_ukb_nikpay_t300000),
  task = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task
  ),
  reg = reg_training_ukb_t300000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_ukb_t300000
)

batchtools::waitForJobs(reg = reg_training_ukb_t300000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_ukb_t300000, fun = function(res) list(res))

ukb_t300000_model <- models[job.id == 1, result[[1]]]
saveRDS(ukb_t300000_model, ukb_t300000_model_file)
ukb_t300000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(ukb_t300000_nikpay_model, ukb_t300000_nikpay_model_file)

## UKB - Only top 3000000 ----
reg_training_ukb_t3000000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_ukb_t3000000"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_ukb_t3000000 <- prsice_tuning_learner
prsice_tuning_learner_ukb_t3000000$opt.pars$pars$pval.level$lower <- 4.904758e-01
prsice_tuning_learner_ukb_t3000000$opt.pars$pars$pval.level$upper <- 4.904758e-01
prsice_tuning_learner_ukb_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_t3000000.rds")
prsice_tuning_learner_ukb_nikpay_t3000000 <- prsice_tuning_learner_ukb_t3000000
prsice_tuning_learner_ukb_nikpay_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_ukb_nikpay_t3000000.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "ukb_t3000000" = prsice_tuning_learner_ukb_t3000000, 
    "ukb_nikpay_t3000000" = prsice_tuning_learner_ukb_nikpay_t3000000),
  task = list(
    "ukb" = ukb_task,
    "ukb_nikpay" = ukb_nikpay_task
  ),
  reg = reg_training_ukb_t3000000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    name = sprintf("%s_Training", project_name),
    ntasks = 1,
    ncpus = tuning_cv_iters,
    partition = "prio",
    account = "GO-3269-1-1",
    walltime = 0L,
    memory = "100G",
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_ukb_t3000000
)

batchtools::waitForJobs(reg = reg_training_ukb_t3000000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_ukb_t3000000, fun = function(res) list(res))

ukb_t3000000_model <- models[job.id == 1, result[[1]]]
saveRDS(ukb_t3000000_model, ukb_t3000000_model_file)
ukb_t3000000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(ukb_t3000000_nikpay_model, ukb_t3000000_nikpay_model_file)

# Combined ----
reg_training_combined <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined"),
  work.dir = file.path(proc_dir),
  writeable = TRUE,
  overwrite = FALSE,
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)

prsice_tuning_learner_combined <- prsice_tuning_learner
prsice_tuning_learner_combined$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined.rds")
prsice_tuning_learner_combined_nikpay <- prsice_tuning_learner
prsice_tuning_learner_combined_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined" = prsice_tuning_learner_combined, 
    "combined_nikpay" = prsice_tuning_learner_combined_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 20000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined
)

batchtools::waitForJobs(reg = reg_training_combined)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined, fun = function(res) res)

combined_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_model, combined_model_file)
combined_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_nikpay_model, combined_nikpay_model_file)

## Combined - Only genome-wide significant -----
reg_training_combined_gws <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined_gws"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_combined_gws <- prsice_tuning_learner
prsice_tuning_learner_combined_gws$opt.pars$pars$pval.level$lower <- 5e-8
prsice_tuning_learner_combined_gws$opt.pars$pars$pval.level$upper <- 5e-8
prsice_tuning_learner_combined_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_gws.rds")
prsice_tuning_learner_combined_gws_nikpay <- prsice_tuning_learner_combined_gws
prsice_tuning_learner_combined_gws_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_gws_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined_gws" = prsice_tuning_learner_combined_gws, 
    "combined_gws_nikpay" = prsice_tuning_learner_combined_gws_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined_gws
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 70000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined_gws
)

batchtools::waitForJobs(reg = reg_training_combined_gws)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined_gws, fun = function(res) res)

combined_gws_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_gws_model, combined_gws_model_file)
combined_gws_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_gws_nikpay_model, combined_gws_nikpay_model_file)


## Combined - Only genome-wide significant -----
reg_training_combined_gws <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined_gws"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_combined_gws <- prsice_tuning_learner
prsice_tuning_learner_combined_gws$opt.pars$pars$pval.level$lower <- 5e-8
prsice_tuning_learner_combined_gws$opt.pars$pars$pval.level$upper <- 5e-8
prsice_tuning_learner_combined_gws$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_gws.rds")
prsice_tuning_learner_combined_gws_nikpay <- prsice_tuning_learner_combined_gws
prsice_tuning_learner_combined_gws_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_gws_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined_gws" = prsice_tuning_learner_combined_gws, 
    "combined_gws_nikpay" = prsice_tuning_learner_combined_gws_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined_gws
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 70000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined_gws
)

batchtools::waitForJobs(reg = reg_training_combined_gws)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined_gws, fun = function(res) res)

combined_gws_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_gws_model, combined_gws_model_file)
combined_gws_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_gws_nikpay_model, combined_gws_nikpay_model_file)

## Combined - Only top 3000 -----
reg_training_combined_t3000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined_t3000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_combined_t3000 <- prsice_tuning_learner
prsice_tuning_learner_combined_t3000$opt.pars$pars$pval.level$lower <- 2.54e-06
prsice_tuning_learner_combined_t3000$opt.pars$pars$pval.level$upper <- 2.54e-06
prsice_tuning_learner_combined_t3000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t3000.rds")
prsice_tuning_learner_combined_t3000_nikpay <- prsice_tuning_learner_combined_t3000
prsice_tuning_learner_combined_t3000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t3000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined_t3000" = prsice_tuning_learner_combined_t3000, 
    "combined_t3000_nikpay" = prsice_tuning_learner_combined_t3000_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined_t3000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 70000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined_t3000
)

batchtools::waitForJobs(reg = reg_training_combined_t3000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined_t3000, fun = function(res) res)

combined_t3000_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_t3000_model, combined_t3000_model_file)
combined_t3000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_t3000_nikpay_model, combined_t3000_nikpay_model_file)

## Combined - Only top 300000 -----
reg_training_combined_t300000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined_t300000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_combined_t300000 <- prsice_tuning_learner
prsice_tuning_learner_combined_t300000$opt.pars$pars$pval.level$lower <- 4.12907e-02
prsice_tuning_learner_combined_t300000$opt.pars$pars$pval.level$upper <- 4.12907e-02
prsice_tuning_learner_combined_t300000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t300000.rds")
prsice_tuning_learner_combined_t300000_nikpay <- prsice_tuning_learner_combined_t300000
prsice_tuning_learner_combined_t300000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t300000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined_t300000" = prsice_tuning_learner_combined_t300000, 
    "combined_t300000_nikpay" = prsice_tuning_learner_combined_t300000_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined_t300000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 70000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined_t300000
)

batchtools::waitForJobs(reg = reg_training_combined_t300000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined_t300000, fun = function(res) res)

combined_t300000_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_t300000_model, combined_t300000_model_file)
combined_t300000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_t300000_nikpay_model, combined_t300000_nikpay_model_file)

## Combined - Only top 3000000 -----
reg_training_combined_t3000000 <- imbs::load_or_create_registry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "training_combined_t3000000"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R")),
  writeable = TRUE,
  overwrite = FALSE
)

prsice_tuning_learner$control$mbo.control$save.on.disk.at.time <- 1
prsice_tuning_learner$control$mbo.control$on.surrogate.error <- "warn"
prsice_tuning_learner$control$continue <- TRUE

prsice_tuning_learner_combined_t3000000 <- prsice_tuning_learner
prsice_tuning_learner_combined_t3000000$opt.pars$pars$pval.level$lower <- 4.904758e-01
prsice_tuning_learner_combined_t3000000$opt.pars$pars$pval.level$upper <- 4.904758e-01
prsice_tuning_learner_combined_t3000000$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t3000000.rds")
prsice_tuning_learner_combined_t3000000_nikpay <- prsice_tuning_learner_combined_t3000000
prsice_tuning_learner_combined_t3000000_nikpay$control$mbo.control$save.file.path <- file.path(proc_dir, "prsice_mbo_tuning_combined_t3000000_nikpay.rds")

ids <- batchtools::batchMap(
  fun = mlr::train,
  learner = list(
    "combined_t3000000" = prsice_tuning_learner_combined_t3000000, 
    "combined_t3000000_nikpay" = prsice_tuning_learner_combined_t3000000_nikpay
  ),
  task = list(
    "combined" = combined_task,
    "combined_nikpay" = combined_nikpay_task
  ),
  reg = reg_training_combined_t3000000
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(
    ntasks = 3,
    name = sprintf("%s_Training_Combined", project_name),
    ncpus = tuning_cv_iters,
    partition = "batch",
    walltime = 0L,
    memory = 70000,
    chunks.as.arrayjobs = TRUE,
    pm.backend = "multicore",
    pm.opts = list(
      level = "mlr.resample"
    )
  ),
  reg = reg_training_combined_t3000000
)

batchtools::waitForJobs(reg = reg_training_combined_t3000000)

models <- batchtools::reduceResultsDataTable(reg = reg_training_combined_t3000000, fun = function(res) res)

combined_t3000000_model <- models[job.id == 1, result[[1]]]
saveRDS(combined_t3000000_model, combined_t3000000_model_file)
combined_t3000000_nikpay_model <- models[job.id == 2, result[[1]]]
saveRDS(combined_t3000000_nikpay_model, combined_t3000000_nikpay_model_file)



