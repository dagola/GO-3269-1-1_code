
source("init.R")

# Define measures ----
measures <- list(aucpr, # use Area under Precision-Recall Curve as performance measure during tuning
                 mlr::mmce, mlr::npv, mlr::fpr, mlr::f1, mlr::fnr, mlr::ssr,
                 mlr::tp, mlr::tn, mlr::gpr, mlr::lsr, mlr::acc,
                 mlr::wkappa, mlr::ppv, mlr::logloss, mlr::ber, mlr::tpr,
                 mlr::brier, mlr::gmean, mlr::fdr, mlr::tnr, mlr::qsr, mlr::bac,
                 mlr::brier.scaled, mlr::fp, mlr::fn, mlr::kappa, mlr::auc)

# Load tasks ----
ukb_testing_task <- readRDS(ukb_testing_task_file)
de_testing_task <- readRDS(de_testing_task_file)
de_downsampled_testing_task <- readRDS(de_downsampled_testing_task_file)
g6_testing_task <- readRDS(g6_testing_task_file)
g7_testing_task <- readRDS(g7_testing_task_file)
eb_testing_task <- readRDS(eb_testing_task_file)

# Load models ----
ukb_model <- readRDS(ukb_model_file)
ukb_gws_model <- readRDS(ukb_gws_model_file)
ukb_t3000_model <- readRDS(ukb_t3000_model_file)
ukb_t3000_model <- readRDS(ukb_t300000_model_file)
ukb_t3000_model <- readRDS(ukb_t3000000_model_file)
ukb_nikpay_model <- readRDS(ukb_nikpay_model_file)
ukb_t3000_nikpay_model <- readRDS(ukb_t3000_nikpay_model_file)
ukb_t300000_nikpay_model <- readRDS(ukb_t300000_nikpay_model_file)
ukb_t3000000_nikpay_model <- readRDS(ukb_t3000000_nikpay_model_file)
ukb_gws_nikpay_model <- readRDS(ukb_gws_nikpay_model_file)
eb_model <- readRDS(eb_model_file)
eb_gws_model <- readRDS(eb_gws_model_file)
eb_t3000_model <- readRDS(eb_t3000_model_file)
eb_t300000_model <- readRDS(eb_t300000_model_file)
eb_t3000000_model <- readRDS(eb_t3000000_model_file)
eb_nikpay_model <- readRDS(eb_nikpay_model_file)
eb_gws_nikpay_model <- readRDS(eb_gws_nikpay_model_file)
eb_t3000_nikpay_model <- readRDS(eb_t3000_nikpay_model_file)
eb_t300000_nikpay_model <- readRDS(eb_t300000_nikpay_model_file)
eb_t3000000_nikpay_model <- readRDS(eb_t3000000_nikpay_model_file)
de_model <- readRDS(de_model_file)
de_gws_model <- readRDS(de_gws_model_file)
de_t3000_model <- readRDS(de_t3000_model_file)
de_t300000_model <- readRDS(de_t300000_model_file)
de_t3000000_model <- readRDS(de_t3000000_model_file)
de_nikpay_model <- readRDS(de_nikpay_model_file)
de_gws_nikpay_model <- readRDS(de_gws_nikpay_model_file)
de_t3000_nikpay_model <- readRDS(de_t3000_nikpay_model_file)
de_t300000_nikpay_model <- readRDS(de_t300000_nikpay_model_file)
de_t3000000_nikpay_model <- readRDS(de_t3000000_nikpay_model_file)
de_downsampled_model <- readRDS(de_downsampled_model_file)
de_gws_downsampled_model <- readRDS(de_gws_downsampled_model_file)
de_t3000_downsampled_model <- readRDS(de_t3000_downsampled_model_file)
de_t300000_downsampled_model <- readRDS(de_t300000_downsampled_model_file)
de_t3000000_downsampled_model <- readRDS(de_t3000000_downsampled_model_file)
de_downsampled_nikpay_model <- readRDS(de_downsampled_nikpay_model_file)
de_gws_downsampled_nikpay_model <- readRDS(de_gws_downsampled_nikpay_model_file)
de_t3000_downsampled_nikpay_model <- readRDS(de_t3000_downsampled_nikpay_model_file)
de_t300000_downsampled_nikpay_model <- readRDS(de_t300000_downsampled_nikpay_model_file)
de_t3000000_downsampled_nikpay_model <- readRDS(de_t3000000_downsampled_nikpay_model_file)
combined_model <- readRDS(combined_model_file)
combined_gws_model <- readRDS(combined_gws_model_file)
combined_t3000_model <- readRDS(combined_t3000_model_file)
combined_t300000_model <- readRDS(combined_t300000_model_file)
combined_t3000000_model <- readRDS(combined_t3000000_model_file)
combined_nikpay_model <- readRDS(combined_nikpay_model_file)
combined_gws_nikpay_model <- readRDS(combined_gws_nikpay_model_file)
combined_t3000_nikpay_model <- readRDS(combined_t3000_nikpay_model_file)
combined_t300000_nikpay_model <- readRDS(combined_t300000_nikpay_model_file)
combined_t3000000_nikpay_model <- readRDS(combined_t3000000_nikpay_model_file)

khera_model <- readRDS(khera_model_file)

# Set system dependent learner pars ----
ukb_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_gws_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_gws_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_gws_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_gws_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t3000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t3000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t300000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t300000_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t300000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t300000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t3000000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000000_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t3000000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_gws_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_gws_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_gws_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_gws_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t3000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t3000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t300000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t300000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t300000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t300000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

ukb_t3000000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
ukb_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
ukb_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_gws_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_gws_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_gws_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_gws_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t3000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t3000_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t3000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t3000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t300000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t300000_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t300000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t300000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t3000000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t3000000_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t3000000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t3000000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_gws_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_gws_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_gws_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_gws_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t3000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t3000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t3000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t3000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t300000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t300000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t300000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t300000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

eb_t3000000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
eb_t3000000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
eb_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
eb_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_gws_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_gws_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_gws_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_gws_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t300000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t300000_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t300000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t300000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_gws_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_gws_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_gws_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_gws_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t300000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t300000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t300000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t300000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_downsampled_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_downsampled_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_downsampled_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_downsampled_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_gws_downsampled_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_gws_downsampled_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_gws_downsampled_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_gws_downsampled_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000_downsampled_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000_downsampled_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000_downsampled_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000_downsampled_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t300000_downsampled_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t300000_downsampled_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t300000_downsampled_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t300000_downsampled_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000000_downsampled_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_downsampled_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000000_downsampled_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_downsampled_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_downsampled_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_downsampled_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_gws_downsampled_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_gws_downsampled_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_gws_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_gws_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000_downsampled_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000_downsampled_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t300000_downsampled_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t300000_downsampled_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t300000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t300000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

de_t3000000_downsampled_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_downsampled_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
de_t3000000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
de_t3000000_downsampled_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_gws_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_gws_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_gws_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_gws_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t3000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t3000_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t3000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t3000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t300000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t300000_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t300000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t300000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t3000000_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t3000000_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t3000000_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t3000000_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_gws_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_gws_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_gws_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_gws_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t3000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t3000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t3000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t3000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t300000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t300000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t300000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t300000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

combined_t3000000_nikpay_model$learner$next.learner$par.vals$prsice.exec <- prsice_exec
combined_t3000000_nikpay_model$learner$next.learner$par.vals$scratch.dir <- NULL
combined_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$prsice.exec <- prsice_exec
combined_t3000000_nikpay_model$learner.model$next.model$learner$par.vals$scratch.dir <- NULL

khera_model$learner$par.vals$prsice.exec <- prsice_exec
khera_model$learner$par.vals$scratch.dir <- NULL

# Predictions ----
reg_testing <- batchtools::makeExperimentRegistry(
  packages = c("data.table", "mlr", "checkmate", "imbs"),
  file.dir = file.path(registries_dir, "testing"),
  work.dir = file.path(proc_dir), 
  source = file.path(functions_dir, c("prsice.R", "learner.R", "plink2.R", "aucpr.R"))
)
reg_testing <- batchtools::loadRegistry(
  file.dir = file.path(registries_dir, "testing"), 
  writeable = TRUE
)

batchtools::batchExport(
  export = list(
    ukb_model = ukb_model,
    ukb_gws_model = ukb_gws_model,
    ukb_t3000_model = ukb_t3000_model,
    ukb_t300000_model = ukb_t300000_model,
    ukb_t3000000_model = ukb_t3000000_model,
    ukb_nikpay_model = ukb_nikpay_model,
    ukb_gws_nikpay_model = ukb_gws_nikpay_model,
    ukb_t3000_nikpay_model = ukb_t3000_nikpay_model,
    ukb_t300000_nikpay_model = ukb_t300000_nikpay_model,
    ukb_t3000000_nikpay_model = ukb_t3000000_nikpay_model,
    eb_model = eb_model,
    eb_gws_model = eb_gws_model,
    eb_t3000_model = eb_t3000_model,
    eb_t300000_model = eb_t300000_model,
    eb_t3000000_model = eb_t3000000_model,
    eb_nikpay_model = eb_nikpay_model,
    eb_gws_nikpay_model = eb_gws_nikpay_model,
    eb_t3000_nikpay_model = eb_t3000_nikpay_model,
    eb_t300000_nikpay_model = eb_t300000_nikpay_model,
    eb_t3000000_nikpay_model = eb_t3000000_nikpay_model,
    de_model = de_model,
    de_gws_model = de_gws_model,
    de_t3000_model = de_t3000_model,
    de_t300000_model = de_t300000_model,
    de_t3000000_model = de_t3000000_model,
    de_nikpay_model = de_nikpay_model,
    de_gws_nikpay_model = de_gws_nikpay_model,
    de_t3000_nikpay_model = de_t3000_nikpay_model,
    de_t300000_nikpay_model = de_t300000_nikpay_model,
    de_t3000000_nikpay_model = de_t3000000_nikpay_model,
    de_downsampled_model = de_downsampled_model,
    de_gws_downsampled_model = de_gws_downsampled_model,
    de_t3000_downsampled_model = de_t3000_downsampled_model,
    de_t300000_downsampled_model = de_t300000_downsampled_model,
    de_t3000000_downsampled_model = de_t3000000_downsampled_model,
    de_downsampled_nikpay_model = de_downsampled_nikpay_model,
    de_gws_downsampled_nikpay_model = de_gws_downsampled_nikpay_model,
    de_t3000_downsampled_nikpay_model = de_t3000_downsampled_nikpay_model,
    de_t300000_downsampled_nikpay_model = de_t300000_downsampled_nikpay_model,
    de_t3000000_downsampled_nikpay_model = de_t3000000_downsampled_nikpay_model,
    combined_model = combined_model,
    combined_gws_model = combined_gws_model,
    combined_t3000_model = combined_t3000_model,
    combined_t300000_model = combined_t300000_model,
    combined_t3000000_model = combined_t3000000_model,
    combined_nikpay_model = combined_nikpay_model,
    combined_gws_nikpay_model = combined_gws_nikpay_model,
    combined_t3000_nikpay_model = combined_t3000_nikpay_model,
    combined_t300000_nikpay_model = combined_t300000_nikpay_model,
    combined_t3000000_nikpay_model = combined_t3000000_nikpay_model,
    khera_model = khera_model
  ),
  reg = reg_testing
)

## Add problems ----
batchtools::addProblem(
  name = "UKB_TESTING",
  data = ukb_testing_task,
  reg = reg_testing
)

batchtools::addProblem(
  name = "DE_TESTING",
  data = de_testing_task,
  reg = reg_testing
)

batchtools::addProblem(
  name = "G6_TESTING",
  data = g6_testing_task,
  reg = reg_testing
)

batchtools::addProblem(
  name = "G7_TESTING",
  data = g7_testing_task,
  reg = reg_testing
)

batchtools::addProblem(
  name = "DE_DOWNSAMPLED_TESTING",
  data = de_downsampled_testing_task,
  reg = reg_testing
)

batchtools::addProblem(
  name = "EB_TESTING",
  data = eb_testing_task,
  reg = reg_testing
)

## Add algorithms ----

batchtools::addAlgorithm(
  name = "UKB",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_GWS",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_gws_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T3000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t3000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T300000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t300000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T3000000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t3000000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_GWS_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_gws_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T3000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t3000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T300000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t300000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "UKB_T3000000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = ukb_t3000000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_GWS",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_gws_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T3000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t3000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T300000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t300000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T3000000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t3000000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T3000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t3000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T300000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t300000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "EB_T3000000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = eb_t3000000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_GWS",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_gws_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T300000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t300000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_GWS_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_gws_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T300000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t300000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_DOWNSAMPLED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_downsampled_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_GWS_DOWNSAMPLED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_gws_downsampled_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000_DOWNSAMPLED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000_downsampled_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T300000_DOWNSAMPLED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t300000_downsampled_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000000_DOWNSAMPLED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000000_downsampled_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_DOWNSAMPLED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_downsampled_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_GWS_DOWNSAMPLED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_gws_downsampled_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000_DOWNSAMPLED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000_downsampled_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T300000_DOWNSAMPLED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t300000_downsampled_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "DE_T3000000_DOWNSAMPLED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = de_t3000000_downsampled_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_GWS",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_gws_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T3000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t3000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T300000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t300000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T3000000",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t3000000_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_GWS_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_gws_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T3000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t3000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T300000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t300000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "COMBINED_T3000000_NIKPAY",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = combined_t3000000_nikpay_model)
  },
  reg = reg_testing
)

batchtools::addAlgorithm(
  name = "KHERA",
  fun = function(job, data, instance, ...) {
    mlr:::predict.WrappedModel(task = data, object = khera_model)
  },
  reg = reg_testing
)

batchtools::addExperiments(reg = reg_testing)

batch_ids <- findNotDone(reg = reg_testing)
batch_ids <- ajoin(batch_ids, findRunning(reg = reg_testing))
batch_ids[, chunk := 1]
batchtools::submitJobs(
  ids = batch_ids,
  resources = list(
    name = sprintf("%s_Testing", project_name),
    ntasks = 1,
    ncpus = 1,
    partition = "batch,prio",
    account = project_name,
    walltime = 0L,
    memory = 90000,
    chunks.as.arrayjobs = TRUE
  ),
  reg = reg_testing
)

batchtools::waitForJobs(reg = reg_testing)

testing_performance <- batchtools::ijoin(
  x = unwrap(getJobPars(reg = reg_testing)), 
  y = unwrap(
    batchtools::reduceResultsDataTable(
      reg = reg_testing, 
      fun = function(pred) {
        pr_roc <- PRROC::pr.curve(
          scores.class0 = pred$data$prob.CAD,
          weights.class0 = as.numeric(pred$data$truth) - 1, 
          max.compute = TRUE, 
          min.compute = TRUE, 
          rand.compute = TRUE
        )
        performances <- mlr::performance(
          pred = pred, 
          measures = measures
        )
        aucpr_ci <- ci.aucpr(aucpr = performances["aucpr"], n = pred$task.desc$class.distribution[pred$task.desc$positive], ci.level = 0.95)
        auc_ci <- tryCatch(
          expr = {
            structure(
              pROC::ci.auc(
                formula = as.formula(sprintf("truth ~ prob.%s", pred$task.desc$positive)), 
                data = pred$data,
                conf.level = 0.95,
                method = "delong"
              )[c(1, 3)],
              names = c("auc.lower", "auc.upper")
            )
          },
          error = function(e) {
            structure(
              c(NA_real_, NA_real_),
              names = c("auc.lower", "auc.upper")
            )
          }
        )
        as.list(
          c(
            performances,
            auc_ci,
            aucpr_ci,
            aucpr.min = pr_roc$min$auc.davis.goadrich,
            aucpr.max = pr_roc$max$auc.davis.goadrich,
            aucpr.rand = pr_roc$rand$auc.davis.goadrich
          )
        )
      }
    )
  )
)
setorder(testing_performance, -aucpr, na.last = TRUE)

saveRDS(testing_performance, file = testing_performance_file)

testing_results <- batchtools::ijoin(
  x = unwrap(batchtools::getJobPars(reg = reg_testing)), 
  y = rbindlist(
    l = batchtools::reduceResultsList(
      reg = reg_testing, 
      fun = function(res, job) {
        cbind("job.id" = job$job.id, res$data)
      }
    )
  ),
  by = "job.id"
)
testing_results[, id := as.character(id)]

saveRDS(testing_results, file = testing_results_file)
