
source("init.R")

# Define Khera model ----
de_nikpay_model <- readRDS(de_nikpay_model_file)
khera_model <- de_nikpay_model$learner.model$next.model

khera_scores <- data.table::fread(khera_weights_file)

khera_model$learner.model$model <- khera_scores[, .(
  ID = sprintf("%s:%s", chr, position_hg19), 
  "#CHROM" = chr, 
  POS = position_hg19, 
  A1 = effect_allele, 
  AX = A2,
  BETA = effect_weight, 
  P = 1, 
  A1_FREQ = 0, 
  N = 0
)]
khera_model$learner$par.vals$missing.handling <- "SET_ZERO"
khera_model$learner$par.vals$model <- "add"
khera_model$learner$par.vals$score <- "sum"
khera_model$learner$par.vals$prsice.exec <- prsice_exec
khera_model$learner$par.vals$scratch.dir <- NULL

saveRDS(object = khera_model, file = khera_model_file)
