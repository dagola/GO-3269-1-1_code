
create_prsice_model <- function(
  base.file,
  summary.statistics.chr.col, 
  summary.statistics.pos.col, 
  summary.statistics.snp.id.col, 
  summary.statistics.effect.allele.col,
  summary.statistics.non.effect.allele.col,
  summary.statistics.stat.col,
  summary.statistics.pvalue.col,
  summary.statistics.info.col,
  summary.statistics.info.threshold,
  summary.statistics.maf.cols,
  summary.statistics.maf.thresholds,
  summary.statistics.beta,
  summary.statistics.or,
  summary.statistics.index,
  target.file.prefix,
  target.keep,
  target.type,
  target.geno,
  target.maf,
  target.info,
  target.nonfounders,
  target.hardcall,
  target.dosage.threshold,
  target.hardcall.threshold,
  ld.external,
  ld.file.prefix,
  ld.type,
  ld.geno,
  ld.maf,
  ld.info,
  ld.keep,
  ld.dosage.threshold,
  ld.hardcall.threshold,
  clumping,
  clumping.kb,
  clumping.r2,
  clumping.p,
  clumping.proxy,
  clumping.pearson,
  dosage.allow.intermediate,
  pval.level,
  model,
  missing.handling,
  score,
  exclude,
  exclude.range,
  extract,
  id.delim,
  ignore.fid,
  keep.ambig,
  scratch.dir,
  scratch.dir.files,
  memory,
  num.threads,
  prsice.exec,
  ...
) {
  
  assertions <- checkmate::makeAssertCollection()
  
  args <- list()
  
  if (!missing(base.file)) {
    checkmate::assert_file(base.file, access = "r", add = assertions)
    args$base.file <- base.file
  }
  
  output_prefix <- tempfile()
  on.exit(BBmisc::suppressAll(file.remove(sprintf(c("%s.all.score", "%s.snp", "%s.log"), output_prefix))), add = TRUE)
  
  args$output.file.prefix <- output_prefix
  
  base_options <- list(
    chr.col = "#CHROM",
    pos.col = "POS",
    snp.id.col = "ID",
    effect.allele.col = "A1",
    pvalue.col = "P"
  )
  if (!missing(summary.statistics.index)) {
    checkmate::assert_flag(summary.statistics.index, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$index <- summary.statistics.index
  }
  if (!missing(summary.statistics.chr.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.chr.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.chr.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$chr.col <- summary.statistics.chr.col
  }
  if (!missing(summary.statistics.pos.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.pos.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.pos.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$pos.col <- summary.statistics.pos.col
  }
  if (!missing(summary.statistics.snp.id.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.snp.id.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.snp.id.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$snp.id.col <- summary.statistics.snp.id.col
  }
  if (!missing(summary.statistics.effect.allele.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.effect.allele.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.effect.allele.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$effect.allele.col <- summary.statistics.effect.allele.col
  }
  if (!missing(summary.statistics.non.effect.allele.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.non.effect.allele.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.non.effect.allele.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$non.effect.allele.col <- summary.statistics.non.effect.allele.col
  }
  if (!missing(summary.statistics.stat.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.stat.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.stat.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$stat.col <- summary.statistics.stat.col
  }
  if (!missing(summary.statistics.pvalue.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.pvalue.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.pvalue.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$pvalue.col <- summary.statistics.pvalue.col
  }
  if (!missing(summary.statistics.info.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.info.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.info.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$info.col <- summary.statistics.info.col
  }
  if (!missing(summary.statistics.info.threshold)) {
    checkmate::assert_number(summary.statistics.info.threshold, na.ok = FALSE, null.ok = FALSE, lower = 0, upper = 1, finite = TRUE, add = assertions)
    base_options$info.threshold <- summary.statistics.info.threshold
  }
  if (!missing(summary.statistics.maf.cols)) {
    if (summary.statistics.index) {
      checkmate::assert_integer(summary.statistics.maf.cols, any.missing = FALSE, all.missing = FALSE, min.len = 1, unique = TRUE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.maf.cols, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$maf.cols <- summary.statistics.maf.cols
  }
  if (!missing(summary.statistics.maf.thresholds)) {
    checkmate::assert_numeric(summary.statistics.maf.thresholds, null.ok = FALSE, any.missing = FALSE, all.missing = FALSE, len = length(summary.statistics.maf.cols), lower = 0, upper = 1, finite = TRUE, add = assertions)
    base_options$maf.thresholds <- summary.statistics.maf.thresholds
  }
  if (!missing(summary.statistics.beta)) {
    checkmate::assert_flag(summary.statistics.beta, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$beta <- summary.statistics.beta
  }
  if (!missing(summary.statistics.or)) {
    checkmate::assert_flag(summary.statistics.or, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$or <- summary.statistics.or
  }
  args$base.options <- base_options
  
  if (!missing(target.type)) {
    checkmate::assert_choice(target.type, choices = c("bed", "bgen"), null.ok = FALSE, add = assertions)
  } else {
    target.type <- "bed"
  }
  
  if (target.type == "bed") {
    checkmate::assert_file(sprintf("%s.bed", target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", target.file.prefix), access = "r", add = assertions)
  } else {
    checkmate::assert_file(sprintf("%s.bgen", target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bgen.bgi", target.file.prefix), access = "r", add = assertions)
  }
  
  args$target.file.prefix <- target.file.prefix
  target_options <- list()
  target_options$type <- target.type
  if (!missing(target.geno)) {
    checkmate::assert_number(target.geno, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    target_options$geno <- target.geno
    target_geno <- sprintf("--geno %f", target.geno)
  } else {
    target_geno <- ""
  }
  if (!missing(target.maf)) {
    checkmate::assert_number(target.maf, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    target_options$maf <- target.maf
    target_maf <- sprintf("--maf %f", target.maf)
  } else {
    target_maf <- ""
  }
  if (!missing(target.info)) {
    checkmate::assert_number(target.info, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    target_options$info <- target.info
  }
  if (!missing(target.nonfounders)) {
    checkmate::assert_flag(target.nonfounders, na.ok = FALSE, null.ok = FALSE, add = assertions)
    target_options$nonfounders <- target.nonfounders
  }
  if (!missing(target.keep)) {
    checkmate::assert_file(target.keep, access = "r", add = assertions)
    target_options$keep <- target.keep
    target_keep <- sprintf("--keep %s", target.keep)
    pheno_file <- sprintf("--pheno %s", target.keep)
  } else {
    target_keep <- ""
    pheno_file <- ""
  }
  args$target.options <- target_options
  
  ld_file_options <- list()
  if (ld.external) {
    if (!missing(ld.type)) {
      checkmate::assert_choice(ld.type, choices = c("bed", "bgen"), null.ok = FALSE, add = assertions)
    } else {
      ld.type <- "bed"
    }
    if (!missing(ld.file.prefix)) {
      if (ld.type == "bed") {
        checkmate::assert_file(sprintf("%s.bed", ld.file.prefix), access = "r", add = assertions)
        checkmate::assert_file(sprintf("%s.bim", ld.file.prefix), access = "r", add = assertions)
        checkmate::assert_file(sprintf("%s.fam", ld.file.prefix), access = "r", add = assertions)
      } else {
        checkmate::assert_file(sprintf("%s.bgen", ld.file.prefix), access = "r", add = assertions)
        checkmate::assert_file(sprintf("%s.bgen.bgi", ld.file.prefix), access = "r", add = assertions)
      }
    }
  }
  if (!missing(ld.geno)) {
    checkmate::assert_number(ld.geno, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    ld_file_options$geno <- ld.geno
  }
  if (!missing(ld.maf)) {
    checkmate::assert_number(ld.maf, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    ld_file_options$maf <- ld.maf
  }
  if (!missing(ld.info)) {
    checkmate::assert_number(ld.info, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    ld_file_options$info <- ld.info
  }
  if (!missing(ld.keep)) {
    checkmate::assert_file(ld.keep, access = "r", add = assertions)
    ld_file_options$keep <- ld.keep
  }
  if (!missing(ld.dosage.threshold)) {
    checkmate::assert_number(ld.dosage.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    ld_file_options$dosage.threshold <- ld.dosage.threshold
  }
  if (!missing(ld.hardcall.threshold)) {
    checkmate::assert_number(ld.hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    ld_file_options$hardcall.threshold <- ld.hardcall.threshold
  }
  args$ld.file.options <- ld_file_options
  
  if (!missing(clumping)) {
    checkmate::assert_flag(clumping, na.ok = FALSE, null.ok = FALSE, add = assertions)
    args$clumping <- clumping
  }
  clumping_options <- list()
  if (!missing(clumping.kb)) {
    checkmate::assert_int(clumping.kb, na.ok = FALSE, lower = 1, null.ok = FALSE, add = assertions)
    clumping_options$kb <- clumping.kb
  }
  if (!missing(clumping.r2)) {
    checkmate::assert_number(clumping.r2, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    clumping_options$r2 <- clumping.r2
  }
  if (!missing(clumping.p)) {
    checkmate::assert_number(clumping.p, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    clumping_options$p <- clumping.p
  }
  if (!missing(clumping.proxy)) {
    checkmate::assert_number(clumping.proxy, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    clumping_options$proxy <- clumping.proxy
  }
  if (!missing(clumping.pearson)) {
    checkmate::assert_flag(clumping.pearson, na.ok = FALSE, null.ok = FALSE, add = assertions)
    clumping_options$pearson <- clumping.pearson
  }
  args$clumping.options <- clumping_options
  
  dosage_options <- list()
  if (!missing(target.hardcall)) {
    checkmate::assert_flag(target.hardcall, na.ok = FALSE, null.ok = FALSE, add = assertions)
    dosage_options$hardcall <- target.hardcall
  }
  if (!missing(target.dosage.threshold)) {
    checkmate::assert_number(target.dosage.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    dosage_options$dosage.threshold <- target.dosage.threshold
  }
  if (!missing(target.hardcall.threshold)) {
    checkmate::assert_number(target.hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    dosage_options$hardcall.threshold <- target.hardcall.threshold
  }
  if (!missing(dosage.allow.intermediate)) {
    checkmate::assert_flag(dosage.options$allow.intermediate, na.ok = FALSE, null.ok = FALSE, add = assertions)
    dosage_options$allow.intermediate <- dosage.allow.intermediate
  }
  args$dosage.options <- dosage_options
  
  checkmate::assert_number(pval.level, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  pval_options <- list(
    levels = pval.level,
    full = FALSE
  )
  if (!missing(model)) {
    checkmate::assert_choice(model, choices = c("add", "dom", "rec", "het"), null.ok = FALSE, add = assertions)
    pval_options$model <- model
  }
  args$pval.options <- pval_options
  
  if (!missing(missing.handling)) {
    checkmate::assert_choice(missing.handling, choices = c("IMPUTE", "SET_ZERO", "CENTER"), null.ok = FALSE, add = assertions)
    args$missing.handling <- missing.handling
  }
  
  args$regress <- FALSE
  
  if (!missing(score)) {
    checkmate::assert_choice(score, choices = c("avg", "std", "sum"), null.ok = FALSE, add = assertions)
    args$score <- score
  }
  if (!missing(exclude)) {
    checkmate::assert_file(exclude, access = "r", add = assertions)
    args$exclude <- exclude
    exclude <- sprintf("--exclude %s", exclude)
  } else {
    exclude <- ""
  }
  if (!missing(exclude.range)) {
    checkmate::assert_data_frame(exclude.range, types = c("integer", "numeric"), any.missing = FALSE, all.missing = FALSE, min.rows = 1, min.cols = 3, col.names = "unique", null.ok = FALSE, add = assertions)
    checkmate::assert_set_equal(colnames(exclude.range), c("CHR", "START", "END"), add = assertions)
    args$exclude.range <- exclude.range
  }
  if (!missing(extract)) {
    checkmate::assert_file(extract, access = "r", add = assertions)
    args$extract <- extract
    extract <- sprintf("--extract %s", extract)
  } else {
    extract <- ""
  }
  if (!missing(id.delim)) {
    checkmate::assert_string(id.delim, na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
    args$id.delim <- id.delim
  }
  if (!missing(ignore.fid)) {
    checkmate::assert_flag(ignore.fid, na.ok = FALSE, null.ok = FALSE, add = assertions)
    args$ignore.fid <- ignore.fid
  }
  if (!missing(keep.ambig)) {
    checkmate::assert_flag(keep.ambig, na.ok = FALSE, null.ok = FALSE, add = assertions)
    args$keep.ambig <- keep.ambig
  }
  
  copy_to_scratch_dir <- FALSE
  if (!missing(scratch.dir)) {
    checkmate::assert_directory(scratch.dir, access = "rw", add = assertions)
    copy_to_scratch_dir <- TRUE
  }
  
  if (!missing(scratch.dir.files)) {
    checkmate::assert_character(scratch.dir.files, add = assertions)
    checkmate::assert_subset(scratch.dir.files, choices = c("base.file", "target.file", "ld.file"), add = assertions)
  } else {
    scratch.dir.files <- character()
  }
  
  args$print.snps <- TRUE
  
  args$seed <- sample(.Machine$integer.max, 1)
  
  if (!missing(memory)) {
    checkmate::assert_int(memory, lower = 1000, add = assertions)
    args$memory <- memory
  }
  if (!missing(num.threads)) {
    checkmate::assert_int(num.threads, lower = 1, add = assertions)
    args$num.threads <- num.threads
  }
  checkmate::assert_string(prsice.exec, na.ok = FALSE, null.ok = FALSE, add = assertions)
  prsice.exec <- path.expand(prsice.exec)
  imbs::assert_command(prsice.exec, add = assertions)
  args$exec <- prsice.exec
  
  checkmate::reportAssertions(assertions)
  
  # Copy files to scratch.dir ----
  if (copy_to_scratch_dir) {
    file_list <- character()
    if (!missing(base.file)) {
      if ("base.file" %in% scratch.dir.files) {
        file_list <- c(file_list, base.file)
        args$base.file <- file.path(scratch.dir, basename(base.file))
      }
    }
    if ("target.file" %in% scratch.dir.files) {
      file_list <- c(file_list, sprintf("%s*", target.file.prefix))
      args$target.file.prefix <- file.path(scratch.dir, basename(target.file.prefix))
    }
    if (!missing(ld.file.prefix)) {
      if ("ld.file" %in% scratch.dir.files) {
        file_list <- c(file_list, sprintf("%s*", ld.file.prefix))
        args$ld.file.prefix <- file.path(scratch.dir, basename(ld.file.prefix))
      }
    }
    if (length(file_list) > 1) {
      file_list <- sprintf("{%s}", paste(file_list, collapse = ","))
    } else {
      file_list <- sprintf("%s", paste(file_list, collapse = " "))
    }
    if (imbs::test_command("rsync")) {
      imbs::system_call(bin = "rsync", args = c("-avPzh", file_list, scratch.dir))
    } else {
      imbs::system_call(bin = "cp", args = c(file_list, scratch.dir))
    }
  }
  
  if (missing(base.file)) {
    if (!missing(scratch.dir)) {
      summary_statistics_file <- tempfile(tmpdir = scratch.dir)
    } else {
      summary_statistics_file <- tempfile()
    }
    on.exit(file.remove(sprintf("%s.log", summary_statistics_file)), add = TRUE)
    on.exit(file.remove(sprintf("%s.PHENO1.glm.logistic.hybrid", summary_statistics_file)), add = TRUE)
    
    plink2_gwa(
      input.prefix = args$target.file.prefix, 
      output.prefix = summary_statistics_file, 
      pheno.file = pheno_file,
      pheno.cols = "--pheno-col-nums 3",
      hide.covar = TRUE,
      firth = "firth-fallback",
      no.fam.pheno = "--no-pheno",
      maf = target_maf,
      geno = target_geno,
      extract = extract,
      exclude = exclude,
      keep = target_keep,
      pheno.coding = "--1",
      num.threads = num.threads,
      memory = memory,
      output.cols = c("chrom", "pos", "ref", "alt", "ax", "a1freq", "a1freqcc", "firth", "beta", "se", "tz", "p"),
      ...
    )
    base.file <- sprintf("%s.PHENO1.glm.logistic.hybrid", summary_statistics_file)
    checkmate::assert_file(base.file, access = "r")
    args$base.file <- base.file
  }
  
  do.call(prsice, args)
  
  model_snps <- data.table::fread(file = sprintf("%s.snp", output_prefix))
  summary_statistics <- data.table::fread(file = args$base.file)
  model <- merge(summary_statistics, model_snps[, "SNP"], by.x = base_options$snp.id.col, by.y = "SNP", all = FALSE)
  
  if (summary.statistics.or) {
    model[, BETA := log(get(summary.statistics.stat.col))]
  } else {
    model[, BETA := get(summary.statistics.stat.col)]
  }
  
  max_score <- model[BETA > 0, sum(2*BETA)]
  min_score <- model[BETA < 0, sum(2*BETA)]
  
  output <- list(
    model = model,
    min.score = min_score,
    max.score = max_score
  )
  class(output) <- "prsice"
  
  return(output)
  
}

predict.prsice <- function(
  object,
  new.target.file.prefix, 
  summary.statistics.chr.col,
  summary.statistics.pos.col, 
  summary.statistics.snp.id.col,
  summary.statistics.effect.allele.col, 
  summary.statistics.non.effect.allele.col, 
  summary.statistics.stat.col,
  summary.statistics.pvalue.col, 
  summary.statistics.beta,
  summary.statistics.or, 
  summary.statistics.index,
  target.type,
  target.nonfounders, 
  target.keep, 
  target.hardcall, 
  target.dosage.threshold,
  target.hardcall.threshold, 
  dosage.allow.intermediate,
  score,
  model,
  missing.handling,
  id.delim,
  ignore.fid, 
  keep.ambig, 
  scratch.dir, 
  scratch.dir.files,
  print.snps, 
  prsice.exec,
  memory, 
  num.threads,
  ...
) {
  
  assertions <- checkmate::makeAssertCollection()
  
  args <- list()
  
  model_file <- tempfile()
  on.exit(BBmisc::suppressAll(file.remove(model_file)), add = TRUE)
  output_prefix <- tempfile()
  on.exit(BBmisc::suppressAll(file.remove(sprintf(c("%s.all.score", "%s.snp", "%s.log"), output_prefix))), add = TRUE)
  
  data.table::fwrite(object$model, model_file, sep = " ", quote = FALSE, col.names = TRUE)
  
  args$output.file.prefix <- output_prefix
  
  args$base.file <- model_file
  base_options <- list()
  if (!missing(summary.statistics.index)) {
    checkmate::assert_flag(summary.statistics.index, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$index <- summary.statistics.index
  }
  if (!missing(summary.statistics.chr.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.chr.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.chr.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$chr.col <- summary.statistics.chr.col
  }
  if (!missing(summary.statistics.pos.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.pos.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.pos.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$pos.col <- summary.statistics.pos.col
  }
  if (!missing(summary.statistics.snp.id.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.snp.id.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.snp.id.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$snp.id.col <- summary.statistics.snp.id.col
  }
  if (!missing(summary.statistics.effect.allele.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.effect.allele.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.effect.allele.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$effect.allele.col <- summary.statistics.effect.allele.col
  }
  if (!missing(summary.statistics.non.effect.allele.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.non.effect.allele.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.non.effect.allele.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$non.effect.allele.col <- summary.statistics.non.effect.allele.col
  }
  if (!missing(summary.statistics.stat.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.stat.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.stat.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$stat.col <- summary.statistics.stat.col
  }
  if (!missing(summary.statistics.pvalue.col)) {
    if (summary.statistics.index) {
      checkmate::assert_int(summary.statistics.pvalue.col, na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
    } else {
      checkmate::assert_string(summary.statistics.pvalue.col, add = assertions, na.ok = FALSE, null.ok = FALSE)
    }
    base_options$pvalue.col <- summary.statistics.pvalue.col
  }
  if (!missing(summary.statistics.beta)) {
    checkmate::assert_flag(summary.statistics.beta, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$beta <- summary.statistics.beta
  }
  if (!missing(summary.statistics.or)) {
    checkmate::assert_flag(summary.statistics.or, na.ok = FALSE, null.ok = FALSE, add = assertions)
    base_options$or <- summary.statistics.or
  }
  args$base.options <- base_options
  
  if (!missing(target.type)) {
    checkmate::assert_choice(target.type, choices = c("bed", "bgen"), null.ok = FALSE, add = assertions)
  } else {
    target.type <- "bed"
  }
  
  if (target.type == "bed") {
    checkmate::assert_file(sprintf("%s.bed", new.target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", new.target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", new.target.file.prefix), access = "r", add = assertions)
  } else {
    checkmate::assert_file(sprintf("%s.bgen", new.target.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bgen.bgi", new.target.file.prefix), access = "r", add = assertions)
  }
  args$target.file.prefix <- new.target.file.prefix
  
  target_options <- list()
  target_options$type <- target.type
  if (!missing(target.nonfounders)) {
    checkmate::assert_flag(target.nonfounders, na.ok = FALSE, null.ok = FALSE, add = assertions)
    target_options$nonfounders <- target.nonfounders
  }
  if (!missing(target.keep)) {
    checkmate::assert_file(target.keep, access = "r", add = assertions)
    target_options$keep <- target.keep
  }
  args$target.options <- target_options
  
  args$clumping <- FALSE
  
  dosage_options <- list()
  if (!missing(target.hardcall)) {
    checkmate::assert_flag(target.hardcall, na.ok = FALSE, null.ok = FALSE, add = assertions)
    dosage_options$hardcall <- target.hardcall
  }
  if (!missing(target.dosage.threshold)) {
    checkmate::assert_number(target.dosage.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    dosage_options$dosage.threshold <- target.dosage.threshold
  }
  if (!missing(target.hardcall.threshold)) {
    checkmate::assert_number(target.hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    dosage_options$hardcall.threshold <- target.hardcall.threshold
  }
  if (!missing(dosage.allow.intermediate)) {
    checkmate::assert_flag(dosage.options$allow.intermediate, na.ok = FALSE, null.ok = FALSE, add = assertions)
    dosage_options$allow.intermediate <- dosage.allow.intermediate
  }
  args$dosage.options <- dosage_options
  
  pval_options <- list(
    levels = 1,
    full = FALSE
  )
  if (!missing(model)) {
    checkmate::assert_choice(model, choices = c("add", "dom", "rec", "het"), null.ok = FALSE, add = assertions)
    pval_options$model <- model
  }
  args$pval.options <- pval_options
  
  if (!missing(missing.handling)) {
    checkmate::assert_choice(missing.handling, choices = c("IMPUTE", "SET_ZERO", "CENTER"), null.ok = FALSE, add = assertions)
    args$missing.handling <- missing.handling
  }
  
  if (!missing(score)) {
    checkmate::assert_choice(score, choices = c("avg", "std", "sum"), null.ok = FALSE, add = assertions)
    args$score <- score
  }
  
  args$regress <- FALSE
  
  if (!missing(id.delim)) {
    checkmate::assert_string(id.delim, na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
    args$id.delim <- id.delim
  }
  if (!missing(ignore.fid)) {
    checkmate::assert_flag(ignore.fid, na.ok = FALSE, null.ok = FALSE, add = assertions)
    args$ignore.fid <- ignore.fid
  }
  if (!missing(keep.ambig)) {
    checkmate::assert_flag(keep.ambig, na.ok = FALSE, null.ok = FALSE, add = assertions)
    args$keep.ambig <- keep.ambig
  }
  
  copy_to_scratch_dir <- FALSE
  if (!missing(scratch.dir)) {
    checkmate::assert_directory(scratch.dir, access = "rw", add = assertions)
    copy_to_scratch_dir <- TRUE
  }
  
  if (!missing(scratch.dir.files)) {
    checkmate::assert_character(scratch.dir.files, add = assertions)
    checkmate::assert_subset(scratch.dir.files, choices = c("base.file", "target.file", "ld.file"), add = assertions)
    scratch.dir.files <- scratch.dir.files[scratch.dir.files %in% c("target.file")]
  } else {
    scratch.dir.files <- character()
  }
  args$print.snps <- TRUE
  args$seed <- sample(.Machine$integer.max, 1)
  if (!missing(memory)) {
    checkmate::assert_int(memory, lower = 1000, add = assertions)
    args$memory <- memory
  }
  if (!missing(num.threads)) {
    checkmate::assert_int(num.threads, lower = 1, add = assertions)
    args$num.threads <- num.threads
  }
  checkmate::assert_string(prsice.exec, na.ok = FALSE, null.ok = FALSE, add = assertions)
  prsice.exec <- path.expand(prsice.exec)
  imbs::assert_command(prsice.exec, add = assertions)
  args$exec <- prsice.exec
  
  checkmate::reportAssertions(assertions)
  
  # Copy files to scratch.dir ----
  if (copy_to_scratch_dir) {
    file_list <- character()
    if ("target.file" %in% scratch.dir.files) {
      file_list <- c(file_list, sprintf("%s*", new.target.file.prefix))
      args$target.file.prefix <- file.path(scratch.dir, basename(new.target.file.prefix))
    }
    if (length(file_list) > 1) {
      file_list <- sprintf("{%s}", paste(file_list, collapse = ","))
    } else {
      file_list <- sprintf("%s", paste(file_list, collapse = " "))
    }
    if (imbs::test_command("rsync")) {
      imbs::system_call(bin = "rsync", args = c("-avPzh", file_list, scratch.dir))
    } else {
      imbs::system_call(bin = "cp", args = c(file_list, scratch.dir))
    }
  }
  
  do.call(what = prsice, args = args)
  
  data.table::fread(sprintf("%s.all.score", output_prefix), header = TRUE, col.names = c("FID", "IID", "SCORE"))
  
}

makeRLearner.classif.prsice <- function() {
  
  makeRLearnerClassif(
    cl = "classif.prsice",
    package = "data.table",
    par.set = makeParamSet(
      makeUntypedLearnerParam(id = "summary.statistics.chr.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.pos.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.snp.id.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.effect.allele.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.non.effect.allele.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.stat.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.pvalue.col", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.info.col", when = "train", tunable = FALSE),
      makeNumericLearnerParam(id = "summary.statistics.info.threshold", when = "train", tunable = TRUE, lower = 0, upper = 1, allow.inf = FALSE),
      makeUntypedLearnerParam(id = "summary.statistics.maf.cols", when = "train", tunable = FALSE),
      makeNumericLearnerParam(id = "summary.statistics.maf.thresholds", when = "train", tunable = TRUE, lower = 0, upper = 1, allow.inf = FALSE),
      makeLogicalLearnerParam(id = "summary.statistics.beta", when = "train", tunable = FALSE, default = TRUE),
      makeLogicalLearnerParam(id = "summary.statistics.or", when = "train", tunable = FALSE, default = FALSE),
      makeLogicalLearnerParam(id = "summary.statistics.index", when = "train", tunable = FALSE, default = FALSE),
      makeDiscreteLearnerParam(id = "target.type", values = c("bed", "bgen"), when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "target.keep", when = "both", tunable = FALSE, requires = quote(!is.na(ld.file.prefix))),
      makeNumericLearnerParam(id = "target.geno", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE),
      makeNumericLearnerParam(id = "target.maf", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE),
      makeNumericLearnerParam(id = "target.info", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(target.type == "bgen")),
      makeLogicalLearnerParam(id = "target.nonfounders", when = "train", tunable = TRUE, default = FALSE),
      makeLogicalLearnerParam(id = "target.hardcall", when = "both", requires = quote(target.type == "bgen"), default = FALSE),
      makeNumericLearnerParam(id = "target.dosage.threshold", when = "both", lower = 0, upper = 1, allow.inf = FALSE, requires = quote(target.type == "bgen" & isTRUE(target.hardcall))),
      makeNumericLearnerParam(id = "target.hardcall.threshold", when = "both", lower = 0, upper = 0.5, allow.inf = FALSE, requires = quote(target.type == "bgen" & isTRUE(target.hardcall))),
      makeLogicalLearnerParam(id = "ld.external", default = TRUE, when = "train", tunable = TRUE, requires = quote(clumping == TRUE)),
      makeUntypedLearnerParam(id = "ld.file.prefix", when = "train", tunable = FALSE, requires = quote(ld.external == TRUE)),
      makeDiscreteLearnerParam(id = "ld.type", values = c("bed", "bgen"), when = "train", tunable = FALSE, requires = quote(!is.na(ld.file.prefix))),
      makeNumericLearnerParam(id = "ld.geno", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(!is.na(ld.file.prefix))),
      makeNumericLearnerParam(id = "ld.maf", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(!is.na(ld.file.prefix))),
      makeNumericLearnerParam(id = "ld.info", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(ld.type == "bgen")),
      makeUntypedLearnerParam(id = "ld.keep", when = "train", tunable = FALSE, requires = quote(!is.na(ld.file.prefix))),
      makeNumericLearnerParam(id = "ld.dosage.threshold", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(ld.type == "bgen")),
      makeNumericLearnerParam(id = "ld.hardcall.threshold", lower = 0, upper = 0.5, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(ld.type == "bgen")),
      makeLogicalLearnerParam(id = "clumping", when = "train", tunable = TRUE, default = TRUE),
      makeIntegerLearnerParam(id = "clumping.kb", lower = 0, when = "train", tunable = TRUE, requires = quote(clumping == TRUE), default = 250),
      makeNumericLearnerParam(id = "clumping.r2", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(clumping == TRUE), default = 0.1),
      makeNumericLearnerParam(id = "clumping.p", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(clumping == TRUE), default = 1),
      makeNumericLearnerParam(id = "clumping.proxy", lower = 0, upper = 1, allow.inf = FALSE, when = "train", tunable = TRUE, requires = quote(clumping == TRUE)),
      makeLogicalLearnerParam(id = "clumping.pearson", when = "train", tunable = TRUE, default = FALSE),
      makeLogicalLearnerParam(id = "dosage.allow.intermediate", when = "both", requires = quote(ld.type == "bgen" || target.type == "bgen")),
      makeNumericLearnerParam(id = "pval.level", lower = 0, upper = 1, allow.inf = FALSE, default = 1, when = "train", tunable = TRUE),
      makeDiscreteLearnerParam(id = "model", values = c("add", "dom", "rec", "het"), default = "add", when = "train", tunable = TRUE),
      makeDiscreteLearnerParam(id = "missing.handling", values = c("IMPUTE", "SET_ZERO", "CENTER"), when = "train", tunable = TRUE, default = "IMPUTE"),
      makeDiscreteLearnerParam(id = "score", values = c("avg", "std", "sum"), when = "train", tunable = FALSE, default = "sum"),
      makeUntypedLearnerParam(id = "exclude", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "exclude.range", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "extract", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "id.delim", when = "train", tunable = FALSE),
      makeLogicalLearnerParam(id = "ignore.fid", when = "train", default = FALSE, tunable = FALSE),
      makeLogicalLearnerParam(id = "keep.ambig", when = "train", tunable = TRUE, default = FALSE),
      makeUntypedLearnerParam(id = "scratch.dir", when = "both", tunable = FALSE),
      makeUntypedLearnerParam(id = "scratch.dir.files", when = "both", tunable = FALSE),
      makeIntegerLearnerParam(id = "memory", lower = 1000, when = "both", tunable = FALSE),
      makeIntegerLearnerParam(id = "num.threads", lower = 1, when = "both", tunable = FALSE),
      makeUntypedLearnerParam(id = "prsice.exec", default = "prsice", when = "both", tunable = FALSE),
      makeUntypedLearnerParam(id = "plink2.exec", default = "plink2", when = "train", tunable = FALSE)
    ),
    properties = c("twoclass", "numerics", "factors", "prob"),
    name = "Polygenic Risk Score with PRSice",
    short.name = "prsice",
    note = ""
  )
  
}

trainLearner.classif.prsice <- function(.learner, .task, .subset, .weights = NULL, ...) {
  
  samples <- getTaskData(.task, recode.target = "01")
  samples_file <- tempfile()
  on.exit(BBmisc::suppressAll(file.remove(samples_file)), add = TRUE)
  
  data.table::fwrite(samples, file = samples_file, sep = " ", quote = FALSE, col.names = FALSE)
  
  learner_par_vals <- getLearnerParVals(.learner)
  function_par_vals <- list(...)
  ps <- getLearnerParamSet(.learner)$pars
  ns <- Filter(function(x) ps[[x]]$when %in% c("train", "both"), names(ps))
  par_defaults <- getDefaults(ps[ns])
  train_par_vals <- c(learner_par_vals[setdiff(names(learner_par_vals), names(function_par_vals))], function_par_vals)
  par_vals <- c(par_defaults[setdiff(names(par_defaults), names(train_par_vals))], train_par_vals)
  
  par_vals[["target.file.prefix"]] <- path.expand(sprintf("%s", unique(samples[["TARGET_FILE_PREFIX"]])))
  par_vals[["target.keep"]] <- samples_file
  
  ssf <- samples[["SUMMARY_STATISTICS_FILE"]] %??% FALSE
  if (!checkmate::test_logical(ssf)) {
    # If summary statistis file is given use it, otherwise a GWA using PLINK2 is done
    par_vals[["base.file"]] <- path.expand(sprintf("%s", unique(samples[["SUMMARY_STATISTICS_FILE"]])))
  }
  
  do.call(what = create_prsice_model, args = par_vals)
  
}

predictLearner.classif.prsice <- function(.learner, .model, .newdata, predict.method = "prob", ...) {
  
  model_par_vals <- getLearnerParVals(.model$learner)
  function_par_vals <- list(...)
  ps <- getLearnerParamSet(.learner)$pars
  ns <- Filter(function(x) ps[[x]]$when %in% c("predict", "both"), names(ps))
  learner_par_defaults <- getDefaults(ps[ns])
  predict_par_vals <- c(model_par_vals[setdiff(names(model_par_vals), names(function_par_vals))], function_par_vals)
  args <- c(learner_par_defaults[setdiff(names(learner_par_defaults), names(predict_par_vals))], predict_par_vals)
  
  if (exists("scratch.dir", where = args)) {
    samples_file <- tempfile(tmpdir = args[["scratch.dir"]])
  } else {
    samples_file <- tempfile()
  }
  on.exit(BBmisc::suppressAll(file.remove(samples_file)), add = TRUE)
  
  args[["new.target.file.prefix"]] <- path.expand(sprintf("%s", unique(.newdata[["TARGET_FILE_PREFIX"]])))
  args[["target.keep"]] <- samples_file
  args[["object"]] <- mlr::getLearnerModel(model = .model, more.unwrap = TRUE)
  
  data.table::fwrite(.newdata, file = samples_file, sep = " ", quote = FALSE, col.names = FALSE)
  
  scores <- do.call(predict.prsice, args = args)
  scores[, IID := factor(IID)]
  scores[, FID := factor(FID)]
  
  prob <- scores[.newdata, on = c("FID" = "FID", "IID" = "IID")][, .(good = 1-(SCORE - args[["object"]][["min.score"]])/(args[["object"]][["max.score"]] - args[["object"]][["min.score"]]), bad = (SCORE - args[["object"]][["min.score"]])/(args[["object"]][["max.score"]] - args[["object"]][["min.score"]]))]
  
  if (.learner$predict.type == "response") {
    return(factor(c(mlr::getTaskDesc(.model)[["negative"]], mlr::getTaskDesc(.model)[["positive"]])[1 + round(prob$bad)]))
  } else {
    data.table::setnames(
      x = prob, 
      old = c("good", "bad"), 
      new = c(mlr::getTaskDesc(.model)[["negative"]], mlr::getTaskDesc(.model)[["positive"]])
    )
    return(as.matrix(prob))
  }
  
}

registerS3method(
  genname = "predict", 
  class = "prsice", 
  method = predict.prsice
)
registerS3method(
  genname = "makeRLearner", 
  class = "classif.prsice", 
  method = makeRLearner.classif.prsice
)
registerS3method(
  genname = "trainLearner", 
  class = "classif.prsice", 
  method = trainLearner.classif.prsice
)
registerS3method(
  genname = "predictLearner",
  class = "classif.prsice", 
  method = predictLearner.classif.prsice
)