
#' Convert BGEN to PLINK2 binary file set
#'
#' @param bgen.file                [\code{string}]\cr
#'                                 File path to BGEN file.
#' @param sample.file              [\code{string}]\cr
#'                                 File path to accompanying sample file.
#' @param output.prefix            [\code{string}]\cr
#'                                 File path prefix of PLINK2 binary file set.
#' @param import.dosage.certainty  [\code{number}]\cr
#'                                 Lower limit for highest genotype probability to be converted into a non-missing dosage.
#' @param ...                      Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec              [\code{string}]\cr
#'                                 Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads              [\code{int}]\cr
#'                                 Number of CPUs usable by PLINK.
#'                                 Default is determined by SLURM environment variables and at least 1.
#' @param memory                   [\code{int}]\cr
#'                                 Memory for PLINK in Mb.
#'                                 Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#'
plink2_bgen_conversion <- function(bgen.file, sample.file, output.prefix, import.dosage.certainty, num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(bgen.file, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(import.dosage.certainty)) {
    import.dosage.certainty <- ""
  } else {
    checkmate::assert_number(import.dosage.certainty, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assersionts)
    import.dosage.certainty <- sprintf("--import-dosage-certainty %f", import.dosage.certainty)
  }
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      "--bgen", bgen.file,
      "--sample", sample.file,
      "--make-pgen", "vzs",
      import.dosage.certainty,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads
    )
  )
  
}

#' Convert BGEN to PLINK1 binary file set
#'
#' @param bgen.file                [\code{string}]\cr
#'                                 File path to BGEN file.
#' @param sample.file              [\code{string}]\cr
#'                                 File path to accompanying sample file.
#' @param output.prefix            [\code{string}]\cr
#'                                 File path prefix of PLINK2 binary file set.
#' @param import.dosage.certainty  [\code{number}]\cr
#'                                 Lower limit for highest genotype probability to be converted into a non-missing dosage.
#' @param hardcall.threshold       [\code{string}]\cr
#'                                 How close must a dosage be to the next integer to be converted to a valid hardcall?
#' @param ...                      Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec              [\code{string}]\cr
#'                                 Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads              [\code{int}]\cr
#'                                 Number of CPUs usable by PLINK.
#'                                 Default is determined by SLURM environment variables and at least 1.
#' @param memory                   [\code{int}]\cr
#'                                 Memory for PLINK in Mb.
#'                                 Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#'
plink2_bgen2bed_conversion <- function(bgen.file, sample.file, output.prefix, import.dosage.certainty, hardcall.threshold, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(bgen.file, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(import.dosage.certainty)) {
    import.dosage.certainty <- ""
  } else {
    checkmate::assert_number(import.dosage.certainty, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assersionts)
    import.dosage.certainty <- sprintf("--import-dosage-certainty %f", import.dosage.certainty)
  }
  if (missing(hardcall.threshold)) {
    hardcall.threshold <- ""
  } else {
    checkmate::assert_number(hardcall.threshold, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE, add = assertions)
    hardcall.threshold <- sprintf("--hard-call-threshold %f", hardcall.threshold)
  }
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      "--bgen", bgen.file,
      "--sample", sample.file,
      "--make-bed",
      import.dosage.certainty,
      hardcall.threshold,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' Convert GEN to PLINK1 binary file set
#'
#' @param gen.file                 [\code{string}]\cr
#'                                 File path to BGEN file.
#' @param sample.file              [\code{string}]\cr
#'                                 File path to accompanying sample file.
#' @param output.prefix            [\code{string}]\cr
#'                                 File path prefix of PLINK2 binary file set.
#' @param import.dosage.certainty  [\code{number}]\cr
#'                                 Lower limit for highest genotype probability to be converted into a non-missing dosage.
#' @param hardcall.threshold       [\code{string}]\cr
#'                                 How close must a dosage be to the next integer to be converted to a valid hardcall?
#' @param ...                      Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec              [\code{string}]\cr
#'                                 Path of PLINK2 executable.
#' @param num.threads              [\code{int}]\cr
#'                                 Number of CPUs usable by PLINK.
#'                                 Default is determined by SLURM environment variables and at least 1.
#' @param memory                   [\code{int}]\cr
#'                                 Memory for PLINK in Mb.
#'                                 Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#'
plink2_gen2bed_conversion <- function(gen.file, sample.file, output.prefix, import.dosage.certainty, hardcall.threshold, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(gen.file, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(import.dosage.certainty)) {
    import.dosage.certainty <- ""
  } else {
    checkmate::assert_number(import.dosage.certainty, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assersionts)
    import.dosage.certainty <- sprintf("--import-dosage-certainty %f", import.dosage.certainty)
  }
  if (missing(hardcall.threshold)) {
    hardcall.threshold <- ""
  } else {
    checkmate::assert_number(hardcall.threshold, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE, add = assertions)
    hardcall.threshold <- sprintf("--hard-call-threshold %f", hardcall.threshold)
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      "--gen", gen.file,
      "--sample", sample.file,
      "--make-bed",
      import.dosage.certainty,
      hardcall.threshold,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' Convert VCF to PLINK2 binary file set
#'
#' @param vcf.file       [\code{string}]\cr
#'                       File path to VCF file.
#' @param output.prefix  [\code{string}]\cr
#'                       File path prefix of PLINK2 binary file set.
#' @param ...            Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec    [\code{string}]\cr
#'                       Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads    [\code{int}]\cr
#'                       Number of CPUs usable by PLINK.
#'                       Default is determined by SLURM environment variables and at least 1.
#' @param memory         [\code{int}]\cr
#'                       Memory for PLINK in Mb.
#'                       Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#'
plink2_vcf_conversion <- function(vcf.file, output.prefix, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      "--vcf", vcf.file, "dosage=DS",
      "--make-pgen", "vzs",
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads
    )
  )
  
}

#' Convert VCF to PLINK1 binary file set
#'
#' @param vcf.file                 [\code{string}]\cr
#'                                 File path to VCF file.
#' @param output.prefix            [\code{string}]\cr
#'                                 File path prefix of PLINK2 binary file set.
#' @param hardcall.threshold       [\code{string}]\cr
#'                                 How close must a dosage be to the next integer to be converted to a valid hardcall?
#' @param ...                      Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec              [\code{string}]\cr
#'                                 Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads              [\code{int}]\cr
#'                                 Number of CPUs usable by PLINK.
#'                                 Default is determined by SLURM environment variables and at least 1.
#' @param memory                   [\code{int}]\cr
#'                                 Memory for PLINK in Mb.
#'                                 Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#'
plink2_vcf2bed_conversion <- function(vcf.file, output.prefix, hardcall.threshold, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(hardcall.threshold)) {
    hardcall.threshold <- ""
  } else {
    checkmate::assert_number(hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    hardcall.threshold <- sprintf("--hard-call-threshold %f", hardcall.threshold)
  }
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      "--vcf", vcf.file, "dosage=DS",
      "--make-bed",
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' Convert PLINK2 binary file set to PLINK1 binary file set
#'
#' @param pfile.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK2 binary file set.
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of PLINK1 binary file set.
#' @param hardcall.threshold   [\code{string}]\cr
#'                             How close must a dosage be to the next integer to be converted to a valid hardcall?
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_make_bed <- function(pfile.prefix, vzs, output.prefix, hardcall.threshold, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(pfile.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_file(sprintf("%s.pgen", pfile.prefix), access = "r", add = assertions)
  checkmate::assert_file(sprintf("%s.psam", pfile.prefix), access = "r", add = assertions)
  if (missing(vzs)) {
    checkmate::assert_file(sprintf("%s.pvar", pfile.prefix), access = "r", add = assertions)
    pfile.prefix <- sprintf("--pfile %s", pfile.prefix) 
  } else {
    checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
    checkmate::assert_file(sprintf("%s.pvar.zst", pfile.prefix), access = "r", add = assertions)
    pfile.prefix <- sprintf("--pfile %s vzs", pfile.prefix)
  }
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(hardcall.threshold)) {
    hardcall.threshold <- ""
  } else {
    checkmate::assert_number(hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    hardcall.threshold <- sprintf("--hard-call-threshold %f", hardcall.threshold)
  }
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      pfile.prefix,
      "--make-bed",
      hardcall.threshold,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' Find closely related individuals using a two-step approach of PLINK2
#'
#' @param input.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK2 binary file set.
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of PLINK1 binary file set.
#' @param screen.maf           [\code{number}]\cr
#'                             Remove variants with MAF lower than this threshold during screening step.
#' @param screen.thin          [\code{number}]\cr
#'                             Thin dataset down to this percentage.
#' @param screen.degree        [\code{number}]\cr
#'                             Only report pairs of individuals having at least the given degree of relation during screening step.
#' @param degree               [\code{number}]\cr
#'                             Only report pairs of individuals having at least the given degree of relation during accurate kinship-coefficient computations. 
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details This function starts with a screening step which considers all sample pairs but only a small number of variants scattered across the genome, and follows up with accurate kinship-coefficient computations for just the sample pairs identified as possible close relations during the screening step.
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_fast_related <- function(input.prefix, output.prefix, screen.maf, screen.thin, screen.degree, degree, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(input.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (checkmate::test_file_exists(sprintf("%s.pgen", input.prefix))) {
    checkmate::assert_file(sprintf("%s.pgen", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.psam", input.prefix), access = "r", add = assertions)
    if (missing(vzs)) {
      input.prefix <- sprintf("--pfile %s", input.prefix)
      checkmate::assert_file(sprintf("%s.pvar", input.prefix), access = "r", add = assertions) 
    } else {
      checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
      checkmate::assert_file(sprintf("%s.pvar.zst", input.prefix), access = "r", add = assertions)
      input.prefix <- sprintf("--pfile %s vzs", input.prefix)
    } 
  } else {
    checkmate::assert_file(sprintf("%s.bed", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", input.prefix), access = "r", add = assertiosn)
    input.prefix <- sprintf("--bfile %s", input.prefix)
  }
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  if (missing(screen.maf)) {
    screen.maf <- ""
  } else {
    checkmate::assert_number(screen.maf, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    screen.maf <- sprintf("--maf %f", screen.maf)
  }
  if (missing(screen.thin)) {
    screen.thin <- ""
  } else {
    checkmate::assert_number(screen.thin, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    screen.thin <- sprintf("--thin %f", screen.thin)
  }
  if (missing(screen.degree)) {
    screen.degree <- ""
  } else {
    checkmate::assert_number(screen.degree, na.ok = FALSE, lower = 0, upper = 1, finit = TRUE, null.ok = FALSE, add = assertions)
    screen.degree <- sprintf("--king-table-filter %f", screen.degree)
  }
  if (missing(degree)) {
    degree <- ""
  } else {
    checkmate::assert_number(degree, na.ok = FALSE, lower = 0, upper = 1, finit = TRUE, null.ok = FALSE, add = assertions)
    degree <- sprintf("--king-table-filter %f", degree)
  }
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  screen_log <- imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      "--make-king-table",
      screen.degree,
      "--out", sprintf("%s.screen", output.prefix),
      screen.maf,
      screen.thin,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
  fine_log <- imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      "--make-king-table",
      degree,
      "--king-table-subset", sprintf("%s.screen.kin0", output.prefix),
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
  return(
    list(
      screen_log = screen_log,
      fine_log = fine_log
    )
  )
  
}

#' Apply PLINK2 linear scoring
#'
#' @param input.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK1 or PLINK2 binary file set. 
#' @param score.file           [\code{string}]\cr
#'                             File path to scores per variant.
#' @param var.id.col.idx       [\code{int}]\cr
#'                             Column index of variant IDs in \code{score.file}.
#' @param allele.col.idx       [\code{int}]\cr
#'                             Column index of effect allele code in \code{score.file}.
#' @param coef.col.idx         [\code{int}]\cr
#'                             Column index of score coefficient in \code{score.file}.
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of output files. 
#' @param header               [\code{flag}]\cr
#'                             Has the \code{score.file} a header line?
#' @param header.read          [\code{flag}]\cr
#'                             Should the header of \code{score.file} be included in output?
#' @param transformation       [\code{string}]\cr
#'                             Any of \code{center} or \center{standardize}.
#'                             If \code{center} shift alle genotypes to mean zero, if \code{standardize} linearly transform the genotypes to mean 0, variance 1.
#' @param model                [\code{string}]\cr
#'                             Any of \code{additive}, \code{dominant}, \code{recessive}.
#'                             If \code{additive} allele dosages are used, if \code{dominant} dosages grater than 1 are treated as 1, if \code{recessive} dosages are transformed by min(dosage-1,0).
#' @param mean.imputation      [\code{flag}]\cr
#'                             Missing genotypes contribute an amount proportional to the imputed allele frequency.
#'                             You can use "--read-freq <file>" in \code{...} to load pre-computed minor allele frequencies.
#' @param se                   [\code{flag}]\cr
#'                             Treat score coefficients as independent standard errors; in this case, standard errors for the score average/sum are reported.
#' @param zs                   [\code{flag}]\cr
#'                             Compress output files using Zstandard compression.
#' @param ignore.dup.ids       [\code{flag}]\cr
#'                             If \code{TRUE} variants with duplicated IDs are skipped (with a warning). If \code{FALSE} and duplicated IDs are found, PLINK2 will error out.
#' @param list.variants        [\code{flag}]\cr
#'                             Write out variant IDs used for scoring to <output.prefix>.sscore.vars{|.zst}.
#' @param output.cols          [\code{character}]\cr
#'                             Control the main report (<output.prefix>.sscore) column set:
#'                             \describe{
#'                               \item[maybefid]{FID, if that column was in the input.}
#'                               \item[fid]{Force FID column to be written even when absent in the input.}
#'                               \item[]{(IID is always present, and positioned here.)}
#'                               \item[maybesid]{SID, if that column was in the input.}
#'                               \item[sid]{Force SID column to be written even when absent in the input.}
#'                               \item[pheno1]{First active phenotype.}
#'                               \item[phenos]{All active phenotypes, if any.}
#'                               \item[nmissallele]{Number of nonmissing alleles.}
#'                               \item[denom]{Denominator of score average (equal to nmissallele value when \code{mean.imputation} set to \code{FALSE}).}
#'                               \item[dosagesum]{Sum of named allele dosages.}
#'                               \item[scoreavgs]{Score averages.}
#'                               \item[scoresums]{Score sums.}
#'                             }
#'                             The default is maybefid, maybesid, phenos, nmissallele, dosagesum, scoreavgs.
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#' 
#' @return Captured system outputs as \code{list} of \code{character} vectors.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_scoring <- function(input.prefix, score.file, var.id.col.idx, allele.col.idx, coef.col.idx, output.prefix, header, header.read, transformation, model, mean.imputation, se, zs, ignore.dup.ids, list.variants, output.cols, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(input.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (checkmate::test_file_exists(sprintf("%s.pgen", input.prefix))) {
    checkmate::assert_file(sprintf("%s.pgen", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.psam", input.prefix), access = "r", add = assertions)
    if (missing(vzs)) {
      input.prefix <- sprintf("--pfile %s", input.prefix)
      checkmate::assert_file(sprintf("%s.pvar", input.prefix), access = "r", add = assertions) 
    } else {
      checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
      checkmate::assert_file(sprintf("%s.pvar.zst", input.prefix), access = "r", add = assertions)
      input.prefix <- sprintf("--pfile %s vzs", input.prefix)
    } 
  } else {
    checkmate::assert_file(sprintf("%s.bed", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", input.prefix), access = "r", add = assertiosn)
    input.prefix <- sprintf("--bfile %s", input.prefix)
  }
  checkmate::assert_file(score.prefix, access = "r", add = assertions)
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  checkmate::assert_int(var.id.col.idx, na.ok = FALSE, null.ok = FALSE, add = assertions, lower = 1)
  checkmate::assert_int(coef.col.idx, na.ok = FALSE, null.ok = FALSE, add = assertions, lower = 1)
  checkmate::assert_int(allele.col.idx, na.ok = FALSE, null.ok = FALSE, add = assertions, lower = 1)
  if (!missing(header)) {
    checkmate::assert_flag(header, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (header) {
      header <- "header"
      if (!missing(header.read)) {
        checkmate::assert_flag(header.read, na.ok = FALSE, null.ok = FALSE, add = assertions)
        if (header.read) {
          header.read <- "header-read"
        } else {
          header.read <- ""
        }
      } else {
        header.read <- ""
      }
    } else {
      header <- ""
    }
  } else {
    header <- ""
  }
  if (!missing(transformation)) {
    checkmate::assert_choice(transformation, choices = c("center", "standardize"), null.ok = FALSE, add = assertions)
  } else {
    transformation <- ""
  }
  if (!missing(model)) {
    checkmate::assert_choice(model, choices = c("additive", "dominant", "recessive"), null.ok = FALSE, add = assertions)
    if (model == "additive") {
      model <- ""
    }
  } else {
    model <- ""
  }
  if (!missing(mean.imputation)) {
    checkmate::assert_flag(mean.imputation, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (!mean.imputation) {
      mean.imputation <- "no-mean-imputation"
    } else {
      mean.imputation <- ""
    }
  } else {
    mean.imputation <- ""
  }
  if (!missing(se)) {
    checkmate::assert_flag(se, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (se) {
      se <- "se"
    } else {
      se <- ""
    }
  } else {
    zs <- ""
  }
  if (!missing(zs)) {
    checkmate::assert_flag(zs, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (zs) {
      zs <- "zs"
    } else {
      zs <- ""
    }
  } else {
    zs <- ""
  }
  if (!missing(ignore.dup.ids)) {
    checkmate::assert_flag(ignore.dup.ids, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (ignore.dup.ids) {
      ignore.dup.ids <- "ignore-dup-ids"
    } else {
      ignore.dup.ids <- ""
    }
  } else {
    ignore.dup.ids <- ""
  }
  if (!missing(list.variants)) {
    checkmate::assert_flag(list.variants, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (list.variants) {
      if (zs) {
        list.variants <- "list-variants-zs"
      } else {
        list.variants <- "list-variants"
      }
    } else {
      list.variants <- ""
    } 
  } else {
    list.variants <- ""
  }
  if (!missing(output.cols)) {
    checkmate::assert_character(output.cols, any.missing = FALSE, all.missing = FALSE, unique = TRUE, null.ok = FALSE, min.len = 1, add = assertions)
    checkmate::assert_subset(output.cols, choices = c("maybefid", "fid", "maybesid", "sid", "pheno1", "phenos", "nmissallele", "denom", "dosagesum", "scoreavgs", "scoresums"), add = assertions)
    output.cols <- sprintf("cols=%s", paste(output.cols, collapse = ","))
  } else {
    output.cols <- ""
  }
  
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      "--score", score.file, var.id.col.idx, allele.col.idx, coef.col.idx, header, header.read, transformation, model, mean.imputation, se, zs, ignore.dup.ids, list.variants, output.cols,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' Rund a GWA (linear/logistic regression) using PLINK2
#'
#' @param input.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK1 or PLINK2 binary file set.  
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of output files.  
#' @param zs                   [\code{flag}]\cr
#'                             Compress output files using Zstandard compression. 
#' @param omit.ref             [\code{flag}]\cr
#'                             There is usually an additive effect line for every nonmajor allele, and no such line for the major allele.  
#'                             To omit REF alleles instead of major alleles, add the \code{omit.ref} modifier.  
#'                             (When performing interaction testing, this tends to cause the multicollinearity check to fail for low-ref-frequency variants.)
#' @param sex                  [\code{flag}]\cr 
#'                             By default, sex (male = 1, female = 2; note that this is a change from PLINK 1.x) is automatically added as a predictor for X chromosome variants, and no others.  
#'                             The \code{sex} modifier causes it to be added everywhere (except chrY).
#' @param no.x.sex             [\code{flag}]\cr 
#'                             By default, sex (male = 1, female = 2; note that this is a change from PLINK 1.x) is automatically added as a predictor for X chromosome variants, and no others.  
#'                             The \code{no.x.sex} modifier excludes it entirely.
#' @param log10                [\code{flag}]\cr 
#'                             The \code{log10} modifier causes p-values to be reported in \eqn{-log10(p)} form.
#' @param pheno.ids            [\code{flag}]\cr
#'                             The \code{pheno.ids} modifier causes the samples used in each set of regressions to be written to an .id file. 
#'                             (When the samples differ on chrX or chrY, .x.id and/or .y.id files are also written.)
#' @param model                [\code{string}]\cr
#'                             Default is the \code{additive} model.
#'                             \code{dominant} and \code{recessive} specify a model assuming full dominance or recessiveness, respectively, for the ref allele.  
#'                             I.e. the genotype column is recoded as 0..1..1 or 0..0..1, respectively.
#'                             The \code{genotypic} model adds an additive effect/dominance deviation 2df joint test (0-2 and 0..1..0 coding), while 'hethom' uses 0..0..1 and 0..1..0 coding instead.
#' @param interaction          [\code{flag}]\cr
#'                             \code{interaction} adds genotype x covariate interactions to the model.  
#'                             Note that this tends to produce 'NA' results (due to the multicollinearity check) when the reference allele is 'wrong'; "--maj-ref" can be used to enable analysis of those variants.   
#' @param intercept            [\code{flag}]\cr
#'                             ??? 
#' @param hide.covar           [\code{flag}]\cr
#'                             ???  
#' @param firth                [\code{string}]\cr
#'                             For logistic regression, when the phenotype {,quasi-}separates the genotype, an NA result is currently reported by default.  
#'                             To fall back on Firth logistic regression instead when the basic logistic regression fails to converge, use the \code{firth = "firth-fallback"}.
#'                             To eliminate the special case and use Firth logistic regression everywhere, use \code{"firth"}.
#'                             \code{"no-firth"} can be used to prevent Firth regression from being attempted in a way that'll still work after alpha testing completes.
#' @param local.covar.file     [\code{string}]\cr
#'                             File path to variant specific covariates.
#'                             Normally, the \code{local.covar.file} should have \eqn{c * n} real-valued columns, where the first \eqn{c} columns correspond to the first sample in the \code{local.psam.file}, columns \eqn{(c+1)} to \eqn{2c} correspond to the second sample, etc.; and the \eqn{m}th line correspond to the \eqn{m}th nonheader line of the \code{local.pvar.file}. 
#'                             (Variants outside of the \code{local.pvar.file} are excluded from the regression.)  
#'                             The local covariates are assigned the names LOCAL1, LOCAL2, etc.
#' @param local.pvar.file      [\code{string}]\cr
#'                             File path to samples with variant specific covariates. 
#' @param local.psam.file      [\code{string}]\cr
#'                             File path to variants with specific covariates. 
#' @param local.omit.last      [\code{flag}]\cr
#'                             To exclude the last local covariate from the regression (necessary if they are e.g. local ancestry coefficients which sum to 1), use this modifier. 
#' @param local.cats           [\code{flag}]\cr 
#'                             Alternatively, with \code{local.cats=k}, the \code{local.covar.file} is expected to have \eqn{n} columns with integer-valued entries in \eqn{[1, k]}.  
#'                             These category assignments are expanded into \eqn{(k-1)} local covariates in the usual manner.
#' @param output.cols          [\code{character}]\cr
#'                             Control the main report (<output.prefix>.glm{.logistic}) column set:
#'                             \describe{
#'                               \item[chrom]{Chromosome ID.}
#'                               \item[pos]{Base-pair coordinate.}
#'                               \item[]{(ID is always present, and positioned here.)}
#'                               \item[ref]{Reference allele.}
#'                               \item[alt1]{Alternate allele 1.}
#'                               \item[alt]{All alternate alleles, comma-separated.}
#'                               \item[]{(A1 is always present, and positioned here. For multiallelic variants, this column may contain multiple comma-separated alleles when the result doesn't depend on which allele is A1.)}
#'                               \item[ax]{Non-A1 alleles, comma-separated.}
#'                               \item[a1count]{A1 allele count (can be decimal with dosage data).}
#'                               \item[totallele]{Allele observation count (can be higher than "--freq" value, due to inclusion of het haploids and chrX model).}
#'                               \item[a1countcc]{A1 count in cases, then controls (case/control only).}
#'                               \item[totallelecc]{Case and control allele observation counts.}
#'                               \item[gcountcc]{Genotype hardcall counts (neither-A1, het-A1, A1-A1) in cases, then controls (case/control only).}
#'                               \item[a1freq]{A1 allele frequency.}
#'                               \item[a1freqcc]{A1 frequency in cases, then controls (case/control only).}
#'                               \item[machr2]{Unphased MaCH imputation quality (frequently labeled 'INFO').}
#'                               \item[firth]{Reports whether Firth regression was used (firth-fallback only).}
#'                               \item[test]{Test identifier.  (Required unless only one test is run.)}
#'                               \item[nobs]{Number of samples in the regression.}
#'                               \item[beta]{Regression coefficient (for A1 if additive test).}
#'                               \item[orbeta]{Odds ratio for case/control, beta for quantitative traits. (Ignored if 'beta' column set included.)}
#'                               \item[se]{Standard error of beta.}
#'                               \item[ci]{Bounds of symmetric approximate confidence interval (requires --ci).}
#'                               \item[tz]{T-statistic for linear regression, Wald Z-score for logistic/Firth.}
#'                               \item[p]{Asymptotic \eqn{p}-value (or \eqn{-log10(p)}) for T/Z-statistic.}
#'                             }
#'                             The default is chrom,pos,ref,alt,firth,test,nobs,orbeta,se,ci,tz,p.
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#' 
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_gwa <- function(input.prefix, output.prefix, zs, omit.ref, sex, no.x.sex, log10, pheno.ids, model, interaction, intercept, hide.covar, firth, local.covar.file, local.pvar.file, local.psam.file, local.omit.last, local.cats, output.cols, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(input.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (checkmate::test_file_exists(sprintf("%s.pgen", input.prefix))) {
    checkmate::assert_file(sprintf("%s.pgen", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.psam", input.prefix), access = "r", add = assertions)
    if (missing(vzs)) {
      input.prefix <- sprintf("--pfile %s", input.prefix)
      checkmate::assert_file(sprintf("%s.pvar", input.prefix), access = "r", add = assertions) 
    } else {
      checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
      checkmate::assert_file(sprintf("%s.pvar.zst", input.prefix), access = "r", add = assertions)
      input.prefix <- sprintf("--pfile %s vzs", input.prefix)
    } 
  } else {
    checkmate::assert_file(sprintf("%s.bed", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", input.prefix), access = "r", add = assertiosn)
    input.prefix <- sprintf("--bfile %s", input.prefix)
  }
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  
  if (!missing(zs)) {
    checkmate::assert_flag(zs, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (zs) {
      zs <- "zs"
    } else {
      zs <- ""
    }
  } else {
    zs <- ""
  }
  if (!missing(omit.ref)) {
    checkmate::assert_flag(omit.ref, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (omit.ref) {
      omit.ref <- "omit-ref"
    } else {
      omit.ref <- ""
    }
  } else {
    omit.ref <- ""
  }
  if (!missing(sex)) {
    checkmate::assert_flag(sex, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (sex) {
      sex <- "sex"
    } else {
      sex <- ""
    }
  } else {
    sex <- ""
  }
  if (!missing(no.x.sex)) {
    checkmate::assert_flag(no.x.sex, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (no.x.sex) {
      no.x.sex <- "no-x-sex"
    } else {
      no.x.sex <- ""
    }
  } else {
    no.x.sex <- ""
  }
  if (!missing(log10)) {
    checkmate::assert_flag(log10, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (log10) {
      log10 <- "log10"
    } else {
      log10 <- ""
    }
  } else {
    log10 <- ""
  }
  if (!missing(pheno.ids)) {
    checkmate::assert_flag(pheno.ids, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (pheno.ids) {
      pheno.ids <- "pheno-ids"
    } else {
      pheno.ids <- ""
    }
  } else {
    pheno.ids <- ""
  }
  if (!missing(model)) {
    checkmate::assert_choice(model, choices = c("genotypic", "hethom", "dominant", "recessive", "additive"), null.ok = FALSE, add = assertions)
    if (model == "additive") {
      model <- ""
    } else {
      model <- "model"
    }
  } else {
    model <- ""
  }
  if (!missing(interaction)) {
    checkmate::assert_flag(interaction, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (interaction) {
      interaction <- "interaction"
    } else {
      interaction <- ""
    }
  } else {
    interaction <- ""
  }
  if (!missing(intercept)) {
    checkmate::assert_flag(intercept, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (intercept) {
      intercept <- "intercept"
    } else {
      intercept <- ""
    }
  } else {
    intercept <- ""
  }
  if (!missing(hide.covar)) {
    checkmate::assert_flag(hide.covar, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (hide.covar) {
      hide.covar <- "hide-covar"
    } else {
      hide.covar <- ""
    }
  } else {
    hide.covar <- ""
  }
  if (!missing(firth)) {
    checkmate::assert_choice(firth, choices = c("no-firth", "firth-fallback", "firth"), null.ok = FALSE, add = assertions)
  } else {
    firth <- ""
  }
  if (!missing(local.covar.file)) {
    checkmate::assert_file(local.covar.file, access = "r", add = assertions)
  } else {
    local.covar.file <- ""
  }
  if (!missing(local.pvar.file)) {
    checkmate::assert_file(local.pvar.file, access = "r", add = assertions)
  } else {
    local.pvar.file <- ""
  }
  if (!missing(local.psam.file)) {
    checkmate::assert_file(local.psam.file, access = "r", add = assertions)
  } else {
    local.psam.file <- ""
  }
  if (!missing(local.omit.last)) {
    checkmate::assert_flag(local.omit.last, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (local.omit.last) {
      local.omit.last <- "local-omit-last"
    } else {
      local.omit.last <- ""
    }
  } else {
    local.omit.last <- ""
  }
  if (!missing(local.cats)) {
    checkmate::assert_int(local.cats, na.ok = FALSE, lower = 1, null.ok = FALSE, add = assertions)
    local.cats <- sprintf("local-cats=%d", local.cats)
  } else {
    local.cats <- ""
  }
  
  if (!missing(output.cols)) {
    checkmate::assert_character(output.cols, any.missing = FALSE, all.missing = FALSE, unique = TRUE, null.ok = FALSE, min.len = 1, add = assertions)
    checkmate::assert_subset(output.cols, choices = c("chrom", "pos", "ref", "alt1", "alt", "ax", "a1count", "totallele", "a1countcc", "totallelecc", "gcountcc", "a1freq", "a1freqcc", "machr2", "firth", "test", "nobs", "beta", "orbeta", "se", "ci", "tz", "p"), add = assertions)
    output.cols <- sprintf("cols=%s", paste(output.cols, collapse = ","))
  } else {
    output.cols <- ""
  }
  
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      "--glm", zs, omit.ref, sex, no.x.sex, log10, pheno.ids, model, interaction, intercept, hide.covar, firth, local.covar.file, local.pvar.file, local.psam.file, local.omit.last, local.cats, output.cols,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

#' PCA with PLINK2
#'
#' @param input.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK1 or PLINK2 binary file set.  
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of output files.  
#' @param var.wts              [\code{flag}]\cr
#'                             The \code{var.wts} modifier requests an additional <output.prefix>.eigenvec.var file with PCs expressed as variant weights instead of sample weights.
#' @param num.pcs              [\code{int}]\cr
#'                             By default, 10 PCs are extracted; you can adjust this by passing a numeric parameter. 
#'                             (Note that 10 is lower than the PLINK 1.9 default of 20; this is due to the randomized algorithm's memory footprint growing quadratically w.r.t. the PC count.)
#' @param approx               [\code{flag}]\cr
#'                             The \code{approx} modifier causes the standard deterministic computation to be replaced with the randomized algorithm originally implemented for Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ, Price AL (2016) Fast Principal-Component Analysis Reveals Convergent Evolution of ADH1B in Europe and East Asia.  
#'                             This can be a good idea when you have >5k samples.
#' @param meanimpute           [\code{flag}]\cr
#'                             The randomized algorithm always uses mean imputation for missing genotype calls.  F
#'                             or comparison purposes, you can use the \code{meanimpute} modifier to request this behavior for the standard computation.
#' @param vzs                  [\code{flag}]\cr
#'                             The \code{vzs} modifier causes the <output.prefix>.eigenvec.var file to be Zstd-compressed.
#' @param output.cols          [\code{character}]\cr
#'                             Control the main report (<output.prefix>.glm{.logistic}) column set:
#'                             \describe{
#'                               \item[maybefid]{FID, if that column was in the input.}
#'                               \item[fid]{Force FID column to be written even when absent in the input.}
#'                               \item[maybesid]{SID, if that column was in the input.}
#'                               \item[sid]{Force SID column to be written even when absent in the input.}
#'                               \item[chrom]{Chromosome ID.}
#'                               \item[pos]{Base-pair coordinate.}
#'                               \item[]{(ID is always present, and positioned here.)}
#'                               \item[ref]{Reference allele.}
#'                               \item[alt1]{Alternate allele 1.}
#'                               \item[alt]{All alternate alleles, comma-separated.}
#'                               \item[maj]{Major allele.}
#'                               \item[nonmaj]{All nonmajor alleles, comma-separated.}
#'                               \item[]{(PCs are always present, and positioned here. Signs are w.r.t. the major, not necessarily reference, allele.)}
#'                             }
#'                             The default with \code{var.wts = FALSE} is maybefid,maybesid.
#'                             With \code{var.wts = TRUE} is maybefid,maybesid,chrom,maj,nonmaj.
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#' 
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_pca <- function(input.prefix, output.prefix, var.wts, num.pcs, approx, meanimpute, vzs, output.cols, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(input.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (checkmate::test_file_exists(sprintf("%s.pgen", input.prefix))) {
    checkmate::assert_file(sprintf("%s.pgen", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.psam", input.prefix), access = "r", add = assertions)
    if (missing(vzs)) {
      input.prefix <- sprintf("--pfile %s", input.prefix)
      checkmate::assert_file(sprintf("%s.pvar", input.prefix), access = "r", add = assertions) 
    } else {
      checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
      checkmate::assert_file(sprintf("%s.pvar.zst", input.prefix), access = "r", add = assertions)
      input.prefix <- sprintf("--pfile %s vzs", input.prefix)
    } 
  } else {
    checkmate::assert_file(sprintf("%s.bed", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", input.prefix), access = "r", add = assertiosn)
    input.prefix <- sprintf("--bfile %s", input.prefix)
  }
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  
  if (!missing(num.pcs)) {
    checkmate::assert_int(num.pcs, na.ok = FALSE, lower = 1, null.ok = FALSE, add = assertions)
  } else {
    num.pcs <- ""
  }
  
  if (!missing(approx)) {
    checkmate::assert_flag(approx, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (approx) {
      approx <- "approx"
    } else {
      approx <- ""
    }
  } else {
    approx <- ""
  }
  
  if (!missing(meanimpute)) {
    checkmate::assert_flag(meanimpute, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (meanimpute) {
      meanimpute <- "meanimpute"
    } else {
      meanimpute <- ""
    }
  } else {
    meanimpute <- ""
  }
  
  if (!missing(var.wts)) {
    checkmate::assert_flag(var.wts, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (var.wts) {
      var.wts <- "var-wts"
      if (!missing(vzs)) {
        checkmate::assert_flag(vzs, na.ok = FALSE, null.ok = FALSE, add = assertions)
        if (vzs) {
          vzs <- "vzs"
        } else {
          vzs <- ""
        }
      }
    } else {
      var.wts <- ""
      vzs <- ""
    }
  } else {
    var.wts <- ""
    vzs <- ""
  }
  
  if (!missing(output.cols)) {
    checkmate::assert_character(output.cols, any.missing = FALSE, all.missing = FALSE, unique = TRUE, null.ok = FALSE, min.len = 1, add = assertions)
    if (var.wts == "var-wts") {
      checkmate::assert_subset(output.cols, choices = c("maybefid", "fid", "maybesid", "sid", "chrom", "pos", "ref", "alt", "alt1", "maj", "nonmaj"), add = assertions)
      vcols <- sprintf("vcols=%s", paste(output.cols[output.cols %in% c("chrom", "pos", "ref", "alt", "alt1", "maj", "nonmaj")], collapse = ","))
    } else {
      checkmate::assert_subset(output.cols, choices = c("maybefid", "fid", "maybesid", "sid"), add = assertions)
      vcols <- ""
    }
    scols <- sprintf("scols=%s", paste(output.cols[output.cols %in% c("maybefid", "fid", "maybesid", "sid")], collapse = ","))
  } else {
    scols <- ""
    vcols <- ""
  }
  
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      "--pca", var.wts, num.pcs, approx, meanimpute, vzs, scols, vcols,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}


#' Subset samples and variants using PLINK2
#'
#' @param input.prefix         [\code{string}]\cr
#'                             File path prefix of PLINK1 or PLINK2 binary file set.  
#' @param output.prefix        [\code{string}]\cr
#'                             File path prefix of output files.  
#' @param keep                 [\code{string}]\cr 
#'                             File name of a list of sample FIDs and IIDs to keep.
#' @param remove               [\code{string}]\cr 
#'                             File name of a list of sample FIDs and IIDs to remove. 
#' @param extract              [\code{string}]\cr 
#'                             File name of a list of variant names to keep. 
#' @param exclude              [\code{string}]\cr 
#'                             File name of a list of variant names to remove.  
#' @param ...                  Additional arguments passed directly to PLINK2 without any checks.
#' @param plink2.exec          [\code{string}]\cr
#'                             Path of PLINK2 executable. Defaults to \code{"plink2"}.
#' @param num.threads          [\code{int}]\cr
#'                             Number of CPUs usable by PLINK.
#'                             Default is determined by SLURM environment variables and at least 1.
#' @param memory               [\code{int}]\cr
#'                             Memory for PLINK in Mb.
#'                             Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#' 
#' @return Captured system output as \code{character} vector.
#' @export
#' 
#' @import checkmate imbs
#' 
plink2_subset <- function(input.prefix, output.prefix, keep, remove, extract, exclude, ..., num.threads, memory, plink2.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_string(input.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (checkmate::test_file_exists(sprintf("%s.pgen", input.prefix))) {
    checkmate::assert_file(sprintf("%s.pgen", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.psam", input.prefix), access = "r", add = assertions)
    make <- "--make-pgen"
  } else {
    checkmate::assert_file(sprintf("%s.bed", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.fam", input.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.bim", input.prefix), access = "r", add = assertiosn)
    input.prefix <- sprintf("--bfile %s", input.prefix)
    make <- "--make-bed"
  }
  checkmate::assert_string(output.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), access = "rwx", add = assertions)
  
  if (!missing(remove)) {
    checkmate::assert_file(remove, add = assertions)
    remove <- sprintf("--remove %s", remove)
  } else {
    remove <- ""
  }
  if (!missing(keep)) {
    checkmate::assert_file(keep, add = assertions)
    keep <- sprintf("--keep %s", keep)
  } else {
    keep <- ""
  }
  if (!missing(exclude)) {
    checkmate::assert_file(exclude, add = assertions)
    exclude <- sprintf("--exclude %s", exclude)
  } else {
    exclude <- ""
  }
  if (!missing(extract)) {
    checkmate::assert_file(extract, add = assertions)
    extract <- sprintf("--extract %s", extract)
  } else {
    extract <- ""
  }
  
  if (missing(plink2.exec)) {
    plink2.exec <- "plink2"
  }
  imbs::assert_command(plink2.exec, add = assertions)
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    plink2.exec,
    args = c(
      input.prefix,
      make,
      remove, keep, extract, exclude,
      "--out", output.prefix,
      "--memory", memory,
      "--threads", num.threads,
      ...
    )
  )
  
}

