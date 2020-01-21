
#' Phasing with SHAPEIT2
#'
#' @param input.file.prefix        [\code{string}]\cr
#'                                 File path without extension to dataset to be phased.
#' @param input.type               [\code{string}]\cr
#'                                 Input file type. Can be any of \code{"PLINK_PED"}, \code{"PLINK_BED"}, \code{"VCF"}, or \code{"GEN"}. If \code{"VCF"} or \code{"GEN"}, input files may be gz-compressed.
#' @param map.file                 [\code{string}]\cr
#' @param output.file.prefix       File path to genetic map.
#' @param missing.code             [\code{string}]\cr
#'                                 Code for missing genotypes.
#' @param hardcall.threshold       [\code{number}]\cr
#'                                 Threshold for hard calling genotypes when GEN/SAMPLE files are given as input. 
#'                                 For each individual at each SNP, the program will assume that the most likely genotype is true if its probability exceeds the threshold.
#'                                 Otherwise, the genotype will be considered as missing. Default is \eqn{0.9}.
#' @param reference.file.prefix    [\code{string}]\cr
#'                                 File path without extension to optional reference dataset in IMPUTE2 format.
#' @param no.mcmc                  [\code{flag}]\cr
#'                                 Disables the MCMC iterations. 
#'                                 This needs a reference panel to be specified.
#' @param chrx                     [\code{flag}]\cr
#'                                 Specifies to SHAPEIT2 that the genotypes to be phased come from the non pseudo autosomal region of the X chromosome. 
#'                                 SHAPEIT2 looks therefore at the sex of each individual to determine the ploidy model to use. 
#' @param noped                    [\code{flag}]\cr
#'                                 Discard all pedigree/family information in the input files. 
#'                                 SHAPEIT2 treats all samples as unrelated when this flag is specified.
#' @param replicates               [\code{int}]\cr
#'                                 Number of haplotype pairs sampled per individual. Default is 1,000.
#' @param burn                     [\code{int}]\cr
#'                                 The number of burn-in iterations used by the algorithm to reach a good starting point. 
#'                                 The default is to use \eqn{7} burn-in iterations.
#' @param prune                    [\code{int}]\cr
#'                                 The number of pruning iterations used by the algorithm to find a parsimonious graph for each individual. 
#'                                 The default is to use \eqn{8} pruning iterations.
#' @param main                     [\code{int}]\cr
#'                                 The number of iterations used by the algorithm to compute transition probabilities in the haplotype graphs. 
#'                                 The default is to use \eqn{20} pruning iterations.
#' @param seed                     [\code{int}]\cr
#'                                 The seed of the random number generator. 
#'                                 Default is the value given by stdlib function time(NULL).
#' @param states                   [\code{int}]\cr
#'                                 To specify the number of conditioning haplotypes to be used during the phasing process. 
#'                                 The running time of the algorithm increases only linearly with this number. 
#'                                 The default value is 100. 
#' @param window                   [\code{number}]\cr
#'                                 Mean size of the windows in Mb in which conditioning haplotypes are defined. 
#'                                 This parameter has impact on the accuracy. 
#'                                 The default value is 2 which is ok for most GWAS application. 
#'                                 Use a 0.5 value when dealing with genotypes derived from sequencing. 
#' @param effective.size           [\code{int}]\cr
#'                                 The effective size of the population you want to phase. 
#'                                 Use the following values:
#'                                 \describe{
#'                                  \item[European]{11418}
#'                                  \item[African]{17469}
#'                                  \item[Asian]{14269}
#'                                  \item[Mixed]{between 11418 and 17469, depending on the proportion of each population in the dataset.}
#'                                 }
#' @param force                    [\code{flag}]\cr
#'                                 If set to \code{TRUE}, don't check missing rates etc.
#' @param exec                     [\code{string}]\cr
#'                                 Path of SHAPEIT2 executable.
#' @param num.threads              [\code{int}]\cr
#'                                 Number of CPUs usable by SHAPEIT2
#'                                 Default is determined by SLURM environment variables and at least 1.
#'                                 
#' @details Simple wrapper to SHAPEIT2
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors.
#' @export
#' 
#' @import checkmate imbs
#' 
shapeit <- function(input.file.prefix, input.type, map.file, output.file.prefix, missing.code, hardcall.threshold, reference.file.prefix, no.mcmc, chrx, noped, replicates, burn, prune, main, seed, states, window, effective.size, force, exec, num.threads) {
  
  assertions <- checkmate::makeAssertCollection()
  
  input.types <- list(
    PED_MAP = c("ped", "map"),
    VCF = c("vcf.gz", "vcf"),
    PLINK_BED = c("bed", "bim", "fam"),
    GEN = c("gen", "sample")
  )
  
  # Basic I/O arguments ----
  checkmate::assert_choice(input.type, choices = names(input.types), null.ok = FALSE, add = assertions)
  checkmate::assert_string(input.file.prefix, add = assertions)
  switch (input.type,
          PED_MAP = {
            checkmate::assert_file(sprintf("%s.ped", input.file.prefix), access = "r", add = assertions)
            checkmate::assert_file(sprintf("%s.map", input.file.prefix), access = "r", add = assertions)
            input <- sprintf("--input-ped %s", input.file.prefix)
            input
          },
          VCF = {
            checkmate::assert_true(
              any(
                checkmate::test_file(sprintf("%s.vcf", input.file.prefix), access = "r"),
                checkmate::test_file(sprintf("%s.vcf.gz", input.file.prefix), access = "r")
              ),
              add = assertions
            )
            input <- sprintf("--input-vcf %s", input.file.prefix)
            input
          },
          PLINK_BED = {
            checkmate::assert_file(sprintf("%s.bed", input.file.prefix), access = "r", add = assertions)
            checkmate::assert_file(sprintf("%s.bim", input.file.prefix), access = "r", add = assertions)
            checkmate::assert_file(sprintf("%s.fam", input.file.prefix), access = "r", add = assertions)
            input <- sprintf("--input-bed %s", input.file.prefix)
            input
          },
          GEN = {
            checkmate::assert_true(
              any(
                all(
                  checkmate::test_file(sprintf("%s.gen", input.file.prefix), access = "r"),
                  checkmate::test_file(sprintf("%s.sample", input.file.prefix), access = "r"),
                ),
                all(
                  checkmate::test_file(sprintf("%s.gen.gz", input.file.prefix), access = "r"),
                  checkmate::test_file(sprintf("%s.sample", input.file.prefix), access = "r"),
                )
              ),
              add = assertions
            )
            input <- sprintf("--input-gen %s", input.file.prefix)
            input
          }
  )
  checkmate::assert_file_exists(map.file, access = "r", add = assertions)
  checkmate::assert_string(output.file.prefix, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.file.prefix), access = "rwx", add = assertions)
  
  if (!missing(missing.code)) {
    checkmate::assert_string(missing.code, na.ok = FALSE, null.ok = FALSE, add = assertions)
    missing.code <- sprintf("--missing-code %s", missing.code)
  } else {
    missing.code <- ""
  }
  
  if (!missing(hardcall.threshold)) {
    checkmate::assert_number(hardcall.threshold, na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    hardcall.threshold <- sprintf("--input-thr %f", hardcall.threshold)
  } else {
    hardcall.threshold <- ""
  }
  
  if (!missing(reference.file.prefix)) {
    checkmate::assert_string(reference.file.prefix, add = assertions)
    checkmate::assert_file(sprintf("%s.haps", reference.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.legend", reference.file.prefix), access = "r", add = assertions)
    checkmate::assert_file(sprintf("%s.sample", reference.file.prefix), access = "r", add = assertions)
    reference <- sprintf("--input-ref %.haps %.legend %.sample", reference.file.prefix, reference.file.prefix, reference.file.prefix)
    if (!missing(no.mcmc)) {
      checkmate::assert_flag(no.mcmc, add = assertions)
      no.mcmc <- "--no-mcmc"
    } else {
      no.mcmc <- ""
    }
  } else {
    reference <- ""
    no.mcmc <- ""
  }
  
  if (!missing(chrx)) {
    checkmate::assert_flag(chrx, add = assertions)
    chrx <- "--chrX"
  } else {
    chrx <- ""
  }
  
  if (!missing(noped)) {
    checkmate::assert_flag(noped, add = assertions)
    noped <- "--noped"
  } else {
    noped <- ""
  }
  
  if (!missing(replicates)) {
    checkmate::assert_int(replicates, na.ok = FALSE, lower = 1, null.ok = FALSE, add = assertions)
    replicates <- sprintf("--replicate %d", replicates)
  } else {
    replicates <- ""
  }
  
  if (!missing(burn)) {
    checkmate::assert_int(burn, lower = 0, null.ok = FALSE, add = assertions)
    burn <- sprintf("--burn %d", burn) 
  } else {
    burn <- ""
  }
  if (!missing(prune)) {
    checkmate::assert_int(prune, lower = 0, null.ok = FALSE, add = assertions)
    prune <- sprintf("--prune %d", prune) 
  } else {
    prune <- ""
  }
  if (!missing(main)) {
    checkmate::assert_int(main, lower = 0, null.ok = FALSE, add = assertions)
    main <- sprintf("--main %d", main) 
  } else {
    main <- ""
  }
  
  if (!missing(seed)) {
    checkmate::assert_int(seed, na.ok = FALSE, null.ok = FALSE, lower = 1, add = assertions)
    seed <- sprintf("--seed %d", seed)
  } else {
    seed <- ""
  }
  
  if (!missing(states)) {
    checkmate::assert_int(states, na.ok = FALSE, null.ok = FALSE, lower = 1, add = assertions)
    states <- sprintf("--states %d", states)
  } else {
    states <- ""
  }
  
  if (!missing(window)) {
    checkmate::assert_number(window, na.ok = FALSE, null.ok = FALSE, lower = 1, add = assertions)
    window <- sprintf("--window %d", window)
  } else {
    window <- ""
  }
  
  if (!missing(effective.size)) {
    checkmate::assert_int(effective.size, na.ok = FALSE, null.ok = FALSE, lower = 1, add = assertions)
    effective.size <- sprintf("--effective-size %d", effective.size)
  } else {
    effective.size <- ""
  }
  
  if (!missing(force)) {
    checkmate::assert_flag(force, add = assertions)
    force <- "--force"
  } else {
    force <- ""
  }
  
  imbs::assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c(
      input,
      "-O", output.file.prefix,
      sprintf("--output-log %s.log", output.file.prefix),
      "-M", map.file,
      missing.code,
      hardcall.threshold, 
      reference, 
      no.mcmc,
      chrx,
      noped, 
      replicates,
      burn, 
      prune,
      main,
      seed,
      states, 
      window, 
      effective.size, 
      force,
      sprintf("--thread %d", num.threads)
    )
  )
  
}