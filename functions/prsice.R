
#' Calculate Polygenic Risk Scores with PRSice
#' 
#' Simple wrapper to call PRSice software from within R.
#'
#' @param base.file              [\code{string}]\cr
#'                               Path to summary statistics file. 
#'                               At least a variant identifier, statistic, effect allele and pvalue column are needed.
#' @param base.options           [\code{list}]\cr
#'                               A \code{list} of options regarding the handling of the summary statistics file.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{chr.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing the chromosome.
#'                                   Default: CHR.
#'                                 }
#'                                 \item{pos.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing the variant coordinate.
#'                                   Default: BP.
#'                                 }
#'                                 \item{snp.id.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing the variant identifier.
#'                                   Default: SNP.
#'                                 }
#'                                 \item{effect.allele.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing allele 1 (effective allele).
#'                                   Default: A1
#'                                 }
#'                                 \item{non.effect.allele.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing allele 2 (non-effective allele). 
#'                                   Default: A2
#'                                 }
#'                                 \item{stat.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing the summary statistic/effect.
#'                                   Default: BETA/OR
#'                                 }
#'                                 \item{pvalue.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing the p-value.
#'                                   Default: P.
#'                                 }
#'                                 \item{info.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing imputation info score of variants. 
#'                                   Default: INFO
#'                                 }
#'                                 \item{info.threshold}{
#'                                   [\code{number}]\cr
#'                                   Variants with info score less than \code{info.threshold} will be ignored. 
#'                                   Default: 0.9.
#'                                 }
#'                                 \item{maf.cols}{
#'                                   [\code{character}]\cr
#'                                   Column headers containing containing minor allele frequencies. 
#'                                   Multiple columns are allowed, e.g. to also filter MAF for cases.
#'                                 }
#'                                 \item{maf.thresholds}{
#'                                   [\code{numeric}]\cr
#'                                   Variants with MAF less than \code{maf.thresholds} in the corresponding \code{maf.cols} will be ignored. 
#'                                 }
#'                                 \item{beta}{
#'                                   [\code{flag}]\cr
#'                                   Whether the test statistic is in the form of BETA or OR. 
#'                                   If set, test statistic is assume to be in the form of BETA. 
#'                                   Overwritten by \code{or}.
#'                                 }
#'                                 \item{or}{
#'                                   [\code{flag}]\cr
#'                                   Whether the test statistic is in the form of BETA or OR. 
#'                                   If set, test statistic is assume to be in the form of BETA. 
#'                                   Overwrites \code{beta}.
#'                                 }
#'                                 \item{index}{
#'                                   [\code{flag}]\cr
#'                                   If set to \code{TRUE}, assume the INDEX instead of NAME for the corresponding columns are provided. 
#'                                   Index should be 0-based (start counting from 0).
#'                                 }
#'                               }
#' @param target.file.prefix     [\code{string}]\cr
#'                               Target genotype file. 
#'                               Currently supports both BGEN and binary PLINK format. 
#'                               For multiple chromosome input, simply substitute the chromosome number with #. 
#'                               PRSice will automatically replace # with 1-22.
#' @param target.options         [\code{list}]\cr
#'                               A \code{list} of options regarding the handling of the target file.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{binary}{
#'                                   [\code{flag}]\cr
#'                                   Indicate whether the target phenotype is binary or not. 
#'                                 }
#'                                 \item{geno}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on gentype missingnes. 
#'                                 }
#'                                 \item{maf}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on minor allele frequency. 
#'                                 }
#'                                 \item{info}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on imputation info score. 
#'                                 }
#'                                 \item{keep}{
#'                                   [\code{string}]\cr
#'                                   File containing the sample(s) to be extracted from the target file. 
#'                                   First column should be FID and the second column should be IID. 
#'                                   If \code{ignore.fid} is set to \code{TRUE}, first column should be IID.
#'                                 }
#'                                 \item{remove}{
#'                                   [\code{string}]\cr
#'                                   File containing the sample(s) to be removed from the target file. 
#'                                   First column should be FID and the second column should be IID. 
#'                                   If \code{ignore.fid} is set to \code{TRUE}, first column should be IID.
#'                                 }
#'                                 \item{type}{
#'                                   [\code{string}]\cr
#'                                   File type of the target file. Support \emph{bed} (binary plink) and \emph{bgen} format. 
#'                                   Default: bed
#'                                 }
#'                                 \item{pheno.file}{
#'                                   [\code{string}]\cr
#'                                   Phenotype file containing the phenotype(s). 
#'                                   First column must be FID of the samples and the second column must be IID of the samples. 
#'                                   When \code{ignore.fid} is set to \code{TRUE}, first column must be the IID of the samples. 
#'                                   Must contain a header if \code{pheno.col} is specified.
#'                                 }
#'                                 \item{pheno.col}{
#'                                   [\code{string}]\cr
#'                                   Column header containing phenotype information. 
#'                                 }
#'                                 \item{nonfounders}{
#'                                   [\code{flag}]\cr
#'                                   Keep the nonfounders in the analysis.
#'                                   Note: They will still be excluded from LD calculation.
#'                                 }
#'                               }
#' @param output.file.prefix     [\code{string}]\cr
#'                               Prefix for all file output.
#' @param ld.file.prefix         [\code{string}]\cr
#'                               LD reference file. Use for LD calculation. 
#'                               If not provided, will use the post-filtered target genotype for LD calculation. 
#'                               Support multiple chromosome input. Please see \code{target.file.prefix} for more information.
#' @param ld.file.options        [\code{list}]\cr
#'                               A \code{list} of options regarding the handling of the LD reference file.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{geno}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on genotype missingness. 
#'                                 }
#'                                 \item{maf}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on minor allele frequency. 
#'                                 }
#'                                 \item{info}{
#'                                   [\code{number}]\cr
#'                                   Filter variants based on imputation info score. 
#'                                 }
#'                                 \item{keep}{
#'                                   [\code{string}]\cr
#'                                   File containing the sample(s) to be extracted from the target file. 
#'                                   First column should be FID and the second column should be IID. 
#'                                   If \code{ignore.fid} is set to \code{TRUE}, first column should be IID.
#'                                 }
#'                                 \item{remove}{
#'                                   [\code{string}]\cr
#'                                   File containing the sample(s) to be removed from the target file. 
#'                                   First column should be FID and the second column should be IID. 
#'                                   If \code{ignore.fid} is set to \code{TRUE}, first column should be IID.
#'                                 }
#'                                 \item{type}{
#'                                   [\code{string}]\cr
#'                                   File type of the LD reference file. 
#'                                   Supports \emph{bed} (binary plink) and \emph{bgen} format. 
#'                                   Default: bed.
#'                                 }
#'                                 \item{dosage.threshold}{
#'                                   [\code{number}]\cr
#'                                   Translate any SNPs with highest genotype probability less than this threshold to missing call.
#'                                 }
#'                                 \item{hardcall.threshold}{
#'                                   [\code{number}]\cr
#'                                   A hardcall is saved when the distance to the nearest hardcall is less than the hardcall threshold. 
#'                                   Otherwise a missing code is saved.
#'                                   Default: 0.1.
#'                                 }
#'                               }
#' @param clumping               [\code{flag}]\cr
#'                               Clumping is disabled if set to \code{FALSE}.
#' @param clumping.options       [\code{list}]\cr
#'                               A \code{list} of options regarding the clumping procedure.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{kb}{
#'                                   [\code{int}]\cr
#'                                   The distance for clumping in kb.
#'                                   Default: 250.
#'                                 }
#'                                 \item{r2}{
#'                                   [\code{number}]\cr
#'                                   The \eqn{r^2} threshold for clumping.
#'                                   Default: 0.1.
#'                                 }
#'                                 \item{p}{
#'                                   [\code{number}]\cr
#'                                   The p-value threshold for clumping.
#'                                   Default: 1.
#'                                 }
#'                                 \item{proxy}{
#'                                   [\code{number}]\cr
#'                                   Proxy threshold for index SNP to be considered as part of the region represented by the clumped SNP(s). 
#'                                   E.g. 0.8 means the index SNP will represent region of any clumped SNP(s) that has a \eqn{r^2>=0.8} even if the index SNP does not physically locate within the region.
#'                                 }
#'                                 \item{pearson}{
#'                                   [\code{flag}]\cr
#'                                   Use Pearson Correlation for LD calculation instead of the maximum likelihood haplotype frequency estimates if set to \code{TRUE}. 
#'                                   This will slightly decrease the accuracy of LD estimates, but should increase the speed of clumping.
#'                                 }
#'                               }
#' @param dosage.options         [\code{list}]\cr
#'                               A \code{list} of options regarding the dosage conversion if no LD reference file is specified.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{allow.intermediate}{
#'                                   [\code{flag}]\cr
#'                                   Allow the generate of intermediate file. 
#'                                   This will speed up PRSice when using dosage data as clumping reference and for hard coding PRS calculation.
#'                                 }
#'                                 \item{dosage.threshold}{
#'                                   [\code{number}]\cr
#'                                   Translate any SNPs with highest genotype probability less than this threshold to missing call.
#'                                 }
#'                                 \item{hardcall.threshold}{
#'                                   [\code{number}]\cr
#'                                   A hardcall is saved when the distance to the nearest hardcall is less than the hardcall threshold. 
#'                                   Otherwise a missing code is saved.
#'                                   Default: 0.1.
#'                                 }
#'                                 \item{hardcall}{
#'                                   [\code{flag}]\cr
#'                                   Use hard coding instead of dosage for PRS construction.
#'                                   Default is to use dosage instead of hard coding.
#'                                 }
#'                               }
#' @param pval.options           [\code{list}]\cr
#'                               A \code{list} of options regarding the p-value thresholding.
#'                               The following options are allowed and checked:
#'                               \describe{
#'                                 \item{levels}{
#'                                   [\code{numeric}]\cr
#'                                   PRSice will only calculate the PRS for the given thresholds.
#'                                 }
#'                                 \item{full}{
#'                                   [\code{flag}]\cr
#'                                   By default, PRSice will include the full model,  i.e. p-value threshold = 1. Setting this flag to \code{TRUE} will disable that behaviour.
#'                                 }
#'                                 \item{lower}{
#'                                   [\code{number}]\cr
#'                                   The starting p-value threshold. 
#'                                   Default: 5e-08.
#'                                 }
#'                                 \item{upper}{
#'                                   [\code{number}]\cr
#'                                   The final p-value threshold. 
#'                                   Default: 0.5.
#'                                 }
#'                                 \item{step.size}{
#'                                   [\code{flag}]\cr
#'                                   The step size of the p-value thresholds. 
#'                                   Default: 5e-05.
#'                                 }
#'                                 \item{model}{
#'                                   [\code{string}]\cr
#'                                   Genetic model used for PRS calculation. 
#'                                   The genetic encoding is based on the target data where the encoding represent number of the coding allele.
#'                                   Available models include:
#'                                   \describe{
#'                                     \item[add]{
#'                                       Additive model, code as 0/1/2 (default).
#'                                     }
#'                                     \item[dom]{
#'                                       Dominant model, code as 0/1/1.
#'                                     }
#'                                     \item[rec]{
#'                                       Recessive model, code as 0/0/1.
#'                                     } 
#'                                     \item[het]{
#'                                      Heterozygous only model, code as 0/1/0.
#'                                     }
#'                                 }
#'                               }
#' @param missing.handling       [\code{string}]\cr
#'                               Method to handle missing genotypes. 
#'                               By default, final scores are averages of valid per-allele scores with missing genotypes contribute an amount proportional to imputed allele frequency (\code{IMPUTE}). 
#'                               To throw out missing observations instead (decreasing the denominator in the final average when this happens), use \code{SET_ZERO}.
#'                               Alternatively, you can use \code{CENTER} to shift all scores to mean zero. 
#' @param regress                [\code{flag}]\cr
#'                               Whether to perform the regression analysis (\code{TRUE}) or simply output all PRS (\code{FALSE}).
#' @param score                  [\code{string}]\cr
#'                               Method to calculate the polygenic score.
#'                               Available methods include:
#'                               \description{
#'                                 \item[avg]{
#'                                   Take the average effect size (default).
#'                                 }
#'                                 \item[std]{
#'                                   Standardize the effect size.
#'                                 }
#'                                 \item[sum]{
#'                                   Direct summation of the effect size.
#'                                 }
#'                               }
#' @param exclude                [\code{string}]\cr
#'                               Path to file which contains variants to be excluded from the analysis.
#' @param exclude.range          [\code{data.frame}]\cr
#'                               A \code{data.frame} with columns CHR, START, END.
#'                               Range of SNPs to be excluded from the whole analysis.
#' @param extract                [\code{string}]\cr
#'                               Path to file which contains variants to be included in the analysis.
#' @param id.delim               [\code{string}]\cr
#'                               This parameter causes sample IDs to be parsed as <FID><delimiter><IID>. 
#'                               Default: '_'. 
#' @param ignore.fid             [\code{flag}]\cr
#'                               If set to \code{TRUE} ignore FID for all input. 
#'                               When this is set, first column of all file will be assume to be IID instead of FID.
#' @param keep.ambig             [\code{flag}]\cr
#'                               Keep ambiguous SNPs. 
#'                               Only use this option if you are certain that the base and target has the same A1 and A2 alleles.
#' @param seed                   [\code{int}]\cr
#'                               Seed used for permutation. If not provided, system time will be used as seed. When same seed and same input is provided, same results can be generated.
#' @param print.snps             [\code{flag}]\cr
#'                               Print all SNPs used to construct the best PRS.
#' @param ...                    Additional arguments passed directly to the system call arguments.
#' @param memory                 [\code{int}]\cr
#'                               Memory for PRSice in Mb. 
#'                               Default is determined by minimum of SLURM environment variables SLURM_MEM_PER_NODE and num.threads * SLURM_MEM_PER_CPU and at least 5000.
#'                               PRSice will try its best to honor this setting but it's not guaranteed.
#' @param num.threads            [\code{int}]\cr
#'                               Number of threads to use.
#' @param exec                   [\code{string}]\cr
#'                               Path to PRSice executable.
#' 
#' @details Not all options are included. For full list of options call \emph{prsice -h}.
#' For a detailed description see \url{https://choishingwan.github.io/PRSice/}.
#'
#' @return The system call output as character vector.
#' 
#' @export
prsice <- function(
  base.file, base.options, 
  target.file.prefix, target.options, 
  output.file.prefix, 
  ld.file.prefix, ld.file.options,
  clumping, clumping.options, 
  dosage.options, 
  pval.options, 
  missing.handling,
  regress,
  score,
  exclude, exclude.range, extract, 
  id.delim, 
  ignore.fid, 
  keep.ambig, 
  seed,
  print.snps,
  ...,
  memory, num.threads, exec
) {
  
  assertions <- checkmate::makeAssertCollection()
  
  # Check base file options ----
  checkmate::assert_file(base.file, add = assertions) 
  if (!missing(base.options)) {
    checkmate::assert_list(base.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = FALSE, add = assertions)
    checkmate::assert_subset(names(base.options), choices = c("effect.allele.col", "non.effect.allele.col", "info.col", "info.threshold", "maf.cols", "maf.thresholds", "beta", "pos.col", "chr.col", "index", "or", "pvalue.col", "snp.id.col", "stat.col"), add = assertions)
    
    if (exists("index", where = base.options)) {
      checkmate::assert_flag(base.options[["index"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      index <- base.options[["index"]]
      if (index) {
        base.options[["index"]] <- "--index"
      } else {
        base.options[["index"]] <- ""
      }
    } else {
      index <- FALSE
      base.options[["index"]] <- ""
    }
    if (exists("effect.allele.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["effect.allele.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["effect.allele.col"]], na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
      }
      base.options[["effect.allele.col"]] <- sprintf("--A1 %s", base.options[["effect.allele.col"]])
    } else {
      base.options[["effect.allele.col"]] <- ""
    }
    if (exists("non.effect.allele.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["non.effect.allele.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["non.effect.allele.col"]], na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
      }
      base.options[["non.effect.allele.col"]] <- sprintf("--A2 %s", base.options[["non.effect.allele.col"]])
    } else {
      base.options[["non.effect.allele.col"]] <- ""
    }
    if (exists("info.threshold", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["info.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["info.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      checkmate::assert_number(base.options[["info.threshold"]], na.ok = FALSE, null.ok = FALSE, lower = 0, upper = 1, finite = TRUE, add = assertions)
      base.options[["info"]] <- sprintf("--base-info %s,%g", base.options[["info.col"]], base.options[["info.threshold"]])
    } else {
      base.options[["info"]] <- ""
    }
    base.options[["info.col"]] <- NULL
    base.options[["info.threshold"]] <- NULL
    if (exists("maf.thresholds", where = base.options)) {
      if (index) {
        checkmate::assert_integer(base.options[["maf.cols"]], any.missing = FALSE, all.missing = FALSE, min.len = 1, unique = TRUE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_character(base.options[["maf.cols"]], any.missing = FALSE, all.missing = FALSE, min.len = 1, min.chars = 1, null.ok = FALSE, add = assertions)
      }
      checkmate::assert_numeric(base.options[["maf.thresholds"]], null.ok = FALSE, any.missing = FALSE, all.missing = FALSE, len = length(base.options[["maf.cols"]]), lower = 0, upper = 1, finite = TRUE, add = assertions)
      base.options[["maf"]] <- sprintf("--base-maf %s", paste(base.options[["maf.cols"]], base.options[["maf.thresholds"]], sep = ",", collapse = ":"))
    } else {
      base.options[["maf"]] <- ""
    }
    base.options[["maf.cols"]] <- NULL
    base.options[["maf.thresholds"]] <- NULL
    if (exists("beta", where = base.options)) {
      checkmate::assert_flag(base.options[["beta"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (base.options[["beta"]]) {
        base.options[["scale"]] <- "--beta"
      }
    } else {
      base.options[["scale"]] <- ""
    }
    base.options[["beta"]] <- NULL
    if (exists("or", where = base.options)) {
      checkmate::assert_flag(base.options[["or"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (base.options[["or"]]) {
        base.options[["scale"]] <- "--or"
      }
    } else {
      base.options[["scale"]] <- ""
    }
    base.options[["or"]] <- NULL
    if (exists("pos.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["pos.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["pos.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      base.options[["pos.col"]] <- sprintf("--bp %s", base.options[["pos.col"]])
    } else {
      base.options[["pos.col"]] <- ""
    }
    if (exists("chr.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["chr.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["chr.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      base.options[["chr.col"]] <- sprintf("--chr %s", base.options[["chr.col"]])
    } else {
      base.options[["chr.col"]] <- ""
    }
    if (exists("pvalue.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["pvalue.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["pvalue.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      base.options[["pvalue.col"]] <- sprintf("--pvalue %s", base.options[["pvalue.col"]])
    } else {
      base.options[["pvalue.col"]] <- ""
    }
    if (exists("snp.id.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["snp.id.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["snp.id.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      base.options[["snp.id.col"]] <- sprintf("--snp %s", base.options[["snp.id.col"]])
    } else {
      base.options[["snp.id.col"]] <- ""
    }
    if (exists("stat.col", where = base.options)) {
      if (index) {
        checkmate::assert_int(base.options[["stat.col"]], na.ok = FALSE, lower = 0, null.ok = FALSE, add = assertions)
      } else {
        checkmate::assert_string(base.options[["stat.col"]], na.ok = FALSE, null.ok = FALSE, min.chars = 1, add = assertions)
      }
      base.options[["stat.col"]] <- sprintf("--stat %s", base.options[["stat.col"]])
    } else {
      base.options[["stat.col"]] <- ""
    }
  } else {
    base.options <- list()
  }
  
  # Check target file options ----
  if (!missing(target.options)) {
    checkmate::assert_list(target.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = FALSE, add = assertions)
    checkmate::assert_subset(names(target.options), choices = c("binary", "geno", "info", "keep", "maf", "nonfounders", "pheno.file", "pheno.col", "remove", "type"), add = assertions)
    
    if (exists("binary", where = target.options)) {
      checkmate::assert_flag(target.options[["binary"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (target.options[["binary"]]) {
        target.options[["binary"]] <- "--binary-target T"
      } else {
        target.options[["binary"]] <- ""
      }
    } else {
      target.options[["binary"]] <- ""
    }
    if (exists("nonfounders", where = target.options)) {
      checkmate::assert_flag(target.options[["nonfounders"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (target.options[["nonfounders"]]) {
        target.options[["nonfounders"]] <- "--nonfounders"
      } else {
        target.options[["nonfounders"]] <- ""
      }
    } else {
      target.options[["nonfounders"]] <- ""
    }
    if (exists("geno", where = target.options)) {
      checkmate::assert_number(target.options[["geno"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      target.options[["geno"]] <- sprintf("--geno %g", target.options[["geno"]])
    } else {
      target.options[["geno"]] <- ""
    }
    if (exists("info", where = target.options)) {
      checkmate::assert_number(target.options[["info"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      target.options[["info"]] <- sprintf("--geno %g", target.options[["info"]])
    } else {
      target.options[["info"]] <- ""
    }
    if (exists("maf", where = target.options)) {
      checkmate::assert_number(target.options[["maf"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      target.options[["maf"]] <- sprintf("--maf %g", target.options[["maf"]])
    } else {
      target.options[["maf"]] <- ""
    }
    if (exists("keep", where = target.options)) {
      checkmate::assert_file(target.options[["keep"]], access = "r", add = assertions)
      target.options[["keep"]] <- sprintf("--keep %s", target.options[["keep"]])
    } else {
      target.options[["keep"]] <- ""
    }
    if (exists("remove", where = target.options)) {
      checkmate::assert_file(target.options[["remove"]], access = "r", add = assertions)
      target.options[["remove"]] <- sprintf("--remove %s", target.options[["remove"]])
    } else {
      target.options[["remove"]] <- ""
    }
    if (exists("pheno.file", where = target.options)) {
      checkmate::assert_file(target.options[["pheno.file"]], access = "r", add = assertions)
      target.options[["pheno.file"]] <- sprintf("--pheno %s", target.options[["pheno.file"]])
    } else {
      target.options[["pheno.file"]] <- ""
    }
    if (exists("pheno.col", where = target.options)) {
      checkmate::assert_string(target.options[["pheno.col"]], na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
      target.options[["pheno.col"]] <- sprintf("--pheno-col %s", target.options[["pheno.col"]])
    } else {
      target.options[["pheno.col"]] <- ""
    }
    if (exists("type", where = target.options)) {
      checkmate::assert_choice(target.options[["type"]], choices = c("bed", "bgen"), null.ok = FALSE, add = assertions)
      target.type <- target.options[["type"]]
      target.options[["type"]] <- sprintf("--type %s", target.options[["type"]])
    } else {
      target.type <- "bed"
      target.options[["type"]] <- ""
    }
  } else {
    target.options <- list()
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
  
  # Check output file prefix ----
  checkmate::assert_string(output.file.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.file.prefix), add = assertions)
  
  # Check LD file options ----
  if (!missing(ld.file.prefix)) {
    checkmate::assert_string(ld.file.prefix, null.ok = FALSE, na.ok = FALSE, add = assertions)
    if (!missing(ld.file.options)) {
      checkmate::assert_list(ld.file.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = FALSE, add = assertions)
      checkmate::assert_subset(names(ld.file.options), choices = c("dosage.threshold", "geno", "hardcall.threshold", "info", "keep", "remove", "maf", "type"), add = assertions)
      
      if (exists("dosage.threshold", where = ld.file.options)) {
        checkmate::assert_number(ld.file.options[["dosage.threshold"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        ld.file.options[["dosage.threshold"]] <- sprintf("--ld-dose-thres %g", ld.file.options[["dosage.threshold"]])
      } else {
        ld.file.options[["dosage.threshold"]] <- ""
      }
      if (exists("geno", where = ld.file.options)) {
        checkmate::assert_number(ld.file.options[["geno"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        ld.file.options[["geno"]] <- sprintf("--ld-geno %g", ld.file.options[["geno"]])
      } else {
        ld.file.options[["geno"]] <- ""
      }
      if (exists("hardcall.theshold", where = ld.file.options)) {
        checkmate::assert_number(ld.file.options[["hardcall.theshold"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        ld.file.options[["hardcall.theshold"]] <- sprintf("--ld-hard-thres %g", ld.file.options[["hardcall.theshold"]])
      } else {
        ld.file.options[["hardcall.theshold"]] <- ""
      }
      if (exists("info", where = ld.file.options)) {
        checkmate::assert_number(ld.file.options[["info"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        ld.file.options[["info"]] <- sprintf("--ld-info %g", ld.file.options[["info"]])
      } else {
        ld.file.options[["info"]] <- ""
      }
      if (exists("keep", where = ld.file.options)) {
        checkmate::assert_file(ld.file.options[["keep"]], access = "r", null.ok = FALSE, add = assertions)
        ld.file.options[["keep"]] <- sprintf("--ld-keep %s", ld.file.options[["keep"]])
      } else {
        ld.file.options[["keep"]] <- ""
      }
      if (exists("remove", where = ld.file.options)) {
        checkmate::assert_file(ld.file.options[["remove"]], access = "r", null.ok = FALSE, add = assertions)
        ld.file.options[["remove"]] <- sprintf("--ld-remove %s", ld.file.options[["remove"]])
      } else {
        ld.file.options[["remove"]] <- ""
      }
      if (exists("maf", where = ld.file.options)) {
        checkmate::assert_number(ld.file.options[["maf"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        ld.file.options[["maf"]] <- sprintf("--ld-maf %g", ld.file.options[["maf"]])
      } else {
        ld.file.options[["maf"]] <- ""
      }
      if (exists("type", where = ld.file.options)) {
        checkmate::assert_choice(ld.file.options[["type"]], choices = c("bed", "bgen"), null.ok = FALSE, add = assertions)
        ld.type <- ld.file.options[["type"]]
        ld.file.options[["type"]] <- sprintf("--ld-type %s", ld.file.options[["type"]])
      } else {
        ld.type <- "bed"
        ld.file.options[["type"]] <- ""
      }
    } else {
      ld.file.options <- list()
      ld.type <- "bed"
    }
    if (ld.type == "bed") {
      checkmate::assert_file(sprintf("%s.bed", ld.file.prefix), access = "r", add = assertions)
      checkmate::assert_file(sprintf("%s.bim", ld.file.prefix), access = "r", add = assertions)
      checkmate::assert_file(sprintf("%s.fam", ld.file.prefix), access = "r", add = assertions)
    } else {
      checkmate::assert_file(sprintf("%s.bgen", ld.file.prefix), access = "r", add = assertions)
      checkmate::assert_file(sprintf("%s.bgen.bgi", ld.file.prefix), access = "r", add = assertions)
    }
    ld.file <- sprintf("--ld %s", ld.file.prefix)
  } else {
    ld.file <- ""
    ld.file.options <- list()
  }
  
  # Check clumping options ----
  checkmate::assert_flag(clumping, na.ok = FALSE, null.ok = FALSE, add = assertions)
  if (!clumping) {
    clumping <- "--no-clump"
    clumping.options <- list()
  } else {
    clumping <- ""
    if (!missing(clumping.options)) {
      checkmate::assert_list(clumping.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = FALSE, add = assertions)
      checkmate::assert_subset(names(clumping.options), choices = c("kb", "r2", "p", "proxy", "pearson"), add = assertions)
      
      if (exists("kb", where = clumping.options)) {
        checkmate::assert_int(clumping.options[["kb"]], na.ok = FALSE, lower = 1, null.ok = FALSE, add = assertions)
        clumping.options[["kb"]] <- sprintf("--clump-kb %g", clumping.options[["kb"]])
      } else {
        clumping.options[["kb"]] <- ""
      }
      if (exists("r2", where = clumping.options)) {
        checkmate::assert_number(clumping.options[["r2"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        clumping.options[["r2"]] <- sprintf("--clump-r2 %g", clumping.options[["r2"]])
      } else {
        clumping.options[["r2"]] <- ""
      }
      if (exists("p", where = clumping.options)) {
        checkmate::assert_number(clumping.options[["p"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        clumping.options[["p"]] <- sprintf("--clump-p %g", clumping.options[["p"]])
      } else {
        clumping.options[["p"]] <- ""
      }
      if (exists("proxy", where = clumping.options)) {
        checkmate::assert_number(clumping.options[["proxy"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
        clumping.options[["proxy"]] <- sprintf("--proxy %g", clumping.options[["proxy"]])
      } else {
        clumping.options[["proxy"]] <- ""
      }
      if (exists("pearson", where = clumping.options)) {
        checkmate::assert_flag(clumping.options[["pearson"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
        if (clumping.options[["pearson"]]) {
          clumping.options[["pearson"]] <- "--pearson"
        } else {
          clumping.options[["pearson"]] <- ""
        }
      } else {
        clumping.options[["pearson"]] <- ""
      }
    } else {
      clumping.options <- list()
    }
  }
  
  # Check dosage options ----
  if (!missing(dosage.options)) {
    checkmate::assert_list(dosage.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = TRUE, add = assertions)
    checkmate::assert_subset(names(dosage.options), choices = c("allow.intermediate", "dosage.threshold", "hardcall.threshold", "hardcall"), add = assertions)
    
    if (exists("allow.intermediate", where = dosage.options)) {
      checkmate::assert_flag(dosage.options[["allow.intermediate"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (dosage.options[["allow.intermediate"]]) {
        dosage.options[["allow.intermediate"]] <- "--allow.intermediate"
      } else {
        dosage.options[["allow.intermediate"]] <- ""
      }
    } else {
      dosage.options[["allow.intermediate"]] <- ""
    }
    if (exists("dosage.threshold", where = dosage.options)) {
      checkmate::assert_number(dosage.options[["dosage.threshold"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      dosage.options[["dosage.threshold"]] <- sprintf("--dose-thres %g", dosage.options[["dosage.threshold"]])
    } else {
      dosage.options[["dosage.threshold"]] <- ""
    }
    if (exists("hardcall.theshold", where = dosage.options)) {
      checkmate::assert_number(dosage.options[["hardcall.theshold"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      dosage.options[["hardcall.theshold"]] <- sprintf("--hard-thres %g", dosage.options[["hardcall.theshold"]])
    } else {
      dosage.options[["hardcall.theshold"]] <- ""
    }
    if (exists("hardcall", where = dosage.options)) {
      checkmate::assert_flag(dosage.options[["hardcall"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (dosage.options[["hardcall"]]) {
        dosage.options[["hardcall"]] <- "--hard"
      } else {
        dosage.options[["hardcall"]] <- ""
      }
    } else {
      dosage.options[["hardcall"]] <- ""
    }
  } else {
    dosage.options <- list()
  }
  
  # Check pval options ----
  if (!missing(pval.options)) {
    checkmate::assert_list(pval.options, types = c("character", "int", "integer", "number", "numeric", "logical"), names = "unique", any.missing = TRUE, all.missing = TRUE, null.ok = FALSE, add = assertions)
    checkmate::assert_subset(names(pval.options), choices = c("levels", "full", "step.size", "lower", "model", "upper"), add = assertions)
    
    if (exists("levels", where = pval.options)) {
      checkmate::assert_numeric(pval.options[["levels"]], lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, all.missing = FALSE, min.len = 1, unique = TRUE, null.ok = FALSE, add = assertions)
      pval.options[["levels"]] <- sprintf("--bar-levels %s", paste(sprintf("%g", pval.options[["levels"]]), collapse = ","))
      pval.options[["fastscore"]] <- "--fastscore"
    } else {
      pval.options[["levels"]] <- ""
      pval.options[["fastscore"]] <- ""
    }
    if (exists("full", where = pval.options)) {
      checkmate::assert_flag(pval.options[["full"]], na.ok = FALSE, null.ok = FALSE, add = assertions)
      if (pval.options[["full"]]) {
        pval.options[["full"]] <- ""
      } else {
        pval.options[["full"]] <- "--no-full"
      }
    } else {
      pval.options[["full"]] <- ""
    }
    if (exists("step.size", where = pval.options)) {
      checkmate::assert_number(pval.options[["step.size"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      pval.options[["step.size"]] <- sprintf("--interval %g", pval.options[["step.size"]])
    } else {
      pval.options[["step.size"]] <- ""
    }
    if (exists("lower", where = pval.options)) {
      checkmate::assert_number(pval.options[["lower"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      pval.options[["lower"]] <- sprintf("--lower %g", pval.options[["lower"]])
    } else {
      pval.options[["lower"]] <- ""
    }
    if (exists("upper", where = pval.options)) {
      checkmate::assert_number(pval.options[["upper"]], na.ok = FALSE, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
      pval.options[["upper"]] <- sprintf("--upper %g", pval.options[["upper"]])
    } else {
      pval.options[["upper"]] <- ""
    }
    if (exists("model", where = pval.options)) {
      checkmate::assert_choice(pval.options[["model"]], choices = c("add", "dom", "rec", "het"), null.ok = FALSE, add = assertions)
      pval.options[["model"]] <- sprintf("--model %s", pval.options[["model"]])
    } else {
      pval.options[["model"]] <- ""
    }
  } else {
    pval.options <- list()
  }
  
  # Check misc options ----
  if (!missing(missing.handling)) {
    checkmate::assert_choice(missing.handling, choices = c("IMPUTE", "SET_ZERO", "CENTER"), null.ok = FALSE, add = assertions)
    if (missing.handling != "IMPUTE") {
      missing.handling <- sprintf("--missing %s", missing.handling)
    }
  } else {
    missing.handling <- ""
  }
  if (!missing(regress)) {
    checkmate::assert_flag(regress, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (regress) {
      regress <- ""
    } else {
      regress <- "--no-regress"
    }
  } else {
    regress <- ""
  }
  if (!missing(score)) {
    checkmate::assert_choice(score, choices = c("avg", "std", "sum"), null.ok = FALSE, add = assertions)
    score <- sprintf("--score %s", score)
  }
  if (!missing(exclude)) {
    checkmate::assert_file(exclude, access = "r", add = assertions)
    exclude <- sprintf("--exclude %s", exclude)
  } else {
    exclude <- ""
  }
  if (!missing(extract)) {
    checkmate::assert_file(extract, access = "r", add = assertions)
    extract <- sprintf("--extract %s", extract)
  } else {
    extract <- ""
  }
  if (!missing(exclude.range)) {
    checkmate::assert_data_frame(exclude.range, types = c("integer", "numeric"), any.missing = FALSE, all.missing = FALSE, min.rows = 1, min.cols = 3, col.names = "unique", null.ok = FALSE, add = assertions)
    checkmate::assert_set_equal(colnames(exclude.range), c("CHR", "START", "END"), add = assertions)
    exclude.range <- sprintf("--x-range %s", paste(sprintf("%s:%s-%s", exclude.range[["CHR"]], exclude.range[["START"]], exclude.range[["END"]]), collapse = ","))
  } else {
    exclude.range <- ""
  }
  if (!missing(id.delim)) {
    checkmate::assert_string(id.delim, na.ok = FALSE, min.chars = 1, null.ok = FALSE, add = assertions)
    id.delim <- sprintf("--id-delim %s", id.delim)
  } else {
    id.delim <- ""
  }
  if (!missing(ignore.fid)) {
    checkmate::assert_flag(ignore.fid, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (ignore.fid) {
      ignore.fid <- "--ignore-fid"
    } else {
      ignore.fid <- ""
    }
  } else {
    ignore.fid <- ""
  }
  if (!missing(keep.ambig)) {
    checkmate::assert_flag(keep.ambig, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (keep.ambig) {
      keep.ambig <- "--keep-ambig"
    } else {
      keep.ambig <- ""
    }
  } else {
    keep.ambig <- ""
  }
  if (!missing(seed)) {
    checkmate::assert_int(seed, na.ok = FALSE, null.ok = FALSE, add = assertions)
    seed <- sprintf("--seed %d", seed)
  } else {
    seed <- ""
  }
  if (!missing(print.snps)) {
    checkmate::assert_flag(print.snps, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (print.snps) {
      print.snps <- "--print-snp"
    } else {
      print.snps <- ""
    }
  } else {
    print.snps <- ""
  }
  
  checkmate::assert_string(exec, na.ok = FALSE, null.ok = FALSE, add = assertions)
  exec <- path.expand(exec)
  imbs::assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  threads <- sprintf("--thread %d", num.threads)
  
  if (missing(memory)) {
    memory <- BBmisc::suppressAll(max(5000, -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), na.rm = TRUE))
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  memory <- sprintf("--memory %dMb", memory)
  
  checkmate::reportAssertions(assertions)
  
  # Call PRSice ----
  imbs::system_call(
    bin = exec, 
    args = c(
      "--base", base.file, unlist(base.options), 
      "--target", target.file.prefix, unlist(target.options), 
      "--out", output.file.prefix, 
      ld.file, unlist(ld.file.options),
      clumping, unlist(clumping.options), 
      unlist(dosage.options), 
      unlist(pval.options), 
      missing.handling,
      regress,
      score,
      exclude, exclude.range, extract, 
      id.delim, 
      ignore.fid, 
      keep.ambig, 
      seed,
      print.snps,
      ...,
      memory, threads
    )
  )
  
}
