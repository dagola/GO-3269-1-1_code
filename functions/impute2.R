
#' Imputation with IMPUTE2
#'
#' @param study.haplotypes.file             [\code{string}]\cr
#'                                          File containing known haplotypes for the study cohort.
#' @param genetic.map.file                  [\code{string}]\cr
#'                                          Fine-scale recombination map for the region to be analyzed. 
#'                                          This file should have three columns: physical position (in base pairs), recombination rate between current position and next position in map (in cM/Mb), and genetic map position (in cM). 
#'                                          The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"). 
#' @param interval                          [\code{integer}]\cr
#'                                          Genomic interval to use for inference, as specified by this integer vector of length two in base pair position. 
#'                                          The boundaries can be expressed either in long form (e.g., -int 5420000 10420000) or in exponential notation (e.g., -int 5.42e6 10.42e6). 
#'                                          This option is particularly useful for restricting test jobs to small regions or splitting whole-chromosome analyses into manageable chunks. 
#' @param reference.haplotypes.file         [\code{string}]\cr
#'                                          File of known haplotypes (reference), with one row per SNP and one column per haplotype. 
#'                                          All alleles must be coded as 0 or 1, and this file must be provided with a corresponding legend file.
#' @param reference.legend.file             [\code{string}]\cr
#'                                          Legend file(s) with information about the SNPs in the -h file(s). 
#'                                          Each file should have four columns: rsID, physical position (in base pairs), allele 0, and allele 1. 
#'                                          The last two columns specify the alleles underlying the 0/1 coding in the corresponding haplotypes file; these alleles can take values in {A,C,G,T}. 
#'                                          Each legend file should also have a header line with an unbroken character string for each column (e.g., "rsID position a0 a1").
#' @param output.file                       [\code{string}]\cr
#'                                          Output file path. Additional output files will be created with '<output.file>_[suffix]'.
#' @param output.snps                       [\code{integer}]\cr
#'                                          "Output SNPs": specifies the SNP types that will be printed to the output file. 
#'                                          By default, all imputed and genotyped SNPs are included in the output, i.e., \code{c(0,1,2,3)}.
#' @param gz                                [\code{flag}]\cr
#'                                          Specifies that the main output file should be compressed by the gzip utility; this also applies to some non-standard output files that can become large.
#' @param output.decimal.places             [\code{int}]\cr
#'                                          Specifies the number of decimal places to use for reporting genotype probabilities in the main output file.
#' @param snp.qc.info                       [\code{flag}]\cr
#'                                          Suppresses printing of info_typeX, concord_typeX, and r2_typeX columns in the info file if set to \code{FALSE}.
#' @param sample.qc.info                    [\code{flag}]\cr
#'                                          Suppresses printing of per-sample quality control metrics file if set to \code{FALSE}.
#' @param predict.genotyped.snps            [\code{flag}]\cr
#'                                          Tells IMPUTE2 to replace the input genotypes from the \code{study.haplotypes.file} with imputed genotypes (applies to Type 2 SNPs only).
#' @param predict.missing.genotyped.snps    [\code{flag}]\cr
#'                                          Unlike \code{predict.genotyped.snps}, which replaces all input genotypes with imputed genotypes, this option tells IMPUTE2 to replace only the missing genotypes at typed SNPs. 
#'                                          That is, any non-missing genotype will simply be reprinted in the output file, whereas input genotypes that have missings will be imputed in the output.
#'                                          
#'                                          WARNING: This is an appealing option that will "fill in" sporadically missing genotypes in your input data. 
#'                                          However, it is possible that this could cause subtle problems in downstream association testing. 
#'                                          It is therefore suggested that you use caution when applying this option.
#' @param buffer                            [\code{int}]\cr
#'                                          Length of buffer region (in kb) to include on each side of the analysis interval specified by the -int option. SNPs in the buffer regions inform the inference but do not appear in output files (unless you activate the \code{include.buffer.in.output} flag). 
#'                                          
#'                                          Using a buffer region helps prevent imputation quality from deteriorating near the edges of the analysis interval. 
#'                                          Larger buffers may improve accuracy for low-frequency variants (since such variants tend to reside on long haplotype backgrounds) at the cost of longer running times.
#' @param allow.large.regions               [\code{flag}]\cr
#'                                          Allows the analysis of regions larger than 7 Mb if set to \code{TRUE}. 
#'                                          If this flag is not set to \code{TRUE} and the analysis interval plus buffer region exceeds 7 Mb, IMPUTE2 will quit with an error.
#' @param include.buffer.in.output          [\code{flag}]\cr
#'                                          Tells IMPUTE2 to include SNPs from the buffer region in all output files if set to \code{TRUE}. 
#' @param effective.size                    [\code{int}]\cr
#'                                          "Effective size" of the population (commonly denoted as Ne in the population genetics literature) from which your dataset was sampled.
#'                                          This parameter scales the recombination rates that IMPUTE2 uses to guide its model of linkage disequilibrium patterns. 
#'                                          When most imputation runs were conducted with reference panels from HapMap Phase 2, it was suggested to use values of \eqn{11418} for imputation from HapMap CEU, \eqn{17469} for YRI, and \eqn{14269} for CHB+JPT.
#'                                          Modern imputation analyses typically involve reference panels with greater ancestral diversity, which can make it hard to determine the "ideal" \code{effective.size} value for a particular study. 
#'                                          Fortunately, imputation accuracy is highly robust to different \code{effective.size} values; within each of several human populations, identical accuracy levels for values between \eqn{10000} and \eqn{25000} have been obtained. 
#'                                          It is suggested to set \code{effective.size} to \eqn{20000} in the majority of modern imputation analyses.
#' @param verbose                           [\code{flag}]\cr
#'                                          Print detailed output about the progress of imputation if set to \code{TRUE}. 
#'                                          By default, IMPUTE2 prints only the number of the current MCMC iteration when performing imputation, but this flag tells it to print more detailed updates.
#' @param filter.rules                      [\code{character}]\cr
#'                                          This option provides flexible variant filtering in the reference panel via "filter rules", which are based on annotation columns in a legend file. 
#'                                          Each column should be labeled by a contiguous string (no whitespace) describing its contents. 
#'                                          
#'                                          For example, a that contains columns named eur.maf and afr.maf and TYPE. 
#'                                          To filter variants based on the numeric annotation values in the legend file, you should combine a column string with a cutoff value and one of these six comparison operators: < <= > >= == != . 
#'                                          For example, giving \code{filter.rules = c("'eur.maf<0.05'")} would tell IMPUTE2 to remove any variants with eur.maf values less than 0.05 from the reference panel. 
#'                                          You can include an arbitrary number of filtering strings, in which case the filtering conditions will be applied in 'or' fashion: if any condition is true, the variant will be removed. 
#'                                          
#'                                          It is very important that you enclose each filtering string in single quotes, as shown above. 
#'                                          Otherwise, the command-line environment may interpret symbols like < and > as linux redirection operators. 
#'                                          There should be no white space within the single quotes.
#' @param ref.template.haplotypes           [\code{int}]\cr
#'                                          Number of reference haplotypes to use as templates when imputing missing genotypes. 
#'                                          
#'                                          As a rule of thumb, you should set \code{ref.template.haplotypes} to the number of reference haplotypes that you expect to be useful for your study population. 
#'                                          If this value is less than the total number of haplotypes in your reference panel, IMPUTE2 will choose a "custom" set of \code{ref.template.haplotypes} haplotypes each time it imputes missing alleles in a study haplotype.  
#'                                          If all of your reference haplotypes have similar ancestry to the subjects in your study, each haplotype is potentially useful for imputation, so the best accuracy can be achieved by setting \code{ref.template.haplotypes} to the total number of reference haplotypes. 
#'                                          Using smaller values will decrease the running time linearly while incurring a slight loss of accuracy.  
#'                                          
#'                                          Conversely, it is now recommend running IMPUTE2 with large reference panels containing haplotypes of diverse ancestry. 
#'                                          In this context, a rule of thumb suggests setting \code{ref.template.haplotypes} to be smaller than the total size of the reference panel.
#'                                          Imputation accuracy is robust to different values of \code{ref.template.haplotypes} within a sensible range, so it should usually be sufficient to choose a value by intuition. 
#'                                          When in doubt, it is suggested that you err on the side of making \code{ref.template.haplotypes} too large, since often it is that diverse reference panels contain more useful haplotypes than one might expect. 
#' @param seed                              [\code{int}]\cr
#'                                          Initial seed for random number generator.
#'                                          The seed is set using \code{sample(.Machine$integer_max, 1)} unless it is manually overridden with this option.
#' @param no.warnings                       [\code{flag}]\cr
#'                                          Turn warnings off if set to \code{TRUE}, so that the warnings file is not created.
#' @param fill.holes                        [\code{flag}]\cr
#'                                          Turns on the "hole-filling" function if set to \code{TRUE}, which allows SNPs that are typed in the \code{study.haplotypes.file} but not in the reference panel to contribute to the inference.
#' @param no.remove                         [\code{flag}]\cr
#'                                          Prevents the program from discarding SNPs whose alleles cannot be aligned across panels. 
#'                                          Such SNPs will be retained in the output, but they will not be used for inference.
#' @param exec                              [\code{string}]\cr
#'                                          Path to IMPUTE2 executable.
#'                                 
#' @details Simple wrapper to IMPUTE2
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors.
#' @export
#' 
#' @import checkmate imbs
#' 
impute <- function(study.haplotypes.file, genetic.map.file, interval, reference.haplotypes.file, reference.legend.file, output.file, output.snps, gz, output.decimal.places, snp.qc.info, sample.qc.info, predict.genotyped.snps, predict.missing.genotyped.snps, buffer, allow.large.regions, include.buffer.in.output, effective.size, verbose, filter.rules, ref.template.haplotypes, seed, no.warnings, fill.holes, no.remove, exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(study.haplotypes.file, access = "r", add = assertions)
  checkmate::assert_file(genetic.map.file, access = "r", add = assertions)
  checkmate::assert_file(reference.haplotypes.file, access = "r", add = assertions)
  checkmate::assert_file(reference.legend.file, access = "r", add = assertions)
  checkmate::assert_string(output.file, na.ok = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_directory(dirname(output.file), access = "rwx", add = assertions)
  
  checkmate::assert_integer(interval, lower = 0, any.missing = FALSE, all.missing = FALSE, len = 2, null.ok = FALSE, add = assertions)
  
  if (!missing(output.snps)) {
  checkmate::assert_integer(output.snps, lower = 0, upper = 3, any.missing = FALSE, all.missing = FALSE, min.len = 1, max.len = 4, null.ok = FALSE, add = assertions)
    output.snps <- sprintf("-os %s", paste(output.snps, sep = " ", collapse = " "))
  } else {
    output.snps <- ""
  }
  if (!missing(gz)) {
    checkmate::assert_flag(gz, na.ok = FALSE, null.ok = FALSE, add = assertions)
    gz <- "-o_gz"
  } else {
    gz <- ""
  }
  if (!missing(output.decimal.places)) {
    checkmate::assert_int(output.decimal.places, lower = 1, null.ok = FALSE, na.ok = FALSE, add = assertions)
    output.decimal.places <- sprintf("-outdp %d", output.decimal.places)
  } else {
    output.decimal.places <- ""
  }
  if (!missing(snp.qc.info)) {
    checkmate::assert_flag(snp.qc.info, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (snp.qc.info) {
      snp.qc.info <- ""
    } else {
      snp.qc.info <- "-no_snp_qc_info"
    }
  } else {
    snp.qc.info <- ""
  }
  if (!missing(sample.qc.info)) {
    checkmate::assert_flag(sample.qc.info, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (sample.qc.info) {
      sample.qc.info <- ""
    } else {
      sample.qc.info <- "-no_sample_qc_info"
    }
  } else {
    sample.qc.info <- ""
  }
  if (!missing(predict.genotyped.snps)) {
    checkmate::assert_flag(predict.genotyped.snps, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (predict.genotyped.snps) {
      predict.genotyped.snps <- "-pgs"
    } else {
      predict.genotyped.snps <- ""
    }
  } else {
    predict.genotyped.snps <- ""
  }
  if (!missing(predict.missing.genotyped.snps)) {
    checkmate::assert_flag(predict.missing.genotyped.snps, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (predict.missing.genotyped.snps) {
      predict.missing.genotyped.snps <- "-pgs_miss"
    } else {
      predict.missing.genotyped.snps <- ""
    }
  } else {
    predict.missing.genotyped.snps <- ""
  }
  
  if (!missing(buffer)) {
    checkmate::assert_int(buffer, lower = 1, na.ok = FALSE, null.ok = FALSE, add = assertions)
    buffer <- sprintf("-buffer %d", buffer)
  } else {
    buffer <- ""
  }
  if (!missing(allow.large.regions)) {
    checkmate::assert_flag(allow.large.regions, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (allow.large.regions) {
      allow.large.regions <- "-allow_large_regions"
    } else {
      allow.large.regions <- ""
    }
  } else {
    allow.large.regions <- ""
  } 
  if (!missing(include.buffer.in.output)) {
    checkmate::assert_flag(include.buffer.in.output, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (include.buffer.in.output) {
      include.buffer.in.output <- "-include_buffer_in_output"
    } else {
      include.buffer.in.output <- ""
    }
  } else {
    include.buffer.in.output <- ""
  } 
  if (!missing(effective.size)) {
    checkmate::assert_int(effective.size, lower = 1, na.ok = FALSE, null.ok = FALSE, add = assertions)
    effective.size <- sprintf("-Ne %d", effective.size)
  } else {
    effective.size <- ""
  }
  if (!missing(verbose)) {
    checkmate::assert_flag(verbose, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (verbose) {
      verbose <- "-verbose"
    } else {
      verbose <- ""
    }
  } else {
    verbose <- ""
  } 
  
  if (!missing(filter.rules)) {
    checkmate::assert_character(filter.rules, any.missing = FALSE, all.missing = FALSE, min.len = 1, null.ok = FALSE, add = assertions)
    filter.rules <- sprintf("-filt_rules_l %s", paste(filter.rules, sep = " ", collapse = " "))
  } else {
    filter.rules <- ""
  }
  
  if (!missing(ref.template.haplotypes)) {
    checkmate::assert_int(ref.template.haplotypes, lower = 1, na.ok = FALSE, null.ok = FALSE, add = assertions)
    ref.template.haplotypes <- sprintf("-k_hap %d", ref.template.haplotypes)
  } else {
    ref.template.haplotypes <- ""
  }
   
  if (!missing(seed)) {
    checkmate::assert_int(seed, lower = 1, na.ok = FALSE, null.ok = FALSE, add = assertions)
    seed <- sprintf("-seed %d", seed)
  } else {
    seed <- sprintf("-seed %d", sample(.Machine$integer.max, 1))
  }
  if (!missing(no.warnings)) {
    checkmate::assert_flag(no.warnings, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (no.warnings) {
      no.warnings <- "-no_warn"
    } else {
      no.warnings <- ""
    }
  } else {
    no.warnings <- ""
  }
  if (!missing(fill.holes)) {
    checkmate::assert_flag(fill.holes, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (fill.holes) {
      fill.holes <- "-fill_holes"
    } else {
      fill.holes <- ""
    }
  } else {
    fill.holes <- ""
  }
  if (!missing(no.remove)) {
    checkmate::assert_flag(no.remove, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (no.remove) {
      no.remove <- "-no_remove"
    } else {
      no.remove <- ""
    }
  } else {
    no.remove <- ""
  }
  
  imbs::assert_command(exec, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c(
      "-known_haps_g", study.haplotypes.file, 
      "-m", genetic.map.file, 
      "-int", interval, 
      "-h", reference.haplotypes.file, 
      "-l", reference.legend.file,
      "-o", output.file, 
      output.snps, 
      gz, 
      output.decimal.places, 
      snp.qc.info, 
      sample.qc.info,
      predict.genotyped.snps, 
      predict.missing.genotyped.snps, 
      buffer, 
      allow.large.regions,
      include.buffer.in.output,
      effective.size,
      verbose, 
      filter.rules,
      ref.template.haplotypes, 
      seed,
      no.warnings,
      fill.holes,
      no.remove
    )
  )
  
}