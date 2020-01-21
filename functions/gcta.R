
#' Create a Genetic Relationship Matrix using GCTA
#'
#' @param bfile         [\code{string}]\cr
#'                      The basename of the binary PLINK files.
#' @param output.prefix [\code{string}]\cr
#'                      The basename of the output files.
#' @param chr           [\code{int}]\cr
#'                      Restrict calculation to the given chromosome.
#' @param maf           [\code{maf}]\cr
#'                      Filter variants with minor allele frequency lower than this threshold from the analysis.
#' @param keep          [\code{string}]\cr
#'                      Path to file containing FID and IID of samples to include in the analysis.
#' @param remove        [\code{string}]\cr
#'                      Path to file containing FID and IID of samples to exclude from the analysis.
#' @param extract       [\code{string}]\cr
#'                      Path to file containing variant IDs to include in the analysis.
#' @param exclude       [\code{string}]\cr
#'                      Path to file containing variant IDs to exclude from the analysis.
#' @param num.threads   [\code{int}]\cr
#'                      Number of CPUs usable by GCTA. Default is determined by SLURM environment variables and at least 1.
#' @param exec          [\code{string}]\cr
#'                      Path to GCTA executable.
#' @param ...           Additional arguments directly passed to the system call arguments.
#'
#' @return The system call output as character vector.
#' 
#' @export
gcta_make_grm <- function(bfile, output.prefix, chr, maf, keep, remove, extract, exclude, num.threads, exec, ...) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
  input <- sprintf("--bfile %s", bfile)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (!missing(chr)) {
    checkmate::assert_integer(chr, lower = 1, upper = 25, any.missing = FALSE, min.len = 1, max.len = 25, unique = TRUE, all.missing = FALSE, add = assertions)
    chr <- sprintf("--chr %d", chr)
  } else {
    chr <- ""
  }
  
  if (!missing(maf)) {
    checkmate::assert_number(maf, lower = 0, upper = 1, finite = TRUE, na.ok = FALSe, null.ok = FALSE, add = assertions)
    maf <- sprintf("--maf %f", maf)
  } else {
    maf <- ""
  }
  
  if (missing(remove)) {
    remove <- ""
  } else {
    checkmate::assert_file(remove, add = assertions)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(keep)) {
    keep <- ""
  } else {
    checkmate::assert_file(keep, add = assertions)
    keep <- sprintf("--keep %s", keep)
  }
  
  if (missing(exclude)) {
    exclude <- ""
  } else {
    checkmate::assert_file(exclude, add = assertions)
    exclude <- sprintf("--exclude %s", exclude)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    checkmate::assert_file(extract, add = assertions)
    extract <- sprintf("--extract %s", extract)
  }
  
  imbs::assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             keep,
             remove,
             extract,
             exclude,
             chr,
             maf,
             "--make-grm",
             "--out", output.prefix, ...)
  )
}

#' Merge multiple Genetic Relationship Matrices with GCTA
#'
#' @param grms.list     [\code{string}]\cr
#'                      Path to file with prefix paths to partial GRMs.
#' @param output.prefix [\code{sting}]\cr
#'                      The basename of the output files.
#' @param exec          [\code{string}]\cr
#'                      Path to GCTA executable.
#'
#' @return The system call output as character vector.
#' 
#' @export
gcta_merge_grms <- function(grms.list, output.prefix, exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(grms.list, add = assertions)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  imbs::assert_command(exec, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c("--mgrm", grms.list,
             "--make-grm",
             "--out", output.prefix)
  )
}

#' Perfom a Principal Component Analysis with GCTA
#'
#' @param bfile         [\code{string}]\cr
#'                      The basename of the binary PLINK files.
#' @param grm.prefix    [\code{string}]\cr
#'                      The basename of the corresponding precomputed GRM.
#' @param num.pcs       [\code{int}]\cr
#'                      Number of principal components to compute.
#' @param output.prefix [\code{string}]\cr
#'                      The basename of output files.
#' @param ...           Additional arguments directly passed to the system call arguments.
#' @param num.threads   [\code{int}]\cr
#'                      Number of CPUs usable by GCTA. Default is determined by SLURM environment variables and at least 1.
#' @param exec          [\code{string}]\cr
#'                      Path to GCTA executable.
#' 
#' @details First a PCA is done using the GRM. Then loadings are produced using the PCA result and the PLINK binary file.
#'
#' @return The system call output as character vector.
#' 
gcta_pca <- function(bfile, grm.prefix, num.pcs, output.prefix, ..., num.threads, exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(sprintf("%s.grm.bin", grm.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.grm.id", grm.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.grm.N.bin", grm.prefix), add = assertions)
  checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
  checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
  checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
  
  checkmate::assertInt(num.pcs, lower = 1, null.ok = FALSE, na.ok = FALSE)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  imbs::assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c("--grm", grm.prefix,
             "--pca", num.pcs,
             "--threads", num.threads,
             "--out", output.prefix, ...)
  )
  
  imbs::system_call(
    bin = exec,
    args = c("--bfile", bfile,
             "--pc-loading", output.prefix,
             "--threads", num.threads,
             "--out", output.prefix, ...)
  )
  
}

#' Projection of Samples in a Pre-computed PCA Space using GCTA
#'
#' @param bfile               [\code{string}]\cr
#'                            The basename of the binary PLINK files of the samples to be projected.
#' @param pca.loadings.prefix [\code{string}]\cr
#'                            The basename of loadings from a GCTA principal component analysis (see \code{\link{gcta_pca}}).
#' @param num.pcs             [\code{int}]\cr
#'                            The number of eigenvectors to procude.
#' @param output.prefix       [\code{string}]\cr
#'                            The basename of output files.
#' @param ...                 Additional arguments directly passed to the system call arguments.
#' @param num.threads         [\code{int}]\cr
#'                            Number of CPUs usable by GCTA. Default is determined by SLURM environment variables and at least 1.
#' @param exec                [\code{string}]\cr
#'                            Path to GCTA executable.
#'                            
#' @details You should ensure that the PLINK file contains only SNPs used in the PCA before!
#'
#' @return The system call output as character vector.
#' 
#' @export
#' 
gcta_pca_projection <- function(bfile, pca.loadings.prefix, num.pcs, output.prefix, ..., num.threads, exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
  checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
  checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.pcl", pca.loadings.prefix), add = assertions)
  
  checkmate::assertInt(num.pcs, lower = 1, null.ok = FALSE, na.ok = FALSE)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  imbs::assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  imbs::system_call(
    bin = exec,
    args = c("--bfile", bfile,
             "--project-loading", pca.loadings.prefix, num.pcs,
             "--threads", num.threads,
             "--out", output.prefix, ...)
  )
  
}