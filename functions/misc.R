
#' Get chromosome, base pair position and alleles from BGEN index file
#'
#' @param bgi.file  [\code{string}]\cr
#'                  File path to BGEN index file.
#' @param keys      [\code{character}]\cr
#'                  Character vector of column names to become indices of returned \code{data.table}.
#'
#' @return \code{data.table} with columns CHR, POS, REF and ALT.
#' @export
#'
#' @import dBI RSQLite data.table checkmate
#' 
query_bgi_chr_pos <- function(bgi.file, keys) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(bgi.file, access = "r", add = assertions)
  checkmate::assert_character(keys, min.len = 0, max.len = 4, unique = TRUE, add = assertions)
  checkmate::assert_subset(keys, choices = c("CHR", "POS", "REF", "ALT"), add = TRUE)
  
  checkmate::reportAssertions(assertions)
  
  index <- DBI::dbConnect(RSQLite::SQLite(), bgi.file)
  on.exit(DBI::dbDisconnect(index), add = TRUE)
  chr_pos <- DBI::dbGetQuery(
    conn = index, 
    statement = "SELECT CAST(chromosome AS INT) AS CHR, position AS POS, allele1 AS REF, allele2 AS ALT FROM Variant"
  )
  data.table::setDT(chr_pos, key = keys)
  return(chr_pos)
}

#' Get chromosome, base pair position and alleles from VCF file
#'
#' @param vcf.file        [\code{string}]\cr
#'                        File path to VCF file.
#' @param col.names       [\code{character}]\cr
#'                        Column names to be assigned to \code{data.table} returned.
#' @param keys            [\code{character}]\cr
#'                        Character vector of column names to become indices of returned \code{data.table}.
#' @param bcftools.exec   [\code{string}]\cr
#'                        File path to bcftools executable.
#'
#' @return  \code{data.table} with column names given by \code{col.names}.
#' @export
#'
#' @import data.table checkmate imbs
#' 
function(vcf.file, col.names, keys, bcftools.exec) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, access = "r", add = assertions)
  checkmate::assert_character(col.names, min.len = 0, max.len = 4, unique = TRUE, add = assertions)
  checkmate::assert_character(keys, min.len = 0, max.len = 4, unique = TRUE, add = assertions)
  checkmate::assert_subset(keys, choices = col.names, add = TRUE)
  imbs::assert_command(bcftools.exec)
  
  checkmate::reportAssertions(assertions)
  
  chr_pos <- data.table::fread(
    cmd = sprintf("%s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s", bcftools.exec, vcf.file), 
    col.names = col.names,
    key = keys
  )
  return(chr_pos)
}

get_EB_phenotype_data <- function() {
  
  env <- new.env()
  nm <- load(file = eb_phenotypes_file, envir = env)[1]
  
  return(data.table::as.data.table(env[[nm]]))
  
}

get_UKB_phenotype_data <- function() {
  
  select_fields <- c(
    "ID" = "f.eid", 
    "SEX" = "f.31.0.0", 
    "BIRTH_YEAR" = "f.34.0.0", 
    "BIRTH_MONTH" = "f.52.0.0", 
    "ETHNICITY_0" = "f.21000.0.0",
    "ETHNICITY_1" = "f.21000.1.0",
    "ETHNICITY_2" = "f.21000.2.0",
    "ICD10_" = sprintf("f.41270.0.%d", 0:212), 
    "MAIN_ICD10_" = sprintf("f.41202.0.%d", 0:65),  
    "SECONDARY_ICD10_" = sprintf("f.41204.0.%d", 0:183), 
    "ICD9_" = sprintf("f.41271.0.%d", 0:46), 
    "MAIN_ICD9_" = sprintf("f.41203.0.%d", 0:27), 
    "SECONDARY_ICD9_" = sprintf("f.41205.0.%d", 0:29), 
    "OPCS4_" = sprintf("f.41272.0.%d", 0:116), 
    "MAIN_OPCS4_" = sprintf("f.41200.0.%d", 0:48), 
    "SECONDARY_OPCS4_" = sprintf("f.41210.0.%d", 0:85),
    "DEATH_PRIMARY_ICD10_0" = "f.40001.0.0",
    "DEATH_PRIMARY_ICD10_1" = "f.40001.1.0",
    "DEATH_SECONDARY_ICD10_0_" = sprintf("f.40002.0.%d", 1:13),
    "DEATH_SECONDARY_ICD10_1_" = sprintf("f.40002.1.%d", 1:13)
  )
  
  ukb_phenotype_data <- data.table::fread(
    file = ukb_phenotypes_file, 
    select = select_fields,
    col.names = names(select_fields)
  )
  
  ukb_phenotype_data[, ID := as.character(ID)]
  ukb_phenotype_data[, SEX := ordered(SEX, levels = c(0, 1), labels = c("Female", "Male"))]
  ukb_phenotype_data[, BIRTH_MONTH := ordered(BIRTH_MONTH, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))]
  ukb_phenotype_data[, ETHNICITY_0 := ordered(ETHNICITY_0, levels = c(-3, -1, 1, 2, 3, 4, 5, 6, 1001, 1002, 1003, 2001, 2002, 2003, 2004, 3001, 3002, 3003, 3004, 4001, 4002, 4003), labels = c("Prefer not to answer", "Do not know", "White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group", "British", "Irish", "Any other white background", "White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background", "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Caribbean", "African", "Any other Black background"))]
  ukb_phenotype_data[, ETHNICITY_1 := ordered(ETHNICITY_1, levels = c(-3, -1, 1, 2, 3, 4, 5, 6, 1001, 1002, 1003, 2001, 2002, 2003, 2004, 3001, 3002, 3003, 3004, 4001, 4002, 4003), labels = c("Prefer not to answer", "Do not know", "White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group", "British", "Irish", "Any other white background", "White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background", "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Caribbean", "African", "Any other Black background"))]
  ukb_phenotype_data[, ETHNICITY_2 := ordered(ETHNICITY_2, levels = c(-3, -1, 1, 2, 3, 4, 5, 6, 1001, 1002, 1003, 2001, 2002, 2003, 2004, 3001, 3002, 3003, 3004, 4001, 4002, 4003), labels = c("Prefer not to answer", "Do not know", "White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group", "British", "Irish", "Any other white background", "White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background", "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Caribbean", "African", "Any other Black background"))]

  return(ukb_phenotype_data)
  
}
