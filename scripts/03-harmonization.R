
source("init.R")

## Define allowed alles
flip_vector <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
single_nucleotides <- 'which(REF %in% c("A", "C", "G", "T") & ALT %in% c("A", "C", "G", "T"))'
key_names <- c("CHR", "POS")

### UKB chromosomal positions
reg_load_bgen_variants <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "load_bgen_variants"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("DBI", "RSQLite", "data.table")
)

batchtools::waitForJobs(reg = reg_load_bgen_variants)

ukb_chr_pos <- data.table::rbindlist(batchtools::reduceResultsList(reg = reg_load_bgen_variants))
data.table::setkeyv(ukb_chr_pos, key_names)

# Harmonize summary statistics alleles to UKB ----
nikpay_summary_statistics <- data.table::fread(
  file = nikpay_summary_statistics_file
)

harmonized_summary_statistics <- ukb_chr_pos[ukb_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
  eval(parse(text = single_nucleotides)) # select SNPs only
  ][
    REF != flip_vector[ALT] # non-palindromic sites only
    ][
    nikpay_summary_statistics[
      which(effect_allele %in% flip_vector & other_allele %in% flip_vector) # select SNPs only
      ], on = c("CHR" = "chromosome", "POS" = "base_pair_location"),
    nomatch = 0L
    ][
      (REF == effect_allele & ALT == other_allele) | 
        (REF == other_allele & ALT == effect_allele) | 
        (REF == flip_vector[effect_allele] & ALT == flip_vector[other_allele]) | 
        (REF == flip_vector[other_allele] & ALT == flip_vector[effect_allele]) # Select orientable sites only
      ]

data.table::fwrite(
  x = harmonized_summary_statistics[, .(
    "#CHROM" = sprintf("%s", CHR), 
    POS, 
    ID = sprintf("%s:%s", CHR, POS), 
    A1 = effect_allele, 
    AX = other_allele, 
    BETA = beta, 
    P = pvalue, 
    A1_FREQ = effect_allele_frequency, 
    N = nikpay_sample_size
  )],
  file = nikpay_summary_statistics_harmonized_file,
  sep = " ", quote = FALSE, na = "NA", col.names = TRUE
)
