
source("init.R")

# Create common SNP sets ----

## Load variant positions
col_names <- c("CHR", "POS", "REF", "ALT")
key_names <- c("CHR", "POS")

### UKB in BGEN files
reg_load_bgen_variants <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "load_bgen_variants"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("DBI", "RSQLite", "data.table")
)

ids <- batchtools::batchMap(
  fun = query_bgi_chr_pos,
  bgi.file = c(sprintf("%s_chr%d_v3.bgen.bgi", ukb_file_prefix, autosomes)),
  more.args = list(
    keys = key_names
  ),
  reg = reg_load_bgen_variants
)
ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_load_bgen_variants, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 9000,
    partition = "batch", walltime = 20,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Load_Variants", project_name)
  )
)

batchtools::waitForJobs(reg = reg_load_bgen_variants)

ukb_chr_pos <- data.table::rbindlist(batchtools::reduceResultsList(reg = reg_load_bgen_variants))
data.table::setkeyv(ukb_chr_pos, key_names)

### EB in VCF files
reg_load_vcf_variants <- imbs::load_or_create_registry(
  file.dir = file.path(registries_dir, "load_vcf_variants"), 
  work.dir = proc_dir,
  writeable = TRUE, 
  overwrite = FALSE,
  packages = c("data.table")
)

ids <- batchtools::batchMap(
  fun = query_vcf_chr_pos,
  vcf.file = c(sprintf("%s_chr%d.vcf.gz", eb_file_prefix, autosomes)),
  more.args = list(
    keys = key_names,
    col.names = col_names,
    bcftools.exec = bcftools_exec
  ),
  reg = reg_load_vcf_variants
)
ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids, 
  reg = reg_load_vcf_variants, 
  resources = list(
    ntasks = 1, ncpus = 1, memory = 4000,
    partition = "batch", walltime = 0,
    chunks.as.arrayjobs = TRUE,
    name = sprintf("%s_Load_Variants", project_name)
  )
)

batchtools::waitForJobs(reg = reg_load_vcf_variants)

eb_chr_pos <- data.table::rbindlist(batchtools::reduceResultsList(reg = reg_load_vcf_variants), fill = TRUE)
data.table::setkeyv(eb_chr_pos, key_names)

### DE in PLINK1 bim file
de_chr_pos <- data.table::fread(
  file = sprintf("%s.bim", de_file_prefix),
  select = c(1, 4, 5, 6),
  col.names = col_names,
  key = key_names
)

### Nikpay summary stats in plain text
nikpay_chr_pos <- data.table::fread(
  file = nikpay_summary_statistics_file,
  select = c(2, 3, 5, 4),
  col.names = col_names, 
  key = key_names
)

### 1kG in PLINK1 bim file
g1k_chr_pos <- data.table::fread(
  file = sprintf("%s.bim", g1k_file_prefix),
  select = c(1, 4, 5, 6),
  col.names = col_names,
  key = key_names
)

## Define allowed alles
flip_vector <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
single_nucleotides <- 'which(REF %in% c("A", "C", "G", "T") & ALT %in% c("A", "C", "G", "T"))'

## Find common SNPs between UKB and EB
ukb_eb_snps <- ukb_chr_pos[ukb_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
  eval(parse(text = single_nucleotides)) # select SNPs only
  ][
    REF != flip_vector[ALT] # non-palindromic sites only
    ][
      eb_chr_pos[eb_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
        eval(parse(text = single_nucleotides)) # select SNPs only
        ][
          REF != flip_vector[ALT] # non-palindromic sites only, 
          ],
      nomatch = 0L
      ][
        (REF == i.REF & ALT == i.ALT) | 
          (REF == i.ALT & ALT == i.REF) | 
          (REF == flip_vector[i.REF] & ALT == flip_vector[i.ALT]) | 
          (REF == flip_vector[i.ALT] & ALT == flip_vector[i.REF]) # Select orientable sites only
        ]
data.table::fwrite(
  x = unique(ukb_eb_snps[, .(CHR, POS)]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_snps_pos_file
)
data.table::fwrite(
  x = unique(ukb_eb_snps[, .(ID = sprintf("%s:%s", CHR, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_snps_file
)
data.table::fwrite(
  x = unique(ukb_eb_snps[, .(CHR, POS, POS, "UKB_EB")]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_plink_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_snps[, .(RANGE = sprintf("%02d:%d-%d", CHR, POS, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_bgenix_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[(REF == flip_vector[i.REF] & ALT == flip_vector[i.ALT]) | 
                                         (REF == flip_vector[i.ALT] & ALT == flip_vector[i.REF]),
                                       .(ID = sprintf("%s:%s", CHR, POS))]),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  file = eb_flip_file
)

## Find common SNPs between UKB, EB and DE
ukb_eb_de_snps <- ukb_eb_snps[
  de_chr_pos[de_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
    eval(parse(text = single_nucleotides)) # select SNPs only
    ][
      REF != flip_vector[ALT] # non-palindromic sites only, 
      ],
  nomatch = 0L
  ][
    (REF == i.REF.1 & ALT == i.ALT.1) | 
      (REF == i.ALT.1 & ALT == i.REF.1) | 
      (REF == flip_vector[i.REF.1] & ALT == flip_vector[i.ALT.1]) | 
      (REF == flip_vector[i.ALT.1] & ALT == flip_vector[i.REF.1]) # Select orientable sites only
    ]
data.table::fwrite(
  x = unique(ukb_eb_de_snps[, .(CHR, POS)]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_snps_pos_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_snps[, .(ID = sprintf("%s:%s", CHR, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_snps_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_snps[, .(CHR, POS, POS, "UKB_EB_Nikpay")]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_plink_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_snps[, .(RANGE = sprintf("%02d:%d-%d", CHR, POS, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_bgenix_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_snps[(REF == flip_vector[i.REF.1] & ALT == flip_vector[i.ALT.1]) | 
                              (REF == flip_vector[i.ALT.1] & ALT == flip_vector[i.REF.1]), 
                            .(ID = sprintf("%s:%s", CHR, POS))]),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  file = de_flip_file
)

## Find common SNPs between UKB, EB and Nikpay summary statistics
ukb_eb_de_nikpay_snps <- ukb_eb_de_snps[
  nikpay_chr_pos[nikpay_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
    eval(parse(text = single_nucleotides)) # select SNPs only
    ][
      REF != flip_vector[ALT] # non-palindromic sites only, 
      ],
  nomatch = 0L
  ][
    (REF == i.REF.2 & ALT == i.ALT.2) | 
      (REF == i.ALT.2 & ALT == i.REF.2) | 
      (REF == flip_vector[i.REF.2] & ALT == flip_vector[i.ALT.2]) | 
      (REF == flip_vector[i.ALT.2] & ALT == flip_vector[i.REF.2]) # Select orientable sites only
    ]
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_snps[, .(CHR, POS)]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_snps_pos_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_snps[, .(ID = sprintf("%s:%s", CHR, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_snps_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_snps[, .(CHR, POS, POS, "UKB_EB_Nikpay")]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_plink_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_snps[, .(RANGE = sprintf("%02d:%d-%d", CHR, POS, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_bgenix_range_file
)

## Find common SNPs between UKB, EB, Nikpay summary statistics and 1000G
ukb_eb_de_nikpay_g1k_snps <- ukb_eb_de_nikpay_snps[
  g1k_chr_pos[g1k_chr_pos[, .N, by = .(CHR, POS)][N==1]][ # select biallelic sites only
    eval(parse(text = single_nucleotides)) # select SNPs only
    ][
      REF != flip_vector[ALT] # non-palindromic sites only
      ],
  nomatch = 0L
  ][
    (REF == i.REF.3 & ALT == i.ALT.3) | 
      (REF == i.ALT.3 & ALT == i.REF.3) | 
      (REF == flip_vector[i.REF.3] & ALT == flip_vector[i.ALT.3]) | 
      (REF == flip_vector[i.ALT.3] & ALT == flip_vector[i.REF.3])
    ]
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[, .(CHR, POS)]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_g1k_snps_pos_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[, .(ID = sprintf("%s:%s", CHR, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_g1k_snps_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[, .(CHR, POS, POS, "UKB_EB_Nikpay_1kG")]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_g1k_plink_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[, .(RANGE = sprintf("%02d:%d-%d", CHR, POS, POS))]), 
  quote = FALSE, 
  sep = "\t",
  col.names = FALSE, 
  file = ukb_eb_de_nikpay_g1k_bgenix_range_file
)
data.table::fwrite(
  x = unique(ukb_eb_de_nikpay_g1k_snps[(REF == flip_vector[i.REF.3] & ALT == flip_vector[i.ALT.3]) | 
                                         (REF == flip_vector[i.ALT.3] & ALT == flip_vector[i.REF.3]),
                                       .(ID = sprintf("%s:%s", CHR, POS))]),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  file = g1k_flip_file
)
