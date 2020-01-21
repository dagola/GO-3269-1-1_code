
source("../init.R", chdir = TRUE)

##============================================================================+
## Packages ----
##============================================================================+
pacman::p_load(data.table)
pacman::p_install_gh("imbs-hl/imbs")
pacman::p_load(batchtools)
pacman::p_load(DBI)
pacman::p_load(RSQLite)
pacman::p_load(cowplot)
pacman::p_load(parallelMap)
pacman::p_load(gtools)
pacman::p_load(PRROC)
pacman::p_load(pROC)

if (system %in% c("damian.novalocal", "n01", "n04")) {
  ramdisk <- "/mnt/ramdisk"
}

##============================================================================+
## Software ----
##============================================================================+
bcftools_exec <- switch(system,
                        "damian.novalocal" = "/home/damian/bin/bcftools1.9",
                        "/software/bin/bcftools1.9"
)
genotype_harmonizer_exec <- switch(system,
                                   "damian.novalocal" = "/home/damian/bin/GenotypeHarmonizer1.4.20",
                                   "/software/bin/GenotypeHarmonizer1.4.20"
) 
plink_exec <- switch(system,
                     "damian.novalocal" = "/home/damian/bin/plink1.90b6.10",
                     "/software/bin/plink1.9b6.9"
) 
plink2_exec <- switch(system,
                      "damian.novalocal" = "/home/damian/bin/plink2.0a20190628",
                      "/software/bin/plink2a20190708"
) 
gcta_exec <- switch(system,
                    "damian.novalocal" = "/home/damian/bin/gcta1.92.2b",
                    "/software/bin/gcta1.92.2b"
)
prsice_exec <- switch(system,
                      "damian.novalocal" = "/home/damian/bin/prsice2.2.2",
                      "/software/bin/prsice2.2.2"
)
shapeit2_exec <- "/software/bin/shapeit.v2.r904.GLIBCv2.17"
impute2_exec <- "/software/bin/impute2.3.2"
