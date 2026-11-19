# Cross-check: our Python PhenoAge vs. the BioAge R package (orig = TRUE uses the Levine 2018 coefficients).
# Run after phenoage.py:  Rscript r/check_bioage.R
# Install once:  remotes::install_github("dayoonkwon/BioAge")
suppressMessages(library(BioAge))
d <- read.csv("data/phenoage.csv")
names(d)[names(d) == "albumin_gl"] <- "albumin_gL"
d$time <- d$permth_exm
d$status <- d$mortstat
bm <- c("albumin_gL", "creat_umol", "glucose_mmol", "lncrp", "lymph", "mcv", "rdw", "alp", "wbc")
x <- phenoage_calc(data = d, biomarkers = bm, orig = TRUE)$data
ok <- is.finite(x$phenoage0)  # BioAge's formula returns Inf when mortality score rounds to 1
cat("rows:", nrow(x), "\n")
cat("rows where BioAge is not finite:", sum(!ok), "\n")
cat("max abs diff (years):", max(abs(x$phenoage0[ok] - d$phenoage[ok])), "\n")
cat("correlation:", cor(x$phenoage0[ok], d$phenoage[ok]), "\n")
if (any(!ok)) cat("non-finite rows: Python PhenoAge range =", range(d$phenoage[!ok]), "\n")
