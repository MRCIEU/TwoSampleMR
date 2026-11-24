# Records the opengwasr output that tests/testthat mocks. Run by hand against the development Data
# API (OPENGWAS_URL, OPENGWAS_DEV_KEY set); the results are committed, tests never call this.
library(opengwasr)
src <- og_api()
dir <- "tests/testthat/fixtures"

analyses <- og_analyses(src, c("ukb-b-19953", "ebi-a-GCST003116"))
saveRDS(analyses, file.path(dir, "analyses.rds"))

# BMI instruments, clumped on the API. Of the strongest 100 index variants we keep the 12 strongest
# and five more: rows 18, 24, 35 and 60 have no rsID (their SNP is the ALID), and 79 is an Imputed
# Association in coronary artery disease, as are 18, 24 and 60.
exposure <- suppressWarnings(og_tophits(src, "ukb-b-19953", limit = 100))
exposure <- exposure[c(1:12, 18, 24, 35, 60, 79), ]
rownames(exposure) <- NULL
attr(exposure, "og_provenance") <- NULL # describes the 100 rows, not these 17
saveRDS(exposure, file.path(dir, "tophits-bmi.rds"))

# The same variants in coronary artery disease.
outcome <- og_associations(src, ifelse(is.na(exposure$rsid), exposure$alid, exposure$rsid),
                           "ebi-a-GCST003116")
saveRDS(outcome, file.path(dir, "associations-cad.rds"))

# Coronary artery disease instruments, a binary exposure (the 10 strongest index variants).
tophits_cad <- suppressWarnings(og_tophits(src, "ebi-a-GCST003116", limit = 10))
attr(tophits_cad, "og_provenance") <- NULL
saveRDS(tophits_cad, file.path(dir, "tophits-cad.rds"))

# Note: on 24 Nov 2026 the dev release of coronary artery disease no longer returned Imputed
# Associations, so associations-cad.rds is the recording of 23 Nov, which has four.
