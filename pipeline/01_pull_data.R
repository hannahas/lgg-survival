# Bioconductor package that talks to NCI Genomic Data Commons (GDC) API to download and prepare data for analysis
library(TCGAbiolinks)
# for path handling
library(here)

# Creates raw folder inside the data folder
raw_dir <- here("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

# builds the query, downloads, and parses into a data frame
clinical <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
# serializes the data frame into a RDS file for later use (preferred to write.csv)
saveRDS(clinical, file.path(raw_dir, "clinical_raw.rds"))

subtypes <- PanCancerAtlas_subtypes()
lgg_subtypes <- subtypes[subtypes$cancer.type == "LGG", ]
saveRDS(lgg_subtypes, file.path(raw_dir, "molecular_subtypes_raw.rds"))

# Your revised hypothesis structure, now that you know what you have:
# Does methylation-based molecular subtype stratify survival better than histologic grade (G2/G3)? — the reproduction
# Does an expression-based model add prognostic value beyond molecular subtype? — your contribution
# Does it hold under nested CV, given ~125 events? — the validation

# Searches GDC file index and returns a manifest, list of files, ids, sizes, barcodes
expr_query <- GDCquery(
  project       = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Slow step, writes each file to data/raw/gcd_downloads/
GDCdownload(expr_query, directory = file.path(raw_dir, "gdc_downloads"))

# Puts all these individual files into a single SummarizedExperiment object, 
# which is a Bioconductor data structure that holds the expression matrix, sample metadata, and feature metadata
expr_se <- GDCprepare(expr_query, directory = file.path(raw_dir, "gdc_downloads"))

# saves into one object
# 60,660 genes × 534 samples. (516 patients)
saveRDS(expr_se, file.path(raw_dir, "expression_raw.rds"))

library(SummarizedExperiment)

pt <- substr(colnames(expr_se), 1, 12)
idx <- match(pt, clinical$submitter_id)
table(expr_se$paper_Grade, clinical$tumor_grade[idx], useNA = "ifany")

missing_idh <- is.na(expr_se$paper_IDH.codel.subtype)
table(missing_idh, clinical$vital_status[idx])
table(missing_idh, clinical$tumor_grade[idx])

# 21 patients lacked molecular subtype calls and were excluded. 
# This missingness was not random — excluded patients had a higher mortality rate
# (62% vs. 24%) despite skewing toward lower histologic grade. 
# This pattern is consistent with these tumors being disproportionately IDH-wildtype, 
# the subgroup with the worst prognosis. Their exclusion likely makes the 
# analysis cohort modestly more favorable than the full LGG population, 
# and may attenuate the apparent survival difference between molecular subtypes.