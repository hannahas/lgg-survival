library(SummarizedExperiment)
library(survival)
library(here)

raw_dir  <- here("data", "raw")
proc_dir <- here("data", "processed")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)

clinical <- readRDS(file.path(raw_dir, "clinical_raw.rds"))
expr_se  <- readRDS(file.path(raw_dir, "expression_raw.rds"))

# Build the patient-level table
# Primary tumors only — drops the 18 recurrences
expr_se <- expr_se[, expr_se$sample_type == "Primary Tumor"]

# Pull molecular calls out of the SummarizedExperiment, keyed by patient
molecular <- data.frame(
  patient_id     = substr(colnames(expr_se), 1, 12),
  sample_barcode = colnames(expr_se),
  idh_codel      = as.character(expr_se$paper_IDH.codel.subtype),
  idh_status     = as.character(expr_se$paper_IDH.status),
  codel_1p19q    = as.character(expr_se$paper_X1p.19q.codeletion),
  meth_subtype   = as.character(expr_se$paper_Supervised.DNA.Methylation.Cluster),
  stringsAsFactors = FALSE
)