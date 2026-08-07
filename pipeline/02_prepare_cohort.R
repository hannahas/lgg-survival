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

cohort <- merge(
  molecular,
  clinical[, c("submitter_id", "vital_status", "days_to_death",
               "days_to_last_follow_up", "tumor_grade",
               "primary_diagnosis", "age_at_index")],
  by.x = "patient_id", by.y = "submitter_id",
  all.x = TRUE
)

# Survival time: death date if dead, last follow-up if alive
cohort$os_days <- ifelse(
  cohort$vital_status == "Dead",
  cohort$days_to_death,
  cohort$days_to_last_follow_up
)

cohort$os_event <- ifelse(cohort$vital_status == "Dead", 1, 0)

ledger <- data.frame(step = "Primary tumor samples", n = nrow(cohort))

cohort <- cohort[!is.na(cohort$idh_codel), ]
ledger <- rbind(ledger, data.frame(step = "Has IDH/codel subtype", n = nrow(cohort)))

cohort <- cohort[cohort$vital_status %in% c("Alive", "Dead"), ]
ledger <- rbind(ledger, data.frame(step = "Vital status known", n = nrow(cohort)))

cohort <- cohort[!is.na(cohort$os_days) & cohort$os_days > 0, ]
ledger <- rbind(ledger, data.frame(step = "Valid follow-up time", n = nrow(cohort)))

print(ledger)
saveRDS(ledger, file.path(proc_dir, "exclusion_ledger.rds"))

cohort$idh_codel <- factor(
  cohort$idh_codel,
  levels = c("IDHwt", "IDHmut-non-codel", "IDHmut-codel")
)

cohort$tumor_grade <- factor(cohort$tumor_grade, levels = c("G2", "G3"))

cohort$meth_subtype <- factor(cohort$meth_subtype)

# Split data
set.seed(42)

strata <- interaction(cohort$idh_codel, cohort$os_event)

train_idx <- unlist(lapply(split(seq_len(nrow(cohort)), strata), function(ix) {
  sample(ix, size = floor(0.8 * length(ix)))
}))

cohort$split <- ifelse(seq_len(nrow(cohort)) %in% train_idx, "train", "test")

saveRDS(cohort, file.path(proc_dir, "cohort.rds"))