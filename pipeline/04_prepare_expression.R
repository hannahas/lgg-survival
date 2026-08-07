library(SummarizedExperiment)
library(here)
library(DESeq2)

# ---- Load cached expression data ----
# Read the raw expression object built in 01_pull_data.R rather than
# re-downloading from GDC.
raw_dir <- here("data", "raw")
expr_se <- readRDS(file.path(raw_dir, "expression_raw.rds"))

# ---- Restrict to primary tumors ----
# The raw pull included 18 recurrent-tumor samples alongside 516 primary
# tumors. Recurrences aren't part of this cohort definition, so drop them
# here to match the same sample_type filter applied to the clinical cohort.
expr_se <- expr_se[, expr_se$sample_type == "Primary Tumor"]

# ---- Restrict to protein-coding genes ----
# Of ~60,660 annotated features, only ~20,000 are protein-coding; the rest
# is lncRNA, pseudogenes, and other non-coding annotation that adds noise
# and multiple-testing burden without adding signal for this analysis.
# This filter uses only static gene annotation (rowData), not expression
# values or outcome data, so it's safe to apply on the full dataset before
# any train/test split logic.
coding_genes <- rowData(expr_se)$gene_type == "protein_coding"
expr_se <- expr_se[coding_genes, ]
dim(expr_se)

# ---- Filter to expressed genes ----
# Even among protein-coding genes, many are silent or near-silent in
# glioma tissue specifically. Keep only genes with at least 10 raw counts
# in at least 10% of samples — a standard low-expression filter that
# removes genes contributing mostly noise. Like the filter above, this is
# outcome-blind: it only asks whether a gene is expressed at all, not
# whether it relates to survival, so it's safe to compute on all samples
# rather than training data only.
counts <- assay(expr_se, "unstranded")
# Keep genes with at least 10 counts in at least 10% of samples
n_samples <- ncol(expr_se)
keep <- rowSums(counts >= 10) >= (0.10 * n_samples)
expr_se <- expr_se[keep, ]

# ---- Restrict to the final analysis cohort ----
# expr_se still contains all 516 primary-tumor samples, but cohort.rds
# reflects the 508 patients remaining after clinical exclusions (missing
# molecular subtype, unknown vital status, invalid follow-up time — see
# 02_prepare_cohort.R). Align expression data to that same cohort so every
# downstream sample has a matching row in cohort.
cohort <- readRDS(here("data", "processed", "cohort.rds"))
patient_ids <- substr(colnames(expr_se), 1, 12)
keep_samples <- patient_ids %in% cohort$patient_id
expr_se <- expr_se[, keep_samples]
dim(expr_se)

# ---- Normalize via DESeq2's variance-stabilizing transformation ----
# Raw counts aren't directly comparable across samples because sequencing
# depth varies per sample. VST puts all samples on a comparable scale
# while stabilizing variance across the expression range, producing
# values suitable for downstream modeling (as opposed to raw counts,
# which most models, including glmnet, assume are not).
#
# design = ~ 1 (no covariates): DESeq2 is built for differential
# expression testing between groups, but it's being used here purely for
# its normalization machinery, not for any group comparison — a null
# design is the correct choice for that use case.
#
# blind = TRUE: dispersion estimation ignores any design/grouping
# information, which is the conservative choice since this data will be
# reused for modeling afterward.
dds <- DESeqDataSetFromMatrix(
  countData = round(assay(expr_se, "unstranded")),
  colData   = as.data.frame(colData(expr_se)),
  design    = ~ 1
)
# vst: variance-stabilizing transformation: flattens variance between samples; returns a DESeqTransform object
# assay() extracts just the numeric values out of it, genes × samples
vst_mat <- assay(vst(dds, blind = TRUE))

# ---- Align expression matrix to cohort patient order ----
# vst_mat's columns are sample barcodes; truncate to patient IDs to match
# cohort$patient_id, then reorder columns to exactly match cohort's row
# order. This lets downstream code index vst_mat and cohort by position
# rather than re-joining every time, but only if the alignment is exact —
# hence the check.
colnames(vst_mat) <- substr(colnames(vst_mat), 1, 12)
# Confirm every cohort patient has an expression column, in the same order
vst_mat <- vst_mat[, cohort$patient_id]

# ---- Pre-filter to most variable genes, using TRAINING DATA ONLY ----
# 16,587 genes is still far too many for ~97 training events to support.
# Reduce to the 5,000 most variable genes before elastic-net modeling.
#
# This filter is computed only on training samples, not the full cohort —
# critical for avoiding leakage. Unlike filtering by correlation with
# survival (which would require the outcome and must happen inside each
# CV fold), filtering by variance alone is outcome-blind: it never looks
# at os_days or os_event, only at how much a gene varies across patients.
# That makes it safe to compute once here rather than inside every fold,
# but it's still restricted to training data only, so no information from
# the test set's expression values leaks into gene selection.
train_ids <- cohort$patient_id[cohort$split == "train"]
vst_train <- vst_mat[, train_ids]

gene_var <- apply(vst_train, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]

vst_train_filtered <- vst_train[top_genes, ]
dim(vst_train_filtered)

# ---- Save outputs ----
# vst_mat: full normalized matrix (16,587 genes x 508 patients, both
#   train and test) — kept in full in case a different gene filter or
#   filter size is wanted later without re-running VST from scratch.
# top_genes: the 5,000 gene IDs selected by the training-only variance
#   filter — needed to subset the test set identically at evaluation time.
# vst_train_filtered: the actual training matrix (5,000 genes x 404
#   patients) that the elastic-net Cox model will be fit on.
saveRDS(vst_mat, here("data", "processed", "vst_mat.rds"))
saveRDS(top_genes, here("data", "processed", "top_genes.rds"))
saveRDS(vst_train_filtered, here("data", "processed", "vst_train_filtered.rds"))
