# Molecular vs. histologic risk stratification in low-grade glioma

## Motivation

Low-grade glioma (LGG) was historically graded by histology alone. Ceccarelli
et al. (2016) showed that molecular features *IDH1/IDH2* mutation status and
1p/19q codeletion stratify survival better than histologic grade, a
finding substantial enough that the WHO revised its official classification
of glioma in 2016 to require molecular testing alongside histology.

This project asks three nested questions:

1. Do IDH status and 1p/19q codeletion stratify overall survival better than
   histologic grade alone, in this cohort? *(Reproduces the known finding —
   the sanity check that the data pipeline is correct.)*
2. Does a Cox model built on RNA-seq expression features add prognostic value
   on top of IDH/1p19q status?
3. Does that combined model hold up under nested cross-validation and,
   ideally, external validation?

## Data

TCGA-LGG, pulled via `TCGAbiolinks` from GDC: clinical/survival data, curated
molecular subtype calls (IDH, 1p/19q), and RNA-seq expression
(STAR-Counts, harmonized GDC pipeline).

## Methods

**Cohort.** Primary tumor samples only (n=516). Patients were excluded for:
missing IDH/1p19q molecular subtype call (n=3), unknown vital status (n=1),
missing death date among deceased patients (n=1), and invalid follow-up time
(n=1 negative value), yielding a final analysis cohort of **508 patients,
123 deaths (24% event rate)**.

Grade and survival fields were taken from GDC's harmonized clinical table;
molecular subtype calls (IDH status, 1p/19q codeletion) were taken from the
curated PanCancer Atlas publication fields, to match the original
classification literature. Grade concordance between the two sources was
verified at 100% (456/456 patients with grade in both).

**Split.** Patients were split 80/20 into train/test, stratified jointly on
IDH/codel subtype and event status (seed = 42). Train: n=404, 97 events.
Test: n=104, 26 events. The test set was not examined prior to or during
model development.

**Survival analysis.** Kaplan-Meier curves and log-rank tests compared
overall survival across histologic grade (G2/G3) and molecular subtype.
Cox proportional hazards models quantified effect sizes as hazard ratios
with 95% CIs; discrimination was assessed via Harrell's concordance index.
The proportional hazards assumption was tested via Schoenfeld residuals
(`cox.zph`). All modeling in this section was performed on the training
set only.

## Results — Hypothesis 1: Molecular subtype vs. histologic grade

Consistent with Ceccarelli et al. (2016), molecular subtype stratified
overall survival far more sharply than histologic grade.

| Model | C-index | SE |
|---|---|---|
| Histologic grade (G2/G3) | 0.653 | 0.025 |
| Molecular subtype (IDH/codel) | 0.741 | 0.032 |
| Subtype + grade + age | 0.817 | 0.023 |

Log-rank test for molecular subtype: χ² = 106, df = 2, p < 2×10⁻¹⁶.
Log-rank test for grade: χ² = 27.9, df = 1, p = 1×10⁻⁷.

Relative to IDH-wildtype (reference), hazard ratios from the subtype-only
Cox model:

| Subtype | HR | 95% CI | Events / N |
|---|---|---|---|
| IDH-mutant, non-codeleted | 0.168 | (0.106, 0.265) | 41 / 197 |
| IDH-mutant, codeleted | 0.107 | (0.059, 0.193) | 16 / 132 |
| IDH-wildtype (reference) | 1.00 | — | 40 / 75 |

IDH-wildtype tumors — despite being histologically graded as low-grade
(WHO grade II/III) — showed survival behavior consistent with glioblastoma:
53% mortality in this cohort vs. 21% (non-codel) and 13% (codel) in the
IDH-mutant groups.

![Kaplan-Meier by molecular subtype](reports/figures/km_idh_subtype.png)

**Molecular subtype and histologic grade are complementary, not redundant.**
The jump from subtype alone (0.741) to the combined model (0.817) is nearly
as large as the jump from grade alone to subtype alone, indicating grade
retains independent prognostic value even after accounting for molecular
status — mirroring current WHO practice of integrated diagnosis.

**Proportional hazards.** The PH assumption was violated for IDH/codel
status (χ² = 14.9, df = 2, p = 6×10⁻⁴), and for all covariates in the
combined model (global p = 0.002). Schoenfeld residuals show the survival
difference between IDH-wildtype and IDH-mutant groups is most pronounced in
the first ~2 years and attenuates thereafter — consistent with
IDH-wildtype's rapid early mortality: by day 2000, only 1 of 75
IDH-wildtype patients remained under observation. Reported hazard ratios
should be read as an average effect over follow-up rather than a constant
instantaneous risk; the qualitative conclusion is well-supported regardless.

## Limitations

Median follow-up in this cohort is under 2 years, despite a small number of
patients followed beyond 15 years. Risk tables accompanying all KM plots
should be consulted before interpreting survival estimates beyond ~3 years,
where patient counts in some subgroups drop into single digits.

## App

*A Shiny app serving the final model will be linked here once deployed.*

## Reproducing this

See [SETUP.md](SETUP.md) for environment setup (`renv` + Bioconductor).