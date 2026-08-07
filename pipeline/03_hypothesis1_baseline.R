library(survival)
library(survminer)
library(here)

proc_dir <- here("data", "processed")
fig_dir  <- here("reports", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cohort <- readRDS(file.path(proc_dir, "cohort.rds"))

# TEST SET IS SEALED — training data only
train <- cohort[cohort$split == "train", ]
nrow(train)

surv_obj <- Surv(time = train$os_days, event = train$os_event)

# Fit a Kaplan-Meier curve
km_idh <- survfit(surv_obj ~ idh_codel, data = train)
print(km_idh)

# Run the log-rank test
# at every timepoint where a death occurs anywhere in the cohort, 
# it compares how many deaths actually happened in each group against 
# how many you'd expect if all three groups had identical survival — 
# pooling everyone together and allocating deaths proportional to group size at risk. 
# Sum that comparison across every event time, and you get a chi-squared statistic. 
# Large deviations between observed and expected deaths, 
# consistently in the same direction, produce a large statistic and a small p-value.

# log-rank result
# tells you the three curves differ, with high confidence, and nothing about magnitude or direction for a specific pair.
survdiff(surv_obj ~ idh_codel, data = train)

p_idh <- ggsurvplot(
  km_idh,
  data        = train,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = TRUE,
  xlab        = "Days from diagnosis",
  ylab        = "Overall survival probability",
  legend.title = "Molecular subtype",
  legend.labs = c("IDH-wildtype", "IDH-mut non-codel", "IDH-mut codel"),
  palette     = c("#D55E00", "#E69F00", "#0072B2")
)

# KM plot
# discrete step every time a patient dies. Denominator is # of patients left in the cohort (this changes over time).
# tick marks are alive at last contact, follow-up ended without observed death

# IDH-WT: no mutation, graded as low-grade even though they behave like a more aggressive glioblastoma
# IDH-mut non-codeleted: astrocytomas
# IDH-mut codeleted: whole-arm loss of chr 1p and 19q. oligodendrogliomas, good prognosis, respond to chemo & radiation

# molecular subtype is better predictor than histology

# Table
# 

print(p_idh)
ggsave(file.path(fig_dir, "km_idh_subtype.png"), plot = p_idh$plot,
       width = 8, height = 7, dpi = 300)

km_grade <- survfit(surv_obj ~ tumor_grade, data = train)
print(km_grade)
survdiff(surv_obj ~ tumor_grade, data = train)

cox_grade <- coxph(surv_obj ~ tumor_grade, data = train)
cox_idh   <- coxph(surv_obj ~ idh_codel,   data = train)
cox_both  <- coxph(surv_obj ~ idh_codel + tumor_grade + age_at_index, data = train)

# Concordance: C-index, the discrimination measure
summary(cox_grade)$concordance # pathologist, better than chance but not by much
summary(cox_idh)$concordance # answer to Hypothesis 1: IDH and 1p/19q status discriminate survival substantially better than histologic grade does, in this cohort
summary(cox_both)$concordance # integrated diagnosis using both molecular and histologic information outperforms either alone

# Hazard ratios (exp(coef)) non-codel patient risk is 17% of the WT patient at any given time
# A codel patient's risk is 11% of the WT patient's

summary(cox_idh)

# "Molecular subtype alone achieves better discrimination than histologic grade (C=0.741 vs. 0.653), 
# but the two are complementary rather than redundant — a combined model incorporating both achieves C=0.817, 
# consistent with current WHO practice of using integrated molecular and histologic diagnosis."


#############################################################################

# The proportional hazards assumption is violated for IDH/codel status (p = 6×10⁻⁴), w
# ith the survival difference between IDH-wildtype and IDH-mutant groups being most pronounced 
# in the first ~2 years and attenuating thereafter. This is consistent with IDH-wildtype tumors' 
# rapid early mortality — by day 2000, only 1 of 75 IDH-wildtype patients remained under observation (see risk table), 
# so late-follow-up comparisons for this group reflect a small, non-representative subset of early survivors. 
# The reported hazard ratios should be interpreted as an average effect over the full follow-up period rather than 
# a constant instantaneous risk; the qualitative conclusion — IDH-wildtype carries substantially worse prognosis — 
# is well-supported by both the log-rank test and the magnitude of early separation, 
# but should not be read as constant across all follow-up times.

zph_idh <- cox.zph(cox_idh)
print(zph_idh)
plot(zph_idh)

zph_both <- cox.zph(cox_both)
print(zph_both)

