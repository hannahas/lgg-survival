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

km_idh <- survfit(surv_obj ~ idh_codel, data = train)
print(km_idh)

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

print(p_idh)
ggsave(file.path(fig_dir, "km_idh_subtype.png"), print(p_idh),
       width = 8, height = 7, dpi = 300)

km_grade <- survfit(surv_obj ~ tumor_grade, data = train)
print(km_grade)
survdiff(surv_obj ~ tumor_grade, data = train)

cox_grade <- coxph(surv_obj ~ tumor_grade, data = train)
cox_idh   <- coxph(surv_obj ~ idh_codel,   data = train)
cox_both  <- coxph(surv_obj ~ idh_codel + tumor_grade + age_at_index, data = train)

summary(cox_grade)$concordance
summary(cox_idh)$concordance
summary(cox_both)$concordance

summary(cox_idh)