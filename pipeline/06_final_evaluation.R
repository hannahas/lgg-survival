library(survival)
library(here)

proc_dir <- here("data", "processed")
cohort   <- readRDS(file.path(proc_dir, "cohort.rds"))

train <- cohort[cohort$split == "train", ]
test  <- cohort[cohort$split == "test", ]

stopifnot(nrow(train) + nrow(test) == nrow(cohort))
stopifnot(!any(train$patient_id %in% test$patient_id))

# ---- Fit the final model on ALL training data ----
# This is the only time the full training set is used for a
# single fit rather than being split further — appropriate now
# because model selection is finished (Hypothesis 2 established
# that expression does not improve on this baseline).
final_fit <- coxph(
  Surv(os_days, os_event) ~ idh_codel + tumor_grade + age_at_index,
  data = train
)
summary(final_fit)

# ---- Evaluate ONCE on the sealed test set ----
test_risk <- predict(final_fit, newdata = test, type = "lp")

c_test <- survival::concordance(
  Surv(test$os_days, test$os_event) ~ I(-test_risk)
)
print(c_test)

# ---- Compare to internal (training) performance, for reference ----
c_train <- summary(final_fit)$concordance
c_train

# Compare subtype composition and event distribution: test vs train
table(train$idh_codel, train$os_event)
table(test$idh_codel, test$os_event)

# Proportions, side by side
round(prop.table(table(train$idh_codel, train$os_event), 1), 3)
round(prop.table(table(test$idh_codel, test$os_event), 1), 3)

# 95% CI for each
c(0.902 - 1.96*0.0232, 0.902 + 1.96*0.0232)   # test
c(0.817 - 1.96*0.0226, 0.817 + 1.96*0.0226)   # train