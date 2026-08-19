# Hypothesis 2: Does an expression-derived risk score add prognostic value
# beyond molecular subtype + grade + age?

# Inner loop tries a range of lambda (coefficients for shriking gene coefficients) 
# values and picks whichever scores best. All done with glmnet in one line

# Outer loop, 5 folds, produces the actual score like a typical cross validation

library(glmnet)
library(survival)
library(here)

proc_dir <- here("data", "processed")

# cohort table has clinical variables and split assignment
# VST matrix has 5,000 most variable genes
cohort              <- readRDS(file.path(proc_dir, "cohort.rds"))
vst_train_filtered  <- readRDS(file.path(proc_dir, "vst_train_filtered.rds"))

train <- cohort[cohort$split == "train", ]

# Guardrail: row i of `train` must correspond to column i of
# `vst_train_filtered`, or every downstream indexing operation
# silently pairs the wrong patient with the wrong expression
# profile. Fail loudly here rather than getting subtly wrong
# results later.
stopifnot(identical(train$patient_id, colnames(vst_train_filtered)))

set.seed(123) # same fold assignment every run
k_outer <- 5

# Stratify fold assignment by event status (death vs censored),
# so each outer fold has a comparable number of events. Without
# this, a fold could randomly end up with very few deaths,
# making that fold's performance estimate unstable.
outer_strata <- train$os_event
outer_folds <- rep(NA_integer_, nrow(train))
for (s in unique(outer_strata)) {
  ix <- which(outer_strata == s)
  outer_folds[ix] <- sample(rep_len(1:k_outer, length(ix)))
}

table(outer_folds, train$os_event)

# The outer loop
outer_results <- vector("list", k_outer)

for (fold in 1:k_outer) {
  message(sprintf("Outer fold %d of %d", fold, k_outer))
  
  test_idx  <- which(outer_folds == fold)
  train_idx <- which(outer_folds != fold)
  
  inner_train      <- train[train_idx, ]
  inner_test       <- train[test_idx, ]
  expr_inner_train <- vst_train_filtered[, train_idx] # 5000 x 82 
  expr_inner_test  <- vst_train_filtered[, test_idx] # 5000 x 322
  
  # ---- Baseline model: subtype + grade + age ----
  # This is the model Hypothesis 2 has to beat — the same
  # combined model that scored C = 0.817 in Hypothesis 1,
  # refit here on this fold's inner-train subset.
  baseline_fit <- coxph(
    Surv(os_days, os_event) ~ idh_codel + tumor_grade + age_at_index,
    data = inner_train
  )
  # Linear predictor (risk score) for the held-out outer-test
  # patients, from a model that never saw them during fitting.
  baseline_risk_test <- predict(baseline_fit, newdata = inner_test, type = "lp")
  
  # ---- Expression model: elastic-net Cox, tuned via INNER cv.glmnet ----
  # glmnet expects samples-as-rows, genes-as-columns — the
  # transpose of how vst_train_filtered is stored (genes x samples).
  x_inner_train <- t(expr_inner_train)   # glmnet wants samples x genes
  y_inner_train <- Surv(inner_train$os_days, inner_train$os_event)
  
  # This single call performs the INNER cross-validation: it
  # internally splits inner_train into 5 folds, fits the
  # elastic-net path across a range of lambda values, and
  # scores each lambda by concordance — returning the
  # best-performing lambda. This is the mechanism that avoids
  # picking lambda based on a biased, self-serving score.
  cv_fit <- cv.glmnet(
    x_inner_train, y_inner_train,
    family = "cox",
    alpha  = 0.5,        # elastic-net: midpoint between ridge (0) and lasso (1)
    nfolds = 5,           # this IS the inner loop
    type.measure = "C"    # optimize for concordance
  )
  
  best_lambda <- cv_fit$lambda.min
  
  # Expression risk score for the held-out outer-test patients
  # Collapse the 5,000-gene elastic-net fit into a single
  # per-patient risk score (the linear predictor), for both
  # inner-train (needed to fit the extended Cox model next)
  # and outer-test (needed to evaluate it).
  x_inner_test <- t(expr_inner_test)
  expr_score_test <- as.numeric(predict(cv_fit, newx = x_inner_test,
                                        s = best_lambda, type = "link"))
  expr_score_train <- as.numeric(predict(cv_fit, newx = x_inner_train,
                                         s = best_lambda, type = "link"))
  
  # Attach the expression score as a new covariate.
  inner_test$expr_score <- expr_score_test
  inner_train$expr_score <- expr_score_train
  
  # ---- Extended model: baseline linear predictor + expression score ----
  # Same baseline formula, plus the single expression-derived
  # score. This is how 5,000 genes get folded into an ordinary
  # coxph() without trying to fit 5,000 individual coefficients
  # against ~80 events.
  extended_fit <- coxph(
    Surv(os_days, os_event) ~ idh_codel + tumor_grade + age_at_index + expr_score,
    data = inner_train %>% dplyr::mutate(expr_score = as.numeric(predict(
      cv_fit, newx = x_inner_train, s = best_lambda, type = "link")))
  )
  extended_risk_test <- predict(extended_fit, newdata = inner_test, type = "lp")
  
  # ---- Score both models on the SAME outer-test fold ----
  # Paired comparison: same patients, same fold, both models
  # evaluated identically. This is what makes "did expression
  # help" a fair question rather than two separate experiments.
  c_baseline <- survival::concordance(Surv(inner_test$os_days, inner_test$os_event) ~ I(-baseline_risk_test))$concordance
  c_extended <- survival::concordance(Surv(inner_test$os_days, inner_test$os_event) ~ I(-extended_risk_test))$concordance
  
  # Record this fold's result, including how many genes the
  # elastic-net penalty actually kept nonzero at the chosen
  # lambda — useful context for interpreting the model later.
  outer_results[[fold]] <- tibble::tibble(
    fold        = fold,
    n_test      = nrow(inner_test),
    events_test = sum(inner_test$os_event),
    lambda      = best_lambda,
    n_genes_selected = sum(coef(cv_fit, s = best_lambda) != 0),
    c_baseline  = c_baseline,
    c_extended  = c_extended
  )
}

# Combine all 5 outer-fold results into one table — this is
# your actual answer to Hypothesis 2. Compare c_baseline vs
# c_extended fold by fold: does expression win consistently,
# or is any apparent lift within noise?
results_df <- dplyr::bind_rows(outer_results)
print(results_df)

mean(results_df$c_baseline)   # ~0.801
mean(results_df$c_extended)   # ~0.804
mean(results_df$c_extended - results_df$c_baseline)   # ~0.003

# Expression adds essentially nothing beyond subtype + grade + age, and the answer to Hypothesis 2 is no.

# Your baseline model is already excellent (C ≈ 0.80–0.82). IDH status alone drives an enormous, 
# well-established biological signal — it's close to a master regulator of this cancer's behavior, 
# reflected genome-wide in methylation and, likely, in much of the bulk expression signal too. 
# It would be somewhat surprising if 5,000 variable genes contained substantial additional 
# prognostic information once IDH/codel status, grade, and age are already accounted for — 
# a lot of the expression variance is plausibly downstream of, and therefore redundant with, 
# the molecular subtype you're already modeling directly.

