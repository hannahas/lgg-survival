library(shiny)
library(survival)
library(ggplot2)

# Load the trained model and reference cohort
# These are 15 column tables with sample, codel status, mutation status, tumor grade, days to event
# Instead of: final_fit <- readRDS("model/final_cox_model.rds")
# Refit directly here so the model's call/environment is self-contained
cohort    <- readRDS("model/cohort.rds")
train     <- cohort[cohort$split == "train", ]

# Surv() bundles a follow-up time (os_days) together with a 0/1 flag 
# (os_event, 1 = death, 0 = still alive/censored at last contact) 
# into one combined outcome object with its own class (Surv)
# This lets a Cox model use both kinds of patients — 
# the ones with a known death day and the ones who just contributed "alive at least until here"
#  underneath the Surv class label it's  just a two-column matrix — os_days and os_event

# Model survival as a function of subtype, grade, and age.

# coxph() function runs a maximum partial likelihood estimation. 
# MLE is the probability of seeing your entire data set as a function of model parameters,
# in order to find the parameters that best estimate the observed outcome. 
# Cox models are semi-parameteric because they don't have a shape of the reference baseline,
# so we can't use full MLE, we need maximum PARTIAL likelihood estimation.
# "partial" deliberately throws away information about absolute timing 
# and uses only the relative ordering of who died when. 

# final_fit is an object that contains the coefficients of the estimation
final_fit <- coxph(
  Surv(os_days, os_event) ~ idh_codel + tumor_grade + age_at_index,
  data = train
)

# Variable      Subgroup                    Role                 HR
# ------------  --------------------------  -------------------  -----
# idh_codel     IDH-wildtype                Reference            1.00
#               IDH-mutant, non-codeleted   vs. reference        0.316
#               IDH-mutant, codeleted       vs. reference        0.155
# tumor_grade   Grade II                    Reference            1.00
#               Grade III                   vs. reference        1.905

# Reported validation numbers, computed in 06_final_evaluation.R —
# hardcoded here since these are fixed, published results, not
# something the app recomputes live.

# Concordance indices for the training and test data. 0.5 = random, 1.0 = perfect.
C_TRAIN <- 0.817
C_TEST  <- 0.902


ui <- fluidPage(
  titlePanel("LGG Survival Risk Estimator"),
  
  fluidRow(
    column(
      width = 4,
      wellPanel(
        h4("Patient profile"),
        selectInput(
          "idh_codel", "Molecular subtype",
          choices = c("IDH-wildtype" = "IDHwt",
                      "IDH-mutant, non-codeleted" = "IDHmut-non-codel",
                      "IDH-mutant, codeleted" = "IDHmut-codel"),
          selected = "IDHmut-non-codel"
        ),
        selectInput(
          "tumor_grade", "Histologic grade",
          choices = c("Grade II" = "G2", "Grade III" = "G3"),
          selected = "G2"
        ),
        numericInput(
          "age", "Age at diagnosis",
          value = 45, min = 18, max = 90
        ),
        actionButton("predict", "Predict", class = "btn-primary")
      ),
      
      wellPanel(
        h4("Model validation"),
        p(strong("Training set concordance: "), textOutput("c_train_display", inline = TRUE)),
        p(strong("Held-out test set concordance: "), textOutput("c_test_display", inline = TRUE)),
        p(em("This is a research demonstration, not a validated clinical
              tool. See the project README for full methodology and
              limitations."))
      )
    ),
    
    column(
      width = 8,
      plotOutput("km_plot", height = "500px")
    )
  )
)


server <- function(input, output, session) {
  
  output$c_train_display <- renderText(sprintf("%.3f", C_TRAIN))
  output$c_test_display  <- renderText(sprintf("%.3f", C_TEST))
  
  # new_patient is a function object. See new_patient() in the next block.
  # eventReactive works whenever predict is clicked.
  # It creates the new_patient function that has the selected variables baked in.
  new_patient <- eventReactive(input$predict, {
    data.frame(
      idh_codel    = factor(input$idh_codel,
                            levels = levels(train$idh_codel)),
      tumor_grade  = factor(input$tumor_grade,
                            levels = levels(train$tumor_grade)),
      age_at_index = input$age
    )
  })
  
  output$km_plot <- renderPlot({
    req(new_patient())
    patient <- new_patient()
    
    # Builds the model's prediction curve
    # It builds a baseline curve from the training data. 
    # The baseline curve is a model-based, weighted construction using every patient in the training set, 
    # where reference-level patients get full weight (exp(0) = 1) and everyone else contributes according 
    # to their own hazard ratio. 
    
    # patient_surv is a fitted survival curve object, a complete, personalized
    # survival estimate for one hypothetical patient, uilt by combining a fixed, 
    # already-trained model with one new patient's covariate values.
    
    # patient_surv is an R list 
    # patient_surv contains a hazard ratio, one number applied to the whole KM curve time.
    patient_surv <- survfit(final_fit, newdata = patient)
    
    # S_patient(t) = S_baseline(t)^HR
    
    # Gray curve:
    # "probability a member of this group survives past time t" — built from, 
    # but not identical to, simply plotting the raw individual events.
    
    # Red curve:
    # The already-existing, already-fit model coxph() (fit once) gets queried with 
    # this specific patient's X/Y/Z values using survfit() (apply per-patient),
    # and returns a prediction, namely, the hazard ratio (one number). 
    # The hazard ratio gets applied to S_baseline(t) — a fixed sequence that already 
    # exists, computed once from the whole training set, completely independent of this patient — 
    # via S_baseline(t)^HR, at every single time point in that sequence. 
    
    
    # Reference group narrowed to matching subtype AND grade (Option A)
    ref_patients <- train[train$idh_codel == input$idh_codel &
                            train$tumor_grade == input$tumor_grade, ]
    n_ref <- nrow(ref_patients)
    events_ref <- sum(ref_patients$os_event)
    
    subtype_km <- survfit(Surv(os_days, os_event) ~ 1, data = ref_patients)
    
    ref_label <- sprintf("Training cohort: %s, %s (n=%d, %d events)",
                         input$idh_codel, input$tumor_grade,
                         n_ref, events_ref)
    
    pred_df <- data.frame(
      time = patient_surv$time,
      surv = patient_surv$surv,
      group = "Predicted (this patient)"
    )
    
    ref_df <- data.frame(
      time = subtype_km$time,
      surv = subtype_km$surv,
      group = ref_label
    )
    
    plot_df <- rbind(pred_df, ref_df)
    
    ggplot(plot_df, aes(x = time, y = surv, color = group)) +
      geom_step(linewidth = 1) +
      scale_color_manual(values = setNames(
        c("#D55E00", "grey60"),
        c("Predicted (this patient)", ref_label)
      )) +
      labs(
        x = "Days from diagnosis",
        y = "Overall survival probability",
        color = NULL,
        title = "Predicted survival for entered patient profile",
        subtitle = if (n_ref < 20) {
          "Note: reference group is small (n<20); interpret with caution"
        } else {
          "Compared against training patients of matching subtype and grade"
        }
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)
