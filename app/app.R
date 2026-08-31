library(shiny)
library(survival)
library(ggplot2)

# Load the trained model and reference cohort
# Instead of: final_fit <- readRDS("model/final_cox_model.rds")
# Refit directly here so the model's call/environment is self-contained
cohort    <- readRDS("model/cohort.rds")
train     <- cohort[cohort$split == "train", ]

final_fit <- coxph(
  Surv(os_days, os_event) ~ idh_codel + tumor_grade + age_at_index,
  data = train
)

# Reported validation numbers, computed in 06_final_evaluation.R —
# hardcoded here since these are fixed, published results, not
# something the app recomputes live.
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
    
    patient_surv <- survfit(final_fit, newdata = patient)
    
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
