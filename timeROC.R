
# ============================================================
# Time-dependent ROC analysis for survival outcomes
# ============================================================
#
# Description:
# Fits Cox proportional hazards models and evaluates their
# discrimination using time-dependent ROC curves and AUC.
#
# Required packages:
#   survival, timeROC
#
# Input data:
#   One row per participant, containing:
#   - follow_up_time: follow-up duration (months)
#   - event_status: event indicator (0 = censored, 1 = event)
#   - predictor_1, predictor_2, etc.: model predictors
#
# Note:
# This script is a template. Replace the example variable
# names with those in your dataset.
# ============================================================

library(survival)
library(timeROC)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

# Example:
# data <- read.csv("data/analysis_data.csv")

# Define survival outcome
time_var  <- "follow_up_time"
event_var <- "event_status"

# Define candidate models
# Each list entry contains the predictors for one Cox model.
model_formulas <- list(
  Model_1 = ~ predictor_1,
  Model_2 = ~ predictor_2,
  Model_3 = ~ predictor_1 + predictor_2,
  Model_4 = ~ predictor_3,
  Model_5 = ~ predictor_1 + predictor_3
)

# Evaluation times (same units as follow-up)
evaluation_times <- c(12, 24, 48, 72, 96, 108, 120)

# Time point for ROC comparison
plot_time <- 72

# ------------------------------------------------------------
# 2. Prepare analysis dataset
# ------------------------------------------------------------

required_vars <- unique(c(
  time_var,
  event_var,
  unlist(lapply(model_formulas, all.vars))
))

# Use the same participants across all models
analysis_data <- data[
  complete.cases(data[, required_vars, drop = FALSE]),
  ,
  drop = FALSE
]

if (nrow(analysis_data) == 0) {
  stop("No complete observations available.")
}

if (!all(analysis_data[[event_var]] %in% c(0, 1))) {
  stop("Event status must be coded as 0 or 1.")
}

if (any(analysis_data[[time_var]] < 0)) {
  stop("Follow-up times must be non-negative.")
}

# ------------------------------------------------------------
# 3. Fit Cox proportional hazards models
# ------------------------------------------------------------

survival_formula <- function(predictor_formula) {
  as.formula(
    paste(
      "Surv(", time_var, ",", event_var, ")",
      deparse(predictor_formula)
    )
  )
}

cox_models <- lapply(model_formulas, function(f) {
  coxph(
    survival_formula(f),
    data = analysis_data,
    x = TRUE,
    y = TRUE
  )
})

# Concordance statistics
concordance_results <- do.call(
  rbind,
  lapply(names(cox_models), function(name) {
    model <- cox_models[[name]]
    c_index <- summary(model)$concordance
    
    data.frame(
      Model = name,
      C_index = unname(c_index[1]),
      SE = unname(c_index[2])
    )
  })
)

print(concordance_results)

# ------------------------------------------------------------
# 4. Generate model predictions
# ------------------------------------------------------------

# Cox linear predictors: higher values indicate higher risk
predictions <- lapply(cox_models, function(model) {
  predict(
    model,
    newdata = analysis_data,
    type = "lp"
  )
})

# ------------------------------------------------------------
# 5. Calculate time-dependent ROC curves
# ------------------------------------------------------------

roc_results <- lapply(predictions, function(marker) {
  timeROC(
    T = analysis_data[[time_var]],
    delta = analysis_data[[event_var]],
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = evaluation_times,
    ROC = TRUE,
    iid = TRUE
  )
})

# Extract AUC values
auc_results <- data.frame(
  Time = roc_results[[1]]$times
)

for (name in names(roc_results)) {
  auc_results[[name]] <- roc_results[[name]]$AUC
}

print(auc_results)

# ------------------------------------------------------------
# 6. Compare time-dependent AUC between models
# ------------------------------------------------------------

# Example: compare Model 1 and Model 3
auc_comparison <- compare(
  roc_results[["Model_1"]],
  roc_results[["Model_3"]]
)

print(auc_comparison)

# ------------------------------------------------------------
# 7. Plot ROC curves at a selected time point
# ------------------------------------------------------------

model_names <- names(roc_results)
plot_colors <- grDevices::hcl.colors(
  length(model_names),
  palette = "Dark 3"
)

pdf("time_dependent_ROC.pdf", width = 8, height = 7)

plot(
  roc_results[[1]],
  time = plot_time,
  col = plot_colors[1],
  lwd = 2,
  title = FALSE
)

if (length(roc_results) > 1) {
  for (i in 2:length(roc_results)) {
    plot(
      roc_results[[i]],
      time = plot_time,
      col = plot_colors[i],
      lwd = 2,
      add = TRUE,
      title = FALSE
    )
  }
}

title(
  main = paste("Time-dependent ROC at", plot_time, "months")
)

legend(
  "bottomright",
  legend = model_names,
  col = plot_colors,
  lwd = 2,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 8. Plot AUC over time
# ------------------------------------------------------------

pdf("time_dependent_AUC.pdf", width = 8, height = 7)

plotAUCcurve(
  roc_results[[1]],
  FP = 2,
  col = plot_colors[1],
  conf.int = FALSE
)

if (length(roc_results) > 1) {
  for (i in 2:length(roc_results)) {
    plotAUCcurve(
      roc_results[[i]],
      FP = 2,
      add = TRUE,
      col = plot_colors[i],
      conf.int = FALSE
    )
  }
}

legend(
  "bottomright",
  legend = model_names,
  col = plot_colors,
  lwd = 2,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 9. Export results
# ------------------------------------------------------------

write.csv(
  concordance_results,
  "concordance_results.csv",
  row.names = FALSE
)

write.csv(
  auc_results,
  "time_dependent_AUC_results.csv",
  row.names = FALSE
)
