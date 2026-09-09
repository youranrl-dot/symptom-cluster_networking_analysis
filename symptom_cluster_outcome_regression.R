## ============================================================
## Section 3.7 -- Predictive Validity of Symptom Clusters (re-run, N=1,496)
##
## Logistic regression of three clinical outcomes on the three Gemini
## symptom-cluster burden scores (count of present symptoms per cluster),
## adjusted for age and sex -- matches the manuscript's Section 3.7 model
## specification, but on the corrected cohort that excludes the 11 patients
## with no detectable CC/HPI text (see symptom_network_analysis.R header).
##
## Requires, in the same working directory:
##   - patient_level_symptom_clusters.csv  (written by symptom_network_analysis.R,
##     the GEMINI version -- re-run that script first so it reflects N=1,496)
##   - crc_patient_outcomes.csv            (in_hospital_mortality, readmit_30d,
##     mortality_1yr, age, gender -- for all 1,507; the merge below restricts
##     to the 1,496 with valid Gemini symptom data)
## ============================================================

# setwd("path/to/your/folder")

clusters <- read.csv("patient_level_symptom_clusters.csv", stringsAsFactors = FALSE)
outcomes <- read.csv("crc_patient_outcomes.csv", stringsAsFactors = FALSE)

## sanity check: clusters should already be N=1,496 (11 excluded) if you
## re-ran symptom_network_analysis.R after the exclusion fix
cat("Patients in cluster file:", nrow(clusters), "\n")

df <- merge(clusters[, c("subject_id", "cluster1_burden", "cluster2_burden", "cluster3_burden")],
            outcomes, by = "subject_id")
cat("Patients after merging with outcomes:", nrow(df), "\n")

df$gender <- relevel(factor(df$gender), ref = "F")  # so the reported OR is for Male vs Female

## helper: fit + print OR / 95% CI / p for each cluster-burden term
run_model <- function(outcome_var, data, label) {
  formula_str <- paste0(outcome_var, " ~ cluster1_burden + cluster2_burden + cluster3_burden + age + gender")
  fit <- glm(as.formula(formula_str), data = data, family = binomial)

  co <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  or_table <- data.frame(
    term = rownames(co),
    OR = round(exp(co[, "Estimate"]), 2),
    CI_low = round(exp(ci[, 1]), 2),
    CI_high = round(exp(ci[, 2]), 2),
    p = signif(co[, "Pr(>|z|)"], 3)
  )

  n_outcome <- sum(data[[outcome_var]])
  pct_outcome <- round(mean(data[[outcome_var]]) * 100, 1)
  cat("\n==============================================\n")
  cat(label, " (n = ", nrow(data), ", events = ", n_outcome,
      " [", pct_outcome, "%])\n", sep = "")
  cat("==============================================\n")
  print(or_table, row.names = FALSE)

  pr2 <- 1 - (fit$deviance / fit$null.deviance)
  cat("Pseudo R-squared (McFadden-ish, 1 - resid/null deviance):", round(pr2, 3), "\n")

  invisible(fit)
}

fit_inhosp   <- run_model("in_hospital_mortality", df, "In-hospital mortality")
fit_readmit  <- run_model("readmit_30d",           df, "30-day readmission")
fit_1yr      <- run_model("mortality_1yr",          df, "1-year mortality")

## ---- Write a combined results table for pasting into the manuscript ------
extract_cluster_rows <- function(fit, outcome_label) {
  co <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  rows <- grep("^cluster", rownames(co))
  data.frame(
    outcome = outcome_label,
    cluster = rownames(co)[rows],
    OR = round(exp(co[rows, "Estimate"]), 2),
    CI_low = round(exp(ci[rows, 1]), 2),
    CI_high = round(exp(ci[rows, 2]), 2),
    p = signif(co[rows, "Pr(>|z|)"], 3)
  )
}

results <- rbind(
  extract_cluster_rows(fit_inhosp,  "In-hospital mortality"),
  extract_cluster_rows(fit_readmit, "30-day readmission"),
  extract_cluster_rows(fit_1yr,     "1-year mortality")
)
cat("\n\n===== Combined cluster-burden results (N=1,496) =====\n")
print(results, row.names = FALSE)

write.csv(results, "section37_cluster_burden_results_N1496.csv", row.names = FALSE)
cat("\nWrote section37_cluster_burden_results_N1496.csv\n")
