## =============================================================================
## Plasma p-tau217 prescreening for Alzheimer's disease trials
## Thresholds -> sample size -> recruitment time -> cost
##
## =============================================================================

library(pROC)

## -----------------------------------------------------------------------------
## 1. THRESHOLDS Calculation
##    liberal      : threshold whose sensitivity is 95%
##    intermediate : Youden index 
##    conservative : threshold whose specificity is 95%
## -----------------------------------------------------------------------------
define_thresholds <- function(score, status, target = 0.95) {
  ok <- is.finite(score) & status %in% c(0, 1)
  score <- score[ok]; status <- status[ok]
  
  roc_obj <- pROC::roc(response = status, predictor = score,
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  co <- pROC::coords(roc_obj, x = "all",
                     ret = c("threshold", "sensitivity", "specificity"),
                     transpose = FALSE)
  co <- co[is.finite(co$threshold), ]
  
  c(liberal      = co$threshold[which.min(abs(co$sensitivity - target))],
    intermediate = co$threshold[which.max(co$sensitivity + co$specificity - 1)],
    conservative = co$threshold[which.min(abs(co$specificity - target))])
}

## Sensitivity and specificity of each threshold (positive if score > threshold)
threshold_performance <- function(score, status, thresholds) {
  ok <- is.finite(score) & status %in% c(0, 1)
  score <- score[ok]; status <- status[ok]
  do.call(rbind, lapply(names(thresholds), function(nm) {
    th <- thresholds[[nm]]
    data.frame(strategy    = nm,
               threshold   = th,
               sensitivity = mean(score[status == 1] >  th),
               specificity = mean(score[status == 0] <= th))
  }))
}

## -----------------------------------------------------------------------------
## 2. ASSUMPTIONS (Prevalence, cost, recruitment time)
## -----------------------------------------------------------------------------

# Target number of Abeta-PET-positive participants to enroll
TP_REQUIRED <- 1000

# Abeta-PET-positivity prevalence in the screened population.
PREVALENCE <- c(CU = 0.30)

# Costs (US$). 
COSTS <- list(pet     = 3000,   # Abeta PET, per scan
              plasma  = 200,    # plasma biomarker, per participant screened
              recruit = 1000)   # recruitment, consenting and clinical assessment,

# Recruitment time.
REFERENCE_TRIALS <- data.frame(
  Screened = c(6757, 6173, 3736, 3664),
  Included = c(1643, 1653,  813,  806),
  Duration = c(2.8,  2.9,   1.8,  1.5)     # years
)

## -----------------------------------------------------------------------------
## 3. SAMPLE SIZE, TIME AND COST
## -----------------------------------------------------------------------------
recruitment_model <- function(sensitivity, specificity, prevalence,
                              tp_required = TP_REQUIRED,
                              trials = REFERENCE_TRIALS,
                              costs = COSTS) {
  
  # Time of recuritment
  mean_rate      <- mean(trials$Included / trials$Duration)  
  mean_inclusion <- mean(trials$Included / trials$Screened)  
  years_per_screened <- mean_inclusion / mean_rate
  
  n_step1 <- ceiling(tp_required / (sensitivity * prevalence))
  n_step2 <- ceiling((sensitivity * prevalence +
                        (1 - specificity) * (1 - prevalence)) * n_step1)
  
  data.frame(
    sensitivity      = sensitivity,
    specificity      = specificity,
    prevalence       = prevalence,
    n_prescreened    = n_step1,                       # plasma tests
    n_pet            = n_step2,                       # confirmatory Abeta PET scans
    time_years       = n_step1 * years_per_screened,
    cost_usd         = n_step1 * costs$plasma +
      n_step2 * (costs$pet + costs$recruit)
  )
}

## Reference: Abeta PET alone
pet_only_model <- function(prevalence,
                           tp_required = TP_REQUIRED,
                           trials = REFERENCE_TRIALS,
                           costs = COSTS) {
  mean_rate      <- mean(trials$Included / trials$Duration)
  mean_inclusion <- mean(trials$Included / trials$Screened)
  years_per_screened <- mean_inclusion / mean_rate
  
  n <- ceiling(tp_required / prevalence)
  data.frame(
    n_prescreened = n,
    n_pet         = n,
    time_years    = n * years_per_screened,
    cost_usd      = n * (costs$pet + costs$recruit)
  )
}

## All three strategies at once, plus the % change versus Abeta PET alone
compare_strategies <- function(thresholds_perf, prevalence, ...) {
  res <- do.call(rbind, lapply(seq_len(nrow(thresholds_perf)), function(i) {
    cbind(strategy = thresholds_perf$strategy[i],
          threshold = thresholds_perf$threshold[i],
          recruitment_model(thresholds_perf$sensitivity[i],
                            thresholds_perf$specificity[i],
                            prevalence, ...))
  }))
  ref <- pet_only_model(prevalence, ...)
  res$sample_size_change_pct <- (res$n_prescreened / ref$n_prescreened - 1) * 100
  res$pet_change_pct         <- (res$n_pet         / ref$n_pet         - 1) * 100
  res$cost_change_pct        <- (res$cost_usd      / ref$cost_usd      - 1) * 100
  res$time_difference_years  <-  res$time_years - ref$time_years
  attr(res, "pet_only") <- ref
  res
}

## Print the results
show_assumptions <- function(prevalence,
                             tp_required = TP_REQUIRED,
                             costs = COSTS,
                             trials = REFERENCE_TRIALS) {
  cat("Assumptions:\n")
  cat("  Target true positives to enroll :", tp_required, "\n")
  cat("  Abeta PET prevalence            :", prevalence, "\n")
  cat("  Cost per Abeta PET scan         : $", costs$pet, "\n", sep = "")
  cat("  Cost per plasma test            : $", costs$plasma, "\n", sep = "")
  cat("  Cost recruitment/assessment     : $", costs$recruit, "\n", sep = "")
  cat("  Reference (n)            :", nrow(trials), "\n\n")
}

