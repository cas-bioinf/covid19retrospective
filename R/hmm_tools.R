compute_counterfactual_mortality_log_OR <- function(hmm_fit, numerator_serie_data, denominator_serie_data, day = 28) {
  if(nrow(numerator_serie_data) != nrow(denominator_serie_data)) {
    stop("Groups to compare must have equal number of rows")
  }

  combined_data <- hmm_fit$data
  combined_data$serie_data = rbind(
    numerator_serie_data %>% mutate(.serie = paste0("__num_", .serie)),
    denominator_serie_data %>% mutate(.serie = paste0("__denom_", .serie))
  )
  combined_data$initial_states = rep(hmm_fit$data$initial_states, 2)

  state_probs <- posterior_epred_rect(hmm_fit, newdata = combined_data, method = posterior_epred_state_prob)

  mortality_state <- which(hmm_fit$data$hidden_state_data$id == "Death")
  numerator_indices <- startsWith(combined_data$serie_data$.serie, "__num_")
  denominator_indices <- startsWith(combined_data$serie_data$.serie, "__denom_")
  if(!all(numerator_indices | denominator_indices)) {
    stop("Bad indexing")
  }

  samples_probs_numerator <-
    rowMeans(state_probs[mortality_state, day + 1, , numerator_indices])
  samples_probs_denominator <-
    rowMeans(state_probs[mortality_state, day + 1, , denominator_indices])

  log_odds_numerator <- log(samples_probs_numerator) - log1p(-samples_probs_numerator)
  log_odds_denominator <- log(samples_probs_denominator) - log1p(-samples_probs_denominator)

  log_odds_numerator - log_odds_denominator
}

evaluate_treatment_hypothesis <- function(fit, hypothesis, treatment_column, model_subgroup, adjusted, model_check) {
  counterfactual_took <- serie_data_28 %>% mutate({{ treatment_column }} := 1)
  counterfactual_no <- serie_data_28 %>% mutate({{treatment_column }} := 0)

  log_OR_samples <- compute_counterfactual_mortality_log_OR(fit, counterfactual_took, counterfactual_no, day = 28)

  bayesian_hypothesis_res_from_draws(
    draws = log_OR_samples,
    model = paste0("HMM", model_subgroup),
    estimand = "log(OR)",
    hypothesis = hypothesis,
    adjusted = adjusted,
    model_check = model_check
  )
}

evaluate_all_treatment_hypotheses <- function(fit, model_subgroup, adjusted, model_check) {
  rbind(
    evaluate_treatment_hypothesis(fit, hypotheses$hcq_reduces_death, took_hcq, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check),
    evaluate_treatment_hypothesis(fit, hypotheses$az_reduces_death, took_az, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check),
    evaluate_treatment_hypothesis(fit, hypotheses$favipiravir_reduces_death, took_favipiravir, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check)
  )
}
