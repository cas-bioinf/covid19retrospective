compute_counterfactual_state_log_OR <- function(hmm_fit, numerator_serie_data, denominator_serie_data, target_state, day = 28,
                                                    cl = NULL, cores = parallel::detectCores()) {
  if(nrow(numerator_serie_data) != nrow(denominator_serie_data)) {
    stop("Groups to compare must have equal number of rows")
  }

  combined_data <- hmm_fit$data
  combined_data$serie_data = rbind(
    numerator_serie_data %>% mutate(.serie = paste0("__num_", .serie), is_numerator = TRUE),
    denominator_serie_data %>% mutate(.serie = paste0("__denom_", .serie), is_numerator = FALSE)
  ) %>%
    mutate(.serie = factor(.serie))

  combined_data$initial_states = rep(hmm_fit$data$initial_states, 2)

  target_state_id <- which(hmm_fit$data$hidden_state_data$id == target_state)
  numerator_indices <- unique(as.integer(combined_data$serie_data$.serie)[combined_data$serie_data$is_numerator])
  denominator_indices <- unique(as.integer(combined_data$serie_data$.serie)[!combined_data$serie_data$is_numerator])
  if(!identical(sort(c(numerator_indices, denominator_indices)), 1:length(unique(combined_data$serie_data$.serie)))) {
    stop("Bad indexing")
  }


  state_probs <- posterior_epred_rect(hmm_fit, newdata = combined_data, method = posterior_epred_state_prob, cl = cl, cores = cores)

  samples_probs_numerator <-
    rowMeans(state_probs[target_state_id, day + 1, , numerator_indices])
  samples_probs_denominator <-
    rowMeans(state_probs[target_state_id, day + 1, , denominator_indices])

  log_odds_numerator <- log(samples_probs_numerator) - log1p(-samples_probs_numerator)
  log_odds_denominator <- log(samples_probs_denominator) - log1p(-samples_probs_denominator)

  log_odds_numerator - log_odds_denominator
}

evaluate_treatment_hypothesis <- function(fit, serie_data_28, hypothesis_mortality, hypothesis_discharged, treatment_column, model_subgroup,
                                          adjusted, model_check, cl = NULL, cores = parallel::detectCores()
                                          ) {
  counterfactual_took <- serie_data_28 %>% mutate({{ treatment_column }} := 1)
  counterfactual_no <- serie_data_28 %>% mutate({{treatment_column }} := 0)

  mortality_log_OR_samples <- compute_counterfactual_state_log_OR(fit, counterfactual_took, counterfactual_no, target_state = "Death", day = 28, cl = cl, cores = cores)
  discharged_log_OR_samples <- compute_counterfactual_state_log_OR(fit, counterfactual_took, counterfactual_no, target_state = "Discharged", day = 28, cl = cl, cores = cores)

  rbind(
    bayesian_hypothesis_res_from_draws(
      draws = mortality_log_OR_samples,
      model = paste0("HMM", model_subgroup),
      estimand = "log(OR)",
      hypothesis = hypothesis_mortality,
      adjusted = adjusted,
      model_check = model_check
    ),
    bayesian_hypothesis_res_from_draws(
      draws = discharged_log_OR_samples,
      model = paste0("HMM", model_subgroup),
      estimand = "log(OR)",
      hypothesis = hypothesis_discharged,
      adjusted = adjusted,
      model_check = model_check
    )
  )
}

evaluate_all_treatment_hypotheses <- function(fit, model_subgroup, adjusted, model_check, serie_data_28 = serie_data_28,
                                              cl = NULL, cores = parallel::detectCores(), cache = "auto") {
  res <- NULL
  if(cache == "auto") {
    all_params <- list(fit, serie_data_28)

    cache_dir <- here::here("local_temp_data", "hmm", "epred_cache")
    if(!dir.exists(cache_dir)) {
      dir.create(cache_dir)
    }
    cache_file <- paste0(cache_dir,"/evaluate_treatments_", digest::digest(all_params), ".rds")
    if(file.exists(cache_file)) {
      res <- readRDS(cache_file)
    }
  }

  if(is.null(res)) {
    res <- rbind(
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$hcq_reduces_death, hypotheses$hcq_increases_discharged, took_hcq, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores),
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$az_reduces_death, hypotheses$az_increases_discharged, took_az, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores),
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$favipiravir_reduces_death, hypotheses$favipiravir_increases_discharged, took_favipiravir, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores)
    )
    if(cache == "auto") {
      saveRDS(res, cache_file)
    }
  }

  res
}

print_my_fit_summary <- function(brmshmmfit) {
  #rstan::check_hmc_diagnostics(brmshmmfit$brmsfit$fit)
  summ <- summary(brmshmmfit$brmsfit)
  cat("\nRate coefficients:\n")
  print(summ$fixed)
  cat("\nPer rate group effects:\n")
  print(brms::ranef(brmshmmfit$brmsfit)$rate_group)
}
