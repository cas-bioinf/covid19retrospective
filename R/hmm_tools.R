compute_counterfactual_outcome_samples <- function(hmm_fit, numerator_serie_data, denominator_serie_data, day = 28,
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

  death_state_id <- which(hmm_fit$data$hidden_state_data$id == "Death")
  discharged_state_id <- which(hmm_fit$data$hidden_state_data$id == "Discharged")

  numerator_indices <- unique(as.integer(combined_data$serie_data$.serie)[combined_data$serie_data$is_numerator])
  denominator_indices <- unique(as.integer(combined_data$serie_data$.serie)[!combined_data$serie_data$is_numerator])
  if(!identical(sort(c(numerator_indices, denominator_indices)), 1:length(unique(combined_data$serie_data$.serie)))) {
    stop("Bad indexing")
  }


  state_probs <- posterior_epred_rect(hmm_fit, newdata = combined_data, method = posterior_epred_state_prob, cl = cl, cores = cores)

  samples_death_probs_numerator <-
    rowMeans(state_probs[death_state_id, day + 1, , numerator_indices])
  samples_death_probs_denominator <-
    rowMeans(state_probs[death_state_id, day + 1, , denominator_indices])

  log_death_odds_numerator <- log(samples_death_probs_numerator) - log1p(-samples_death_probs_numerator)
  log_death_odds_denominator <- log(samples_death_probs_denominator) - log1p(-samples_death_probs_denominator)

  samples_discharged_probs_numerator <-
    rowMeans(state_probs[discharged_state_id, day + 1, , numerator_indices])
  samples_discharged_probs_denominator <-
    rowMeans(state_probs[discharged_state_id, day + 1, , denominator_indices])

  log_hospital_odds_numerator <- log1p(-samples_discharged_probs_numerator - samples_death_probs_numerator) - log(samples_discharged_probs_numerator)
  log_hospital_odds_denominator <- log1p(-samples_discharged_probs_denominator - samples_death_probs_denominator) - log(samples_discharged_probs_denominator)


  mean_from_probs <- function(x) {
    if(length(x) != day + 1) {
      stop("Bad dimension")
    }
    sum(diff(x) * 1:day) / x[day + 1]
  }
  # samples_hospital_mean_numerator <-
  #   rowMeans(
  #      apply(state_probs[discharged_state_id, , , numerator_indices], MARGIN = c(2,3), FUN = mean_from_probs)
  #   )
  # samples_hospital_mean_denominator <-
  #   rowMeans(
  #     apply(state_probs[discharged_state_id, , , denominator_indices], MARGIN = c(2,3), FUN = mean_from_probs)
  #   )
  samples_hospital_numerator <-
      apply(state_probs[discharged_state_id, , , numerator_indices], MARGIN = c(2,3), FUN = mean_from_probs)
  samples_hospital_denominator <-
      apply(state_probs[discharged_state_id, , , denominator_indices], MARGIN = c(2,3), FUN = mean_from_probs)


  samples_duration_ratios <- rowMeans(
    log(samples_hospital_numerator) - log(samples_hospital_denominator)
  )

  list(death_log_OR = log_death_odds_numerator - log_death_odds_denominator,
       hospital_log_OR = log_hospital_odds_numerator - log_hospital_odds_denominator,
       hospital_duration_log_ratio = samples_duration_ratios
  )
}

evaluate_treatment_hypothesis <- function(fit, serie_data_28, hypothesis_mortality, hypothesis_hospital, treatment_column, model_subgroup,
                                          adjusted, model_check, cl = NULL, cores = parallel::detectCores()
                                          ) {
  counterfactual_took <- serie_data_28 %>% mutate({{ treatment_column }} := 1)
  counterfactual_no <- serie_data_28 %>% mutate({{treatment_column }} := 0)

  counterfactual_res <- compute_counterfactual_outcome_samples(fit, counterfactual_took, counterfactual_no, day = 28, cl = cl, cores = cores)

  rbind(
    bayesian_hypothesis_res_from_draws(
      draws = counterfactual_res$death_log_OR,
      model = paste0("HMM", model_subgroup),
      estimand = "log(OR)",
      hypothesis = hypothesis_mortality,
      adjusted = adjusted,
      model_check = model_check
    ),
    bayesian_hypothesis_res_from_draws(
      draws = counterfactual_res$hospital_log_OR,
      model = paste0("HMM", model_subgroup),
      estimand = "log(OR)",
      hypothesis = hypothesis_hospital,
      adjusted = adjusted,
      model_check = model_check
    ),
    bayesian_hypothesis_res_from_draws(
      draws = counterfactual_res$hospital_duration_log_ratio,
      model = paste0("HMM", model_subgroup),
      estimand = "log(duration_ratio)",
      hypothesis = hypothesis_hospital,
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
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$hcq_death, hypotheses$hcq_hospital, took_hcq, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores),
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$az_death, hypotheses$az_hospital, took_az, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores),
      evaluate_treatment_hypothesis(fit, serie_data_28, hypotheses$favipiravir_death, hypotheses$favipiravir_hospital, took_favipiravir, model_subgroup = model_subgroup, adjusted = adjusted, model_check = model_check, cl = cl, cores = cores)
    )
    if(cache == "auto") {
      saveRDS(res, cache_file)
    }
  }

  res
}

print_hmm_fit_summary <- function(brmshmmfit) {
  #rstan::check_hmc_diagnostics(brmshmmfit$brmsfit$fit)
  summ <- summary(brmshmmfit$brmsfit)
  cat("\nRate coefficients:\n")
  fixed_summ <- summ$fixed[,c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")]
  colnames(fixed_summ) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
  print(fixed_summ)
  all_ranef <- brms::ranef(brmshmmfit$brmsfit)
  if(!is.null(all_ranef$rate_group)) {
    cat("\nPer rate group effects:\n")
    print(all_ranef$rate_group)
  }
  if(!is.null(all_ranef$.rate_id)) {
    cat("\nPer rate effects:\n")
    print(all_ranef$.rate_id)
  }

  if(!is.null(summ$random$hospital_id)) {
    cat("\nBetween-site differences:\n")
    hospital_summ <- summ$random$hospital_id[grepl("^sd", rownames(summ$random$hospital_id)),c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")]
    colnames(hospital_summ) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    print(hospital_summ)
  }
}

print_hmm_hypothesis_res <- function(hypo_res) {
  hypo_res %>%
    mutate(hypothesis = if_else(estimand == "log(duration_ratio)", paste0(hypothesis, " (duration)"), hypothesis)) %>%
    select(hypothesis, point_estimate, ci_low, ci_high) %>%
    rename(`95% CI Low` = ci_low, `95% CI High` = ci_high) %>% print()
}
