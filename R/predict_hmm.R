posterior_epred.brmshmmfit <- function(fit, nsamples = NULL, to_corresponding_states = FALSE) {
  validate_brmshmmfit(fit)
  brms:::contains_samples(fit$brmsfit)

  pred_rawdata <- fit$data

  data_hmm <- make_data_hmm(pred_rawdata)

  epred_mu <- brms::posterior_epred(fit$brmsfit, newdata = data_hmm$brmsdata, nsamples = nsamples)
  transition_matrices <- compute_all_transition_matrices(data_hmm, epred_mu)

  if(is.null(nsamples)) {
    nsamples <- dim(epred_mu)[1]
  }

  N_states_hidden <- data_hmm$standata$N_states_hidden

  result_rect <- array(NA_integer_, c(max(pred_rawdata$serie_data$.time), max(pred_rawdata$serie_data$.serie), nsamples))
  for(s in unique(pred_rawdata$serie_data$.serie)) {

    state <- pred_rawdata$initial_states[s]
    result_rect[1, s, ] <- state

    max_time <- fit$data$serie_data %>% filter(.serie == s) %>% pull(.time) %>% max()
    for(i in 1:nsamples) {
      for(t in 2:max_time) {
        ps <- data_hmm$standata$predictor_sets_rect[s, t]
        state <- sample(1:N_states_hidden, size = 1, prob = transition_matrices[state, , ps, i])
        result_rect[t, s, i] <- state
      }
    }
  }

  result <- array(NA_integer_, c(nsamples, nrow(pred_rawdata$serie_data)))

  for(o in 1:nrow(pred_rawdata$serie_data)) {
    result[, o] <- result_rect[pred_rawdata$serie_data$.time[o], pred_rawdata$serie_data$.serie[o], ]
  }

  if(to_corresponding_states) {
    data_hmm$standata$corresponding_observation[result]
  } else {
    result
  }
}

posterior_predict.brmshmmfit <- function(fit, nsamples = NULL) {
  validate_brmshmmfit(fit)
  epred <- posterior_epred.brmshmmfit(fit, nsamples = nsamples, to_corresponding_states = FALSE)

  if(is.null(nsamples)) {
    nsamples = dim(epred)[1]
  }

  observation_probs_samples <- rstan::extract(fit$brmsfit$fit, "observation_probs")$observation_probs

  N_hiden_states = dim(observation_probs_samples)[2]
  N_observed_states = dim(observation_probs_samples)[3]

  result <- array(NA_integer_, dim(epred))
  for(i in 1:nsamples) {
    observation_probs <- observation_probs_samples[i, ,]
    one_sample <- array(NA_integer_, dim(epred)[2])
    for(h in 1:N_hiden_states) {
      indices <- epred[i,] == h
      one_sample[indices] <- sample(1:N_observed_states, size = sum(indices), replace = TRUE,
                                    prob = observation_probs_samples[i, h, ])
    }
    result[i,] <- one_sample
  }

  result
}

prediction_to_wide_format <- function(data, prediction) {
  validate_brmshmmdata(data)
  nsamples <- dim(prediction)[1]

  samples_df <- as.data.frame(t(prediction))
  names(samples_df) <- paste0("__S", 1:nsamples)

  data$serie_data %>% cbind(samples_df) %>%
    pivot_longer(matches("^__S[0-9]*$"), names_prefix = "__S", names_to = ".sample", values_to = "prediction") %>%
    mutate(sample = as.integer(sample))

}

compute_all_transition_matrices <- function(data_hmm, epred_mu) {
  nsamples <- dim(epred_mu)[1]
  epred_mu_exp <- exp(epred_mu)
  standata <- data_hmm$standata
  res <- array(NA_real_, c(standata$N_states_hidden, standata$N_states_hidden, standata$N_predictor_sets, nsamples))
  for(i in 1:nsamples) {
    for(ps in 1:standata$N_predictor_sets) {
      rate_matrix <- matrix(0, standata$N_states_hidden, standata$N_states_hidden)

      rate_values <- epred_mu_exp[i, standata$rate_predictors[ps,]]
      for(r in 1:standata$N_rates) {
        rate_matrix[standata$rates_from[r], standata$rates_to[r]] = rate_values[r]
        rate_matrix[standata$rates_from[r],standata$rates_from[r]] =
          rate_matrix[standata$rates_from[r],standata$rates_from[r]] + rate_values[r]
      }

      res[,, ps, i] <-  expm::expm(rate_matrix)
    }
  }
  res
}
