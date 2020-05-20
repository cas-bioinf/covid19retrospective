posterior_epred.brmshmmfit <- function(fit, nsamples = NULL) {
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

  result_rect <- array(NA_integer_, c(max(pred_rawdata$serie_data$.time), max(as.integer(pred_rawdata$serie_data$.serie)), nsamples))
  for(s in unique(as.integer(pred_rawdata$serie_data$.serie))) {

    result_rect[1, s, ] <- as.integer(pred_rawdata$initial_states[s])

    max_time <- pred_rawdata$serie_data %>% filter(as.integer(.serie) == s) %>% pull(.time) %>% max()

    for(i in 1:nsamples) {
      state <- pred_rawdata$initial_states[s]
      for(t in 2:max_time) {
        ps <- data_hmm$standata$predictor_sets_rect[s, t]
        if(ps == 0) {
          stop(paste0("Element with no predictor at serie ",s,", time ", t))
        }
        #cat(state, ":", ps, ",", i, ":", transition_matrices[state, , ps, i], "\n")
        state <- sample.int(N_states_hidden, size = 1, prob = transition_matrices[state, , ps, i])
        if(is.na(state)) {
          stop("NA state")
        }
        result_rect[t, s, i] <- state
      }
    }
  }

  result <- array(NA_integer_, c(nsamples, nrow(pred_rawdata$serie_data)))

  for(o in 1:nrow(pred_rawdata$serie_data)) {
    result_for_row <- result_rect[pred_rawdata$serie_data$.time[o], as.integer(pred_rawdata$serie_data$.serie[o]), ]
    if(any(is.na(result_for_row))) {
      stop(paste0("NA for row ", o))
    }
    result[, o] <- result_for_row
  }

  if(any(is.na(result))) {
    stop("NA in translation from rect")
  }
  result
}

hidden_to_corresponding_observed <- function(data_processed, x) {
  if(!is.integer(x) || any(x < 1 | x > data_processed$standata$N_states_hidden)) {
    stop("Invalid state data")
  }
  res_corresponding <- data_processed$standata$corresponding_observation[x]
  dim(res_corresponding) <- dim(x)

  res_corresponding

}

posterior_predict.brmshmmfit <- function(fit, nsamples = NULL) {
  validate_brmshmmfit(fit)
  epred <- posterior_epred.brmshmmfit(fit, nsamples = nsamples)

  if(any(is.na(epred))) {
    stop("NA in epred")
  }

  if(is.null(nsamples)) {
    nsamples = dim(epred)[1]
  }

  observation_probs_samples <- rstan::extract(fit$brmsfit$fit, "observation_probs")$observation_probs

  N_hidden_states = dim(observation_probs_samples)[2]
  N_observed_states = dim(observation_probs_samples)[3]

  result <- array(NA_integer_, dim(epred))
  for(i in 1:nsamples) {
    observation_probs <- observation_probs_samples[i, ,]
    one_sample <- array(NA_integer_, dim(epred)[2])
    for(h in 1:N_hidden_states) {
      indices <- epred[i,] == h
      one_sample[indices] <- sample.int(N_observed_states, size = sum(indices), replace = TRUE,
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
    pivot_longer(matches("^__S[0-9]*$"), names_prefix = "__S", names_to = ".sample", values_to = ".predicted") %>%
    mutate(.sample = as.integer(.sample))

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
          rate_matrix[standata$rates_from[r],standata$rates_from[r]] - rate_values[r]
      }

      res[,, ps, i] <-  expm::expm(rate_matrix)
    }
  }
  res
}
