brmhmm <- function(brmshmmdata) {
  d <- validate_brmshmmdata(brmshmmdata)

  prepdata <- make_data_hmm(d)

  bfit <- brms::brm(
    formula = d$formula,
    family = rate_hmm_family,
    data = prepdata$brmsdata,
    prior = d$prior,
    stanvars = rate_hmm_stanvars(prepdata$standata)
    )

  structure(list(
    brmsfit = bfit,
    data = brmshmmdata
  ), class = "brmshmmfit")
}


validate_brmshmmfit <- function(fit) {
  if(!inherits(fit, "brmshmmfit")) {
    stop("Not of class brmshmmfit")
  }
  fit
}

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

  result
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
