hmm_simulator <- function(N_series, N_time, N_mid_states, use_noisy_states = FALSE) {

  # States and transitions
  if(N_mid_states < 1) {
    stop("Need at least one mid state")
  }
  N_states_hidden <- 2 + N_mid_states * 2
  N_states_observed <- 2 + N_mid_states

  s_discharged_hid <- 1
  s_improving_hid <- 2:(N_mid_states + 1)
  s_worsening_hid <- (N_mid_states + 2):(N_mid_states * 2 + 1)
  s_death_hid <- N_states_hidden

  s_discharged_obs <- 1
  s_mid_obs <- 2:(N_mid_states + 1)
  s_death_obs <- N_states_observed

  N_rates <- N_mid_states + #worsening to death
    2 * N_mid_states + #worsening <-> improving
    2 * (N_mid_states - 1) + #worsen/improve by one
    1 #discharge

  rates_from <- array(0, N_rates)
  rates_to <- array(0, N_rates)
  rates_prior_mean <- numeric(N_rates)

  next_rate <- 1
  rates_death <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_death] <- s_worsening_hid
  rates_to[rates_death] <- s_death_hid
  rates_prior_mean[rates_death] <- 0.01
  next_rate <- next_rate + length(rates_death)

  rates_to_improving <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_to_improving] <- s_worsening_hid
  rates_to[rates_to_improving] <- s_improving_hid
  rates_prior_mean[rates_to_improving] <- 0.1
  next_rate <- next_rate + length(rates_to_improving)

  rates_to_worsening <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_to_worsening] <- s_improving_hid
  rates_to[rates_to_worsening] <- s_worsening_hid
  rates_prior_mean[rates_to_worsening] <- 0.1
  next_rate <- next_rate + length(rates_to_worsening)

  rates_worsen_one <- next_rate:(next_rate + N_mid_states - 2)
  rates_from[rates_worsen_one] <- s_worsening_hid[1:(length(s_worsening_hid) - 1)]
  rates_to[rates_worsen_one] <- s_worsening_hid[2:length(s_worsening_hid)]
  rates_prior_mean[rates_worsen_one] <- 0.2
  next_rate <- next_rate + length(rates_worsen_one)

  rates_improve_one <- next_rate:(next_rate + N_mid_states - 2)
  rates_from[rates_improve_one] <- s_improving_hid[2:length(s_improving_hid)]
  rates_to[rates_improve_one] <- s_improving_hid[1:(length(s_improving_hid) - 1)]
  rates_prior_mean[rates_improve_one] <- 0.2
  next_rate <- next_rate + length(rates_improve_one)

  rate_discharge <- next_rate
  rates_from[rate_discharge] <- s_improving_hid[1]
  rates_to[rate_discharge] <- s_discharged_hid
  rates_prior_mean[rate_discharge] <- 0.1


  # Predictors
  N_predictor_sets <- 2
  N <- N_predictor_sets * N_rates
  rate_predictors <- array(1:(N_rates * N_predictor_sets), c(N_predictor_sets,N_rates))

  intercept_prior_logmean <- array(NA_real_, N)
  for(ps in 1:N_predictor_sets) {
    intercept_prior_logmean[rate_predictors[ps,]] <- log(rates_prior_mean)
  }
  intercept_prior_sd <- array(1, N);
  mu <- rnorm(N, intercept_prior_logmean, intercept_prior_sd)

  transition_matrices <- array(NA_real_, c(N_predictor_sets,N_states_hidden, ncol = N_states_hidden))

  for(ps in 1:N_predictor_sets) {
    rate_matrix <- matrix(0, nrow = N_states_hidden, ncol = N_states_hidden)
    for(r in 1:N_rates) {
      rate_value <- exp( mu[ rate_predictors[ps, r] ] )
      rate_matrix[rates_from[r], rates_to[r]] <- rate_value
      rate_matrix[rates_from[r], rates_from[r]] <- rate_matrix[rates_from[r], rates_from[r]] - rate_value
    }
    transition_matrices[ps, ,] <- expm::expm(rate_matrix)
  }


  # Observation model
  sensitivity_low_bound <- 0.8
  sensitivity <- runif(N_mid_states, min = sensitivity_low_bound, max = 1)
  corresponding_observation <- array(NA_integer_, N_states_hidden)
  corresponding_observation[s_discharged_hid] <- s_discharged_obs
  corresponding_observation[s_worsening_hid] <- s_mid_obs
  corresponding_observation[s_improving_hid] <- s_mid_obs
  corresponding_observation[s_death_hid] <-  s_death_obs
  if(use_noisy_states) {
    N_noisy_states <- N_mid_states
    noisy_states <- s_mid_obs
    noisy_state_id <- array(NA_integer_, N_states_observed)
    noisy_state_id[noisy_states] <- 1:N_noisy_states

    if(N_noisy_states == 1) {
      stop("Need at least 2 mid states for noisy states")
    } else if(N_noisy_states == 2) {
      N_other_observations <- 1
      noisy_states_other_obs <- array(NA_integer_, c(N_noisy_states, N_other_observations))
      noisy_states_other_obs[1, ] <- noisy_states[2]
      noisy_states_other_obs[2, ] <- noisy_states[1]
    } else {
      N_other_observations <- 2
      noisy_states_other_obs <- array(NA_integer_, c(N_noisy_states, N_other_observations))
      noisy_states_other_obs[1, ] <- c(noisy_states[1] + 1,noisy_states[1] + 2)
      noisy_states_other_obs[N_noisy_states, ] <- c(noisy_states[N_noisy_states] - 2,noisy_states[N_noisy_states] - 1)
      for(ns in 2:(N_noisy_states - 1)) {
        noisy_states_other_obs[ns, ] <- c(noisy_states[ns] - 1, noisy_states[ns] + 1)
      }
    }

    other_observations_probs <- array(NA_real_, c(N_noisy_states, N_other_observations))
    for(ns in 1:N_noisy_states) {
      other_observations_probs[ns, ] <- MCMCpack::rdirichlet(1, rep(1, N_other_observations))
    }
  } else {
    N_noisy_states <- 0
    noisy_states <- integer(0)
    noisy_state_id <- integer(0)

    N_other_observations <- 1
    noisy_states_other_obs <- array(as.integer(0), c(0, 1))

    other_observations_probs <- array(0, c(0, 1))
  }



  # The actual Markov chain
  N_observations <- N_series * N_time - floor(N_series / 2)
  series <- array(NA_integer_, N_observations)
  times <- array(NA_integer_, N_observations)
  obs_states <- array(NA_integer_, N_observations)
  predictor_sets <- array(NA_integer_, N_observations)

  initial_states_prob <- numeric(N_states_hidden)
  initial_states_prob[s_worsening_hid] <- 1 / N_mid_states

  true_base_states <- array(NA_integer_, N_observations)
  true_improving <- array(FALSE, N_observations)
  next_observation <- 1
  for(p in 1:N_series) {
    if(p %% 2 == 1) {
      max_time <- N_time
    } else {
      max_time <- N_time - 1
    }
    state <- sample(1:N_states_hidden, 1, prob = initial_states_prob)
    for(t in 1:max_time) {
      series[next_observation] <- p
      times[next_observation] <- t

      cor_obs <- corresponding_observation[state]
      true_base_states[next_observation] <- cor_obs
      true_improving[next_observation] <- state %in% s_improving_hid

      if(p %% 2 == 1 && t > 4) {
        ps <- 1
      } else {
        ps <- 2
      }
      predictor_sets[next_observation] <- ps

      if(cor_obs %in% noisy_states) {
        noisy_id <- noisy_state_id[cor_obs]
        if(runif(1) < sensitivity[noisy_id]) {
          obs_states[next_observation] <- cor_obs
        } else if(N_other_observations == 1) {
          obs_states[next_observation] <- noisy_states_other_obs[noisy_id, 1]
        } else {
          obs_states[next_observation] <- sample(noisy_states_other_obs[noisy_id, ], size = 1,
                                       prob = other_observations_probs[noisy_id, ])
        }
      } else {
        obs_states[next_observation] <- cor_obs
      }

      #new state
      state <- sample(1:N_states_hidden, size = 1, prob = transition_matrices[ps, state, ])
      next_observation <- next_observation + 1
    }
  }

  list(
    observed = list(
      N = N,
      intercept_prior_logmean = intercept_prior_logmean,
      intercept_prior_sd = intercept_prior_sd,

      N_states_hidden = N_states_hidden,
      N_states_observed = N_states_observed,

      N_rates = N_rates,
      rates_from = rates_from,
      rates_to = rates_to,

      corresponding_observation = corresponding_observation,

      N_noisy_states = N_noisy_states,
      noisy_states = noisy_states,
      N_other_observations = N_other_observations,
      noisy_states_other_obs = noisy_states_other_obs,

      sensitivity_low_bound = sensitivity_low_bound,

      initial_states_prob = initial_states_prob,

      N_observations = N_observations,
      N_series = N_series,
      N_time = N_time,
      N_predictor_sets = N_predictor_sets,

      series = series,
      times = times,
      obs_states = obs_states,

      predictor_sets = predictor_sets,
      rate_predictors = rate_predictors
    ),
    observed_structured = list(
      formula = ~ .rate_id * treated,
      serie_data = data.frame(
        .serie = series,
        .time = times,
        .observed = obs_states,
        treated = factor(predictor_sets)),
      rate_data = data.frame(
        .from = rates_from,
        .to = rates_to
      ),
      prior = brms::set_prior("normal(0,1)") + brms::set_prior("normal(0,1)", "Intercept"),
      state_data = data.frame(.id = 1:N_states_hidden,
                              .corresponding_obs = corresponding_observation,
                              .initial_prob = initial_states_prob)
    ),
    true = list(
      mu = mu,
      sensitivity = sensitivity,
      other_observations_probs = other_observations_probs,
      transition_matrices = transition_matrices,
      true_base_states = true_base_states,
      true_improving = true_improving
    )
  )
}
