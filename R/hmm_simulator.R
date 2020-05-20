hmm_simulator <- function(N_series, N_time, N_mid_states, use_noisy_states = FALSE, N_treatments = 1) {

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
  rates_group <- character(N_rates)

  next_rate <- 1
  rates_death <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_death] <- s_worsening_hid
  rates_to[rates_death] <- s_death_hid
  rates_group[rates_death] <- "death"
  next_rate <- next_rate + length(rates_death)

  rates_to_improving <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_to_improving] <- s_worsening_hid
  rates_to[rates_to_improving] <- s_improving_hid
  rates_group[rates_to_improving] <- "to_improving"
  next_rate <- next_rate + length(rates_to_improving)

  rates_to_worsening <- next_rate:(next_rate + N_mid_states - 1)
  rates_from[rates_to_worsening] <- s_improving_hid
  rates_to[rates_to_worsening] <- s_worsening_hid
  rates_group[rates_to_worsening] <- "to_worsening"
  next_rate <- next_rate + length(rates_to_worsening)

  rates_worsen_one <- next_rate:(next_rate + N_mid_states - 2)
  rates_from[rates_worsen_one] <- s_worsening_hid[1:(length(s_worsening_hid) - 1)]
  rates_to[rates_worsen_one] <- s_worsening_hid[2:length(s_worsening_hid)]
  rates_group[rates_worsen_one] <- "worsen_one"
  next_rate <- next_rate + length(rates_worsen_one)

  rates_improve_one <- next_rate:(next_rate + N_mid_states - 2)
  rates_from[rates_improve_one] <- s_improving_hid[2:length(s_improving_hid)]
  rates_to[rates_improve_one] <- s_improving_hid[1:(length(s_improving_hid) - 1)]
  rates_group[rates_improve_one] <- "improve_one"
  next_rate <- next_rate + length(rates_improve_one)

  rate_discharge <- next_rate
  rates_from[rate_discharge] <- s_improving_hid[1]
  rates_to[rate_discharge] <- s_discharged_hid
  rates_group[rate_discharge] <- "improve_one"

  rates_group <- factor(rates_group)

  # Model coefficients
  intercept <- rnorm(1, -3, 1)

  group_effects <- numeric(length(levels(rates_group)))
  names(group_effects) <- levels(rates_group)

  group_effects["death"] <- rnorm(1, -2.3, 1)
  group_effects["to_improving"] <- 0 #Using this as the reference category
  group_effects["to_worsening"] <- rnorm(1, 0, 1)
  group_effects["worsen_one"] <- rnorm(1, 0.7, 1)
  group_effects["improve_one"] <- rnorm(1, 0.7, 1)

  rate_id_sd_intercept <- abs(rnorm(1,0,0.5))
  rate_id_intercept <- rnorm(N_rates, 0, rate_id_sd_intercept)

  if(N_treatments > 1) {
    stop("More than one treatment currently not supported")
  }
  group_sd_treated <- abs(rnorm(1,0,1))
  group_treated <- rnorm(length(levels(rates_group)), 0, group_sd_treated)


  # Predictors
  N_predictor_sets <- N_treatments + 1
  N <- N_predictor_sets * N_rates
  rate_predictors <- array(1:(N_rates * N_predictor_sets), c(N_predictor_sets,N_rates))

  mu <- array(NA_real_, N)
  for(ps in 1:N_predictor_sets) {
    mu[rate_predictors[ps,]] <- intercept + group_effects[levels(rates_group)[rates_group]] + rate_id_intercept
    if(ps == 2) {
      mu[rate_predictors[ps,]] <- mu[rate_predictors[ps,]] + group_treated[rates_group]
    }
  }

  # Transition matrices from predictors
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

  sensitivity <- runif(N_noisy_states, min = sensitivity_low_bound, max = 1)


  hidden_state_data <- data.frame(id = 1:N_states_hidden,
                          corresponding_obs = corresponding_observation)


  observed_state_data <- data.frame(id = 1:N_states_observed,
                          is_noisy = (1:N_states_observed) %in% noisy_states)
  if(N_noisy_states > 0) {
    for(other_obs in 1:N_other_observations) {
      observed_state_data[[paste0("other_obs_", other_obs)]] <- noisy_states_other_obs[noisy_state_id[observed_state_data$id], other_obs]
    }
  }


  # The actual Markov chain
  N_observations <- N_series * N_time - floor(N_series / 2)
  series <- array(NA_integer_, N_observations)
  times <- array(NA_integer_, N_observations)
  obs_states <- array(NA_integer_, N_observations)
  predictor_sets <- array(NA_integer_, N_observations)
  treated <- logical(N_observations)

  initial_states <- sample(s_worsening_hid, N_series, replace = TRUE)


  true_base_states <- array(NA_integer_, N_observations)
  true_improving <- array(FALSE, N_observations)
  next_observation <- 1
  for(p in 1:N_series) {
    if(p %% 2 == 1) {
      max_time <- N_time
    } else {
      max_time <- N_time - 1
    }
    state <- initial_states[p]
    for(t in 1:max_time) {
      series[next_observation] <- p
      times[next_observation] <- t

      cor_obs <- corresponding_observation[state]
      true_base_states[next_observation] <- cor_obs
      true_improving[next_observation] <- state %in% s_improving_hid

      if(N_treatments > 0) {
        if(p %% 2 == 1 && t > 4) {
          ps <- 1
          treated[next_observation] <- FALSE
        } else {
          ps <- 2
          treated[next_observation] <- TRUE
        }
      } else {
        ps <- 1
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

  model_formula <- ~ 0 + Intercept + rate_death + rate_to_worsening + rate_improve_one + rate_worsen_one + (1 || .rate_id)
  prior <- brms::set_prior("normal(-2.3,1)", "b", coef = "Intercept") +
    brms::set_prior("normal(-2.3,1)", "b", coef = "rate_death") +
    brms::set_prior("normal(0,1)", "b", coef = "rate_to_worsening") +
    brms::set_prior("normal(0.7,1)", "b", coef = "rate_improve_one") +
    brms::set_prior("normal(0.7,1)", "b", coef = "rate_worsen_one") +
    brms::set_prior("normal(0, 0.5)", "sd", group = ".rate_id")

  if(N_treatments > 0) {
    model_formula <- update.formula(model_formula, ~ . + (0 + treated || group))
    prior <- prior + brms::set_prior("normal(0, 1)", "sd", group = "group")
  }
  list(
    observed = brmshmmdata(
        formula = model_formula,
        prior = prior,
        serie_data = data.frame(
          .serie = as.integer(series),
          .time = as.integer(times),
          .observed = as.integer(obs_states),
          treated = as.numeric(treated)),
        rate_data = data.frame(
          .from = as.integer(rates_from),
          .to = as.integer(rates_to),
          group = rates_group,
          rate_death = as.numeric(rates_group == "death"),
          rate_to_worsening = as.numeric(rates_group == "to_worsening"),
          rate_improve_one = as.numeric(rates_group == "improve_one"),
          rate_worsen_one = as.numeric(rates_group == "worsen_one")
        ),
        hidden_state_data = hidden_state_data,
        initial_states = initial_states,
        sensitivity_low_bound = sensitivity_low_bound,
        observed_state_data = observed_state_data
    ),
    true = loo::nlist(
      # BRMS model
      b = c(intercept, group_effects["death"], group_effects["to_worsening"], group_effects["improve_one"], group_effects["worsen_one"]),
      sd_1 = array(rate_id_sd_intercept, 1),
      sd_2 = array(group_sd_treated, 1),
      r_1_1 = rate_id_intercept,
      r_2_1 = group_treated,

      # HMM model
      sensitivity,
      other_observations_probs,

      # Derived quantitites and helpers
      mu,
      transition_matrices,
      true_base_states,
      true_improving
    )
  )
}
