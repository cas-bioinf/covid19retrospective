make_data_hmm <- function(formula_processed,
                          serie_data, rate_data,
                          hidden_state_data, initial_states,
                          observed_state_data = NULL, sensitivity_low_bound = 0.5) {

  rate_data <- rate_data %>% mutate(.rate_id = factor(1:n()))



  N_states_hidden <- length(unique(c(rate_data$.from, rate_data$.to)))

  N_rates = nrow(rate_data)

  hidden_state_data <- hidden_state_data %>% arrange(id)
  if(!identical(hidden_state_data$id, 1:N_states_hidden)) {
    stop("Hiden state ids need to be consecutive and not duplicated and span the same range as rate_data")
  }

  all_corresponding_states <- hidden_state_data %>% pull(corresponding_obs) %>%
    unique() %>% sort()
  if(is.null(observed_state_data)) {

    if(!identical(all_corresponding_states, 1:length(all_corresponding_states))) {
      stop("If observed_state_data is NULL, corresponding_obs cannot have gaps")
    }

    observed_state_data <- data.frame(id = all_corresponding_states, is_noisy = FALSE)
  }

  N_states_observed <- nrow(observed_state_data)
  observed_state_data <- observed_state_data %>% arrange(id)
  if(!identical(observed_state_data$id, 1:N_states_observed)) {
    stop("Observed state ids need to be consecutive.")
  }

  if(!all(all_corresponding_states %in% observed_state_data$id)) {
    stop("Some corresponding observations are not in observed_state_data")
  }





  N_noisy_states <- sum(observed_state_data$is_noisy)
  if(N_noisy_states > 0) {
    noisy_states <- observed_state_data %>% filter(is_noisy) %>% pull(id)

    N_other_observations <- observed_state_data %>% select(starts_with("other_obs_")) %>% length()
    if(N_other_observations < 1) {
      stop("Noisy states require other_obs_x columns")
    }

    noisy_states_other_obs <- observed_state_data %>%
      filter(is_noisy) %>%
      select(starts_with("other_obs_")) %>%
      as.matrix()
  } else {
    noisy_states = numeric(0)
    N_other_observations = 1
    noisy_states_other_obs = array(0, c(0, 1))
  }


  N_series <- max(serie_data$.serie)
  if(!identical(sort(unique(serie_data$.serie)), 1:N_series)) {
    stop("Patient IDs need to be consecutive and not duplicated")
  }

  N_time <- max(serie_data$.time)

  #TODO be smart about this, using allvars (probably pass something from make_standata)
  all_vars_needed <- brms:::all_vars(formula_processed)

  serie_data_vars <- intersect(all_vars_needed, names(serie_data))

  serie_data_distinct <- serie_data %>% select(one_of(serie_data_vars)) %>%
    dplyr::distinct() %>%
    mutate(.predictor_set = 1:n())

  N_predictor_sets <- nrow(serie_data_distinct)

  serie_data_raw <- serie_data
  serie_data <- serie_data %>% dplyr::left_join(serie_data_distinct, by = serie_data_vars)
  if(any(is.na(serie_data$.predictor_set)) || nrow(serie_data_raw) != nrow(serie_data)) {
    stop("Failed join")
  }

  brmsdata_all <- crossing(serie_data_distinct, rate_data)

  #TODO avoid repetitions (would only apply if some rates are identical)
  # brmsdata_distinct <- brmsdata_all %>% select(one_of(all_vars_needed)) %>%
  #   dplyr::distinct() %>%
  #   mutate(.predictor_id = 1:n())
  #
  # brmsdata <- brmsdata_all %>% left_join(
  brmsdata <- brmsdata_all %>%
    arrange(.predictor_set, .rate_id) %>%
    mutate(.predictor_id = 1:n(), .dummy = rnorm(n()))

  obs_states <- serie_data$.observed
  obs_states[is.na(serie_data$.observed)] <- 0

  rate_predictors <- array(NA_integer_, c(N_predictor_sets, N_rates))
  for(i in 1:nrow(brmsdata)) {
    rate_predictors[brmsdata$.predictor_set[i], brmsdata$.rate_id[i]] <- brmsdata$.predictor_id[i]
  }

  standata <- loo::nlist(
    N_states_hidden,
    N_states_observed,
    N_rates,
    rates_from = rate_data$.from,
    rates_to = rate_data$.to,
    corresponding_observation = hidden_state_data$corresponding_obs,

    N_noisy_states,
    noisy_states,
    N_other_observations,
    noisy_states_other_obs,

    sensitivity_low_bound,

    N_observations = nrow(serie_data),
    N_series,
    N_time,
    N_predictor_sets,

    initial_states,
    series = serie_data$.serie,
    times = serie_data$.time,
    obs_states,
    predictor_sets = serie_data$.predictor_set,
    rate_predictors
  )

  loo::nlist(standata, brmsdata)
}

make_standata_hmm <- function(formula, serie_data, rate_data,
                              hidden_state_data, initial_states,
                              observed_state_data = NULL,
                              sensitivity_low_bound = 0.5) {
  formula <- make_brms_formula_hmm(formula)
  data <- make_data_hmm(formula, serie_data = serie_data, rate_data = rate_data,
                        hidden_state_data = hidden_state_data,
                        initial_states = initial_states,
                        observed_state_data = observed_state_data,
                        sensitivity_low_bound = sensitivity_low_bound)
  c(brms::make_standata(formula, data = data$brmsdata), data$standata)
}
