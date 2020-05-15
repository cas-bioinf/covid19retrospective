make_brms_formula_hmm <- function(formula) {
  if(!is.null(brms:::lhs(formula))) {
    stop("Formula needs to be one-sided")
  }
  update.formula(formula, .dummy ~ .)
}

make_stancode_hmm <- function(formula,
                              serie_data, rate_data,
                              hidden_state_data, initial_states,
                              observed_state_data = NULL, prior = NULL) {
  formula <- make_brms_formula_hmm(formula)
  data <- make_data_hmm(formula, serie_data, rate_data,
                        hidden_state_data, initial_states, observed_state_data)
  brms::make_stancode(formula, family = rate_hmm_family,
                      data = data$brmsdata,
                      stanvars = rate_hmm_stanvars(data$standata), prior = prior)
}

rate_hmm_family <- structure(list(family = "rate_hmm", link = "identity", dpars = "mu",
                           lb = NA, ub = NA, type = "real", vars = NULL, specials = NULL,
                           log_lik = NULL, posterior_predict = NULL, posterior_epred = NULL),
                      class = c("rate_hmm", "brmsfamily"))

.family_rate_hmm <- function() {

}

rate_hmm_stanvars <- function(standata) {

  brms::stanvar(scode = rate_hmm_functions_code, block = "functions") +
  brms::stanvar(x = standata, scode = rate_hmm_data_code, block = "data") +
  brms::stanvar(scode = rate_hmm_tdata_code, block = "tdata") +
  brms::stanvar(scode = rate_hmm_parameters_code, block = "parameters")
}

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

rate_hmm_functions_code <- "
  // Compute a single transition matrix
  matrix compute_transition_matrix(
    int N_states, int N_rates, int[] rates_from, int[] rates_to, vector rates
  ) {
      matrix[N_states, N_states] rate_matrix = rep_matrix(0, N_states, N_states);
      vector[N_states] outgoing_rate_sum = rep_vector(0, N_states);
      for(r in 1:N_rates) {
        rate_matrix[rates_from[r], rates_to[r]] = rates[r];
        outgoing_rate_sum[rates_from[r]] += rates[r];
      }
      for(s in 1:N_states) {
        rate_matrix[s,s] = -outgoing_rate_sum[s];
      }
      return matrix_exp(rate_matrix);
  }
"

rate_hmm_data_code <- '
  // HMM data
  int<lower=1> N_states_hidden;
  int<lower=1> N_states_observed;

  int<lower=1> N_rates;
  int<lower=1, upper=N_states_hidden> rates_from[N_rates];
  int<lower=1, upper=N_states_hidden> rates_to[N_rates];

  int<lower=1, upper=N_states_observed> corresponding_observation[N_states_hidden];

  int<lower=0> N_noisy_states;
  int<lower=1, upper=N_states_observed> noisy_states[N_noisy_states];
  int<lower=1> N_other_observations;
  int<lower=1, upper=N_states_observed> noisy_states_other_obs[N_noisy_states, N_other_observations];

  real<lower=0, upper=1> sensitivity_low_bound;


  // Observations
  int<lower=1> N_observations;
  int<lower=1> N_series;
  int<lower=1> N_time;
  int<lower=1> N_predictor_sets;

  int<lower=1, upper=N_states_hidden> initial_states[N_series];
  int<lower=1, upper=N_series> series[N_observations];
  int<lower=1, upper=N_time> times[N_observations];
  //0 for unobserved states
  int<lower=0, upper=N_states_observed> obs_states[N_observations];

  int<lower=1, upper=N_predictor_sets> predictor_sets[N_observations];
  int<lower=1, upper=N> rate_predictors[N_predictor_sets, N_rates];
'

rate_hmm_tdata_code <- '
  int<lower=0,upper=1> is_state_noisy[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=0,upper=N_noisy_states> noisy_state_id[N_states_observed] = rep_array(0, N_states_observed);

  //Rectangule observations and predictors. 0 for missing data
  int<lower=0, upper=N_states_observed> obs_states_rect[N_series, N_time] = rep_array(0, N_series, N_time);
  int<lower=0, upper=N_predictor_sets> predictor_sets_rect[N_series, N_time] = rep_array(0, N_series, N_time);
  int<lower=0, upper=N_time> max_time[N_series] = rep_array(0,N_series);

  for(s_index in 1:N_noisy_states) {
    is_state_noisy[noisy_states[s_index]] = 1;
    noisy_state_id[noisy_states[s_index]] = s_index;
  }

  for(o in 1:N_observations) {
    obs_states_rect[series[o], times[o]] = obs_states[o];
    predictor_sets_rect[series[o], times[o]] = predictor_sets[o];
    if(obs_states[o] != 0) {
      max_time[series[o]] = max(max_time[series[o]], times[o]);
    }
  }

  //Check data validity
  for(p in 1:N_series) {
    if(max_time[p] == 0) {
      reject("serie ", p, " has no osbservations");
    }
    for(t in 1:(max_time[p] - 1)) {
      if(predictor_sets_rect[p, t] == 0) {
        reject("serie ", p, " has missing predictors for time ", t);
      }
    }
  }
'

rate_hmm_parameters_code <- "
  vector<lower=sensitivity_low_bound, upper=1>[N_noisy_states] sensitivity;
  simplex[N_other_observations] other_observations_probs[N_noisy_states];
"

stan_llh.rate_hmm <- function( ... ) {
"
  matrix[N_states_hidden, N_states_observed] observation_probs = rep_matrix(0, N_states_hidden, N_states_observed);
  matrix[N_states_hidden, N_states_hidden] transition_matrices_t[N_predictor_sets];
  vector[N] exp_mu = exp(mu);

  //Compute all transition matrices
  for(ps in 1:N_predictor_sets) {
    vector[N_rates] rate_values = exp_mu[rate_predictors[ps]];
    transition_matrices_t[ps] =
      transpose(compute_transition_matrix(N_states_hidden,  N_rates, rates_from, rates_to, rate_values));
  }

  for(s_hidden in 1:N_states_hidden) {
    int corresponding_obs = corresponding_observation[s_hidden];
    if(is_state_noisy[corresponding_obs]) {
      int noisy_id = noisy_state_id[corresponding_obs];
      observation_probs[s_hidden, corresponding_obs] = sensitivity[noisy_id];
      for(other_index in 1:N_other_observations) {
        int s_observed = noisy_states_other_obs[noisy_id, other_index];
        observation_probs[s_hidden, s_observed] =
          (1 - sensitivity[noisy_id]) * other_observations_probs[noisy_id, other_index];
      }
    } else {
      observation_probs[s_hidden, corresponding_obs] = 1;
    }
  }

  for(p in 1:N_series) {
    matrix[N_states_hidden, max_time[p]] alpha;
    vector[max_time[p]] alpha_log_norms;
    alpha[, 1] = rep_vector(0, N_states_hidden);
    alpha[initial_states[p], 1] = observation_probs[initial_states[p], obs_states_rect[p, 1]];

    {
      real norm = max(alpha[, 1]);
      alpha[, 1] /= norm;
      alpha_log_norms[1] = log(norm);
    }
    for(t in 2:max_time[p]) {
      real col_norm;
      int ps = predictor_sets_rect[p, t - 1];
      alpha[, t] = observation_probs[, obs_states_rect[p, t]] .* (transition_matrices_t[ps] * alpha[, t - 1]);
      col_norm = max(alpha[, t]);
      alpha[, t] /= col_norm;
      alpha_log_norms[t] = log(col_norm) + alpha_log_norms[t - 1];
    }
    target += log(sum(alpha[,max_time[p]])) + alpha_log_norms[max_time[p]];
  }

"
}

