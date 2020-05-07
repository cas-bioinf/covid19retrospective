make_brms_formula_hmm <- function(formula) {
  if(!is.null(brms:::lhs(formula))) {
    stop("Formula needs to be one-sided")
  }
  update.formula(formula, .dummy ~ .)
}

make_stancode_hmm <- function(formula, serie_data, rate_data, state_data, prior = NULL) {
  formula <- make_brms_formula_hmm(formula)
  data <- make_data_hmm(formula, serie_data, rate_data, state_data)
  brms::make_stancode(formula, family = rate_hmm_family,
                      data = data$brmsdata,
                      stanvars = rate_hmm_stanvars(data$standata))
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

make_data_hmm <- function(formula, serie_data, rate_data, state_data) {

  serie_data <- serie_data %>% mutate(.predictor_set = 1:n())
  rate_data <- rate_data %>% mutate(.rate_id = 1:n())
  brmsdata <- crossing(serie_data, rate_data) %>% mutate(.dummy = 0)



  N_states_hidden <- length(unique(c(rate_data$.from, rate_data$.to)))

  N_rates = nrow(rate_data)

  state_data <- state_data %>% arrange(.id)
  if(!identical(state_data$.id, 1:nrow(state_data))) {
    stop("State ids need to be consecutive and not duplicated")
  }
  #TODO allow for noisy states
  N_states_observed <- length(unique(state_data$.corresponding_obs))

  N_series <- max(serie_data$.serie)
  if(!identical(sort(unique(serie_data$.serie)), 1:N_series)) {
    stop("Patient IDs need to be consecutive and not duplicated")
  }

  N_time <- max(serie_data$.time)

  #TODO be smart about this, using allvars (probably pass something from make_standata)
  N_predictor_sets <- nrow(serie_data)
  predictor_sets <- 1:N_predictor_sets

  obs_states <- serie_data$.observed
  obs_states[is.na(serie_data$.observed)] <- 0

  rate_predictors <- array(NA_integer_, c(N_predictor_sets, N_rates))
  for(i in 1:nrow(brmsdata)) {
    rate_predictors[brmsdata$.predictor_set[i], brmsdata$.rate_id[i]] <- i
  }

  standata <- loo::nlist(
    N_states_hidden,
    N_states_observed,
    N_rates,
    rates_from = rate_data$.from,
    rates_to = rate_data$.to,
    corresponding_observation = state_data$.corresponding_obs,

    N_noisy_states = 0,
    noisy_states = numeric(0),
    N_other_observations = 1,
    noisy_states_other_obs = array(0, c(0, 1)),

    sensitivity_low_bound = 0.5,
    initial_states_prob = state_data$.initial_prob,

    N_observations = nrow(serie_data),
    N_series,
    N_time,
    N_predictor_sets,

    series = serie_data$.serie,
    times = serie_data$.time,
    obs_states,
    predictor_sets = serie_data$.predictor_set,
    rate_predictors
  )

  loo::nlist(standata, brmsdata)
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

  simplex[N_states_hidden] initial_states_prob;

  // Observations
  int<lower=1> N_observations;
  int<lower=1> N_series;
  int<lower=1> N_time;
  int<lower=1> N_predictor_sets;

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
    alpha[, 1] = initial_states_prob .* observation_probs[, obs_states_rect[p, 1]];

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

