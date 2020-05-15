make_brms_formula_hmm <- function(formula) {
  if(!is.null(brms:::lhs(formula))) {
    stop("Formula needs to be one-sided")
  }

  brms::brmsformula(update.formula(formula, .dummy ~ .), family = rate_hmm_family)
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
                           ybounds = c(-Inf, Inf),
                           log_lik = NULL, posterior_predict = NULL, posterior_epred = NULL),
                      class = c("rate_hmm", "brmsfamily"))

.family_rate_hmm <- function() {

}

rate_hmm_stanvars <- function(standata) {

  brms::stanvar(scode = rate_hmm_functions_code, block = "functions") +
  rate_hmm_stanvars_data(standata) +
  brms::stanvar(scode = rate_hmm_tdata_code, block = "tdata") +
  brms::stanvar(scode = rate_hmm_parameters_code, block = "parameters") +
  brms::stanvar(scode = rate_hmm_genquant_code, block = "genquant")

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

  matrix compute_observation_matrix(
    int N_states_hidden, int N_states_observed, int[] corresponding_observation,
    int[] is_state_noisy, int[] noisy_state_id,
    vector sensitivity, int N_other_observations, int[,] noisy_states_other_obs,
    vector[] other_observations_probs
  ) {
    matrix[N_states_hidden, N_states_observed] observation_probs = rep_matrix(0, N_states_hidden, N_states_observed);
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

    return observation_probs;
  }
"

rate_hmm_stanvars_data  <- function(standata) {
  brms::stanvar(x = standata$N_states_hidden, name = "N_states_hidden", scode = "  // HMM data\n  int<lower=1> N_states_hidden;", block = "data") +
  brms::stanvar(x = standata$N_states_observed, name = "N_states_observed", scode = "  int<lower=1> N_states_observed;", block = "data") +

  brms::stanvar(x = standata$N_rates, name = "N_rates", scode = "  int<lower=1> N_rates;", block = "data") +
  brms::stanvar(x = standata$rates_from, name = "rates_from", scode = "  int<lower=1, upper=N_states_hidden> rates_from[N_rates];", block = "data") +
  brms::stanvar(x = standata$rates_to, name = "rates_to", scode = "  int<lower=1, upper=N_states_hidden> rates_to[N_rates];", block = "data") +

  brms::stanvar(x = standata$corresponding_observation, name = "corresponding_observation", scode = "  int<lower=1, upper=N_states_observed> corresponding_observation[N_states_hidden];", block = "data") +

  brms::stanvar(x = standata$N_noisy_states, name = "N_noisy_states", scode = "  int<lower=0> N_noisy_states;", block = "data") +
  brms::stanvar(x = standata$noisy_states, name = "noisy_states", scode = "  int<lower=1, upper=N_states_observed> noisy_states[N_noisy_states];", block = "data") +
  brms::stanvar(x = standata$N_other_observations, name = "N_other_observations", scode = "  int<lower=1> N_other_observations;", block = "data") +
  brms::stanvar(x = standata$noisy_states_other_obs, name = "noisy_states_other_obs", scode = "  int<lower=1, upper=N_states_observed> noisy_states_other_obs[N_noisy_states, N_other_observations];", block = "data") +

  brms::stanvar(x = standata$sensitivity_low_bound, name = "sensitivity_low_bound", scode = "  real<lower=0, upper=1> sensitivity_low_bound;", block = "data") +


  brms::stanvar(x = standata$N_observations, name = "N_observations", scode = "  // Observations\n  int<lower=1> N_observations;", block = "data") +
  brms::stanvar(x = standata$N_series, name = "N_series", scode = "  int<lower=1> N_series;", block = "data") +
  brms::stanvar(x = standata$N_time, name = "N_time", scode = "  int<lower=1> N_time;", block = "data") +
  brms::stanvar(x = standata$N_predictor_sets, name = "N_predictor_sets", scode = "  int<lower=1> N_predictor_sets;", block = "data") +

  brms::stanvar(x = standata$initial_states, name = "initial_states", scode = "  int<lower=1, upper=N_states_hidden> initial_states[N_series];", block = "data") +
  brms::stanvar(x = standata$series, name = "series", scode = "  int<lower=1, upper=N_series> series[N_observations];", block = "data") +
  brms::stanvar(x = standata$times, name = "times", scode = "  int<lower=1, upper=N_time> times[N_observations];", block = "data") +
  brms::stanvar(x = standata$obs_states, name = "obs_states", scode = "      //0 for unobserved states
\n  int<lower=0, upper=N_states_observed> obs_states[N_observations];", block = "data") +

  brms::stanvar(x = standata$predictor_sets, name = "predictor_sets", scode = "  int<lower=1, upper=N_predictor_sets> predictor_sets[N_observations];", block = "data") +
  brms::stanvar(x = standata$rate_predictors, name = "rate_predictors", scode = "  int<lower=1, upper=N> rate_predictors[N_predictor_sets, N_rates];", block = "data")
}

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
  matrix[N_states_hidden, N_states_observed] observation_probs = compute_observation_matrix(
   N_states_hidden, N_states_observed, corresponding_observation,
    is_state_noisy, noisy_state_id, sensitivity, N_other_observations, noisy_states_other_obs,
    other_observations_probs);
  matrix[N_states_hidden, N_states_hidden] transition_matrices_t[N_predictor_sets];
  vector[N] exp_mu = exp(mu);

  //Compute all transition matrices
  for(ps in 1:N_predictor_sets) {
    vector[N_rates] rate_values = exp_mu[rate_predictors[ps]];
    transition_matrices_t[ps] =
      transpose(compute_transition_matrix(N_states_hidden,  N_rates, rates_from, rates_to, rate_values));
  }



  for(s in 1:N_series) {
    matrix[N_states_hidden, max_time[s]] alpha;
    vector[max_time[s]] alpha_log_norms;
    alpha[, 1] = rep_vector(0, N_states_hidden);
    alpha[initial_states[s], 1] = observation_probs[initial_states[s], obs_states_rect[s, 1]];

    {
      real norm = max(alpha[, 1]);
      alpha[, 1] /= norm;
      alpha_log_norms[1] = log(norm);
    }
    for(t in 2:max_time[s]) {
      real col_norm;
      int ps = predictor_sets_rect[s, t - 1];
      alpha[, t] = observation_probs[, obs_states_rect[s, t]] .* (transition_matrices_t[ps] * alpha[, t - 1]);
      col_norm = max(alpha[, t]);
      alpha[, t] /= col_norm;
      alpha_log_norms[t] = log(col_norm) + alpha_log_norms[t - 1];
    }
    target += log(sum(alpha[,max_time[s]])) + alpha_log_norms[max_time[s]];
  }

"
}

rate_hmm_genquant_code <- '
    matrix[N_states_hidden, N_states_observed] observation_probs = compute_observation_matrix(
     N_states_hidden, N_states_observed, corresponding_observation,
      is_state_noisy, noisy_state_id, sensitivity, N_other_observations, noisy_states_other_obs,
      other_observations_probs);
'

