functions {
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

}

data {
  // This is to emulate brms
  int N;
  vector[N] intercept_prior_logmean;
  vector<lower=0>[N] intercept_prior_sd;

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
  int<lower=1> N_patients;
  int<lower=1> N_time;
  int<lower=1> N_predictor_sets;

  int<lower=1, upper=N_patients> patients[N_observations];
  int<lower=1, upper=N_time> times[N_observations];
  //0 for unobserved states
  int<lower=0, upper=N_states_observed> obs_states[N_observations];

  int<lower=1, upper=N_predictor_sets> predictor_sets[N_observations];
  int<lower=1, upper=N> rate_predictors[N_predictor_sets, N_rates];
}

transformed data {
  int<lower=0,upper=1> is_state_noisy[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=0,upper=N_noisy_states> noisy_state_id[N_states_observed] = rep_array(0, N_states_observed);

  //Rectangule observations and predictors. 0 for missing data
  int<lower=0, upper=N_states_observed> obs_states_rect[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0, upper=N_predictor_sets> predictor_sets_rect[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0, upper=N_time> max_time[N_patients] = rep_array(0,N_patients);

  for(s_index in 1:N_noisy_states) {
    is_state_noisy[noisy_states[s_index]] = 1;
    noisy_state_id[noisy_states[s_index]] = s_index;
  }

  for(o in 1:N_observations) {
    obs_states_rect[patients[o], times[o]] = obs_states[o];
    predictor_sets_rect[patients[o], times[o]] = predictor_sets[o];
    if(obs_states[o] != 0) {
      max_time[patients[o]] = max(max_time[patients[o]], times[o]);
    }
  }

  //Check data validity
  for(p in 1:N_patients) {
    if(max_time[p] == 0) {
      reject("Patient ", p, " has no osbservations");
    }
    for(t in 1:(max_time[p] - 1)) {
      if(predictor_sets_rect[p, t] == 0) {
        reject("Patient ", p, " has missing predictors for time ", t);
      }
    }
  }
}

parameters {
  vector<lower=sensitivity_low_bound, upper=1>[N_noisy_states] sensitivity;
  simplex[N_other_observations] other_observations_probs[N_noisy_states];

  //emulating brms
  vector[N] mu;
}

model {
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

  for(p in 1:N_patients) {
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

  // Priors
  mu ~ normal(intercept_prior_logmean, intercept_prior_sd);
}

