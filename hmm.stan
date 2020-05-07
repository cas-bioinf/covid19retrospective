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
  int<lower=1> N_patients;
  int<lower=1> N_states_hidden;
  int<lower=1> N_states_observed;
  int<lower=1> N_time;

  int<lower=1> N_rates;
  int<lower=1, upper=N_states_hidden> rates_from[N_rates];
  int<lower=1, upper=N_states_hidden> rates_to[N_rates];
  vector<lower=0>[N_rates] rates_prior_alpha;
  vector<lower=0>[N_rates] rates_prior_beta;

  int<lower=1, upper=N_states_observed> corresponding_observation[N_states_hidden];

  int<lower=0> N_noisy_states;
  int<lower=1, upper=N_states_observed> noisy_states[N_noisy_states];
  int<lower=1> N_other_observations;
  int<lower=1, upper=N_states_observed> noisy_states_other_obs[N_noisy_states, N_other_observations];

  real<lower=0, upper=1> sensitivity_low_bound;

  simplex[N_states_hidden] initial_states_prob;

  int<lower=1, upper=N_states_observed> observations[N_time, N_patients];
}

transformed data {
  int<lower=0,upper=1> is_state_noisy[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=0,upper=N_noisy_states> noisy_state_id[N_states_observed] = rep_array(0, N_states_observed);

  for(s_index in 1:N_noisy_states) {
    is_state_noisy[noisy_states[s_index]] = 1;
    noisy_state_id[noisy_states[s_index]] = s_index;
  }
}

parameters {
  vector<lower=0>[N_rates] rates;
  vector<lower=sensitivity_low_bound, upper=1>[N_noisy_states] sensitivity;
  simplex[N_other_observations] other_observations_probs[N_noisy_states];
}

model {
  matrix[N_states_hidden, N_states_observed] observation_probs = rep_matrix(0, N_states_hidden, N_states_observed);
  matrix[N_states_hidden, N_states_hidden] transition_matrix_t =
    transpose(compute_transition_matrix(N_states_hidden,  N_rates, rates_from, rates_to, rates));

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
    matrix[N_states_hidden, N_time] alpha;
    vector[N_time] alpha_log_norms;
    alpha[, 1] = initial_states_prob .* observation_probs[, observations[1, p]];

    {
      real norm = max(alpha[, 1]);
      alpha[, 1] /= norm;
      alpha_log_norms[1] = log(norm);
    }
    for(t in 2:N_time) {
      real col_norm; //TODO merge with col_norm
      alpha[, t] = observation_probs[, observations[t, p]] .* (transition_matrix_t * alpha[, t - 1]);
      col_norm = max(alpha[, t]);
      alpha[, t] /= col_norm;
      alpha_log_norms[t] = log(col_norm) + alpha_log_norms[t - 1];
    }
    target += log(sum(alpha[,N_time])) + alpha_log_norms[N_time];
  }

  // Priors
  rates ~ gamma(rates_prior_alpha, rates_prior_beta);
}

