make_brms_formula_hmm <- function(formula) {
  if(!is.null(brms:::lhs(formula))) {
    stop("Formula needs to be one-sided")
  }

  brms::brmsformula(update.formula(formula, .dummy ~ .), family = rate_hmm_family)
}

make_stancode_hmm <- function(brmshmmdata) {
  d <- validate_brmshmmdata(brmshmmdata)

  data <- make_data_hmm(d)
  brms::make_stancode(d$formula, family = rate_hmm_family,
                      data = data$brmsdata,
                      stanvars = rate_hmm_stanvars(data$standata), prior = d$prior)
}

rate_hmm_family <- function() {
  structure(list(family = "rate_hmm",
                                  link = "identity", dpars = "mu",
                           lb = NA, ub = NA, type = "real", vars = NULL, specials = NULL,
                           ybounds = c(-Inf, Inf),
                           log_lik = NULL, posterior_predict = NULL, posterior_epred = NULL),
                      #class = c("rate_hmm", "brmsfamily", "family"))
                      class = c("rate_hmm", "brmsfamily"))
}

summary.rate_hmm <- function(f, link = FALSE) {
  "Rate-based HMM"
}

.family_rate_hmm <- function() {

}

# family_info.rate_hmm <- function(x, y, ...) {
#   NULL
# }

posterior_epred_rate_hmm <- function(prep) {
  brms:::posterior_epred_gaussian(prep)
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
  // The first matrix index is 'from', second is 'to'
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


  brms::stanvar(x = standata$N_series, name = "N_series", scode = "  // Observations\n  int<lower=1> N_series;", block = "data") +
  brms::stanvar(x = standata$N_time, name = "N_time", scode = "  int<lower=1> N_time;", block = "data") +
  brms::stanvar(x = standata$N_predictor_sets, name = "N_predictor_sets", scode = "  int<lower=1> N_predictor_sets;", block = "data") +

  brms::stanvar(x = standata$initial_states, name = "initial_states", scode = "  int<lower=1, upper=N_states_hidden> initial_states[N_series];", block = "data") +
  brms::stanvar(x = standata$serie_max_time, name = "serie_max_time", scode = "  int<lower=1, upper=N_time> serie_max_time[N_series];", block = "data") +
  brms::stanvar(x = standata$obs_states, name = "obs_states_rect", scode = "      //0 for unobserved states
\n  int<lower=0, upper=N_states_observed> obs_states_rect[N_series, N_time];", block = "data") +

  brms::stanvar(x = standata$predictor_sets_rect, name = "predictor_sets_rect", scode = "  int<lower=0, upper=N_predictor_sets> predictor_sets_rect[N_series, N_time];", block = "data") +
  brms::stanvar(x = standata$rate_predictors, name = "rate_predictors", scode = "  int<lower=1, upper=N> rate_predictors[N_predictor_sets, N_rates];", block = "data") +
  brms::stanvar(x = standata$optimize_possible, name = "optimize_possible", scode = "  int<lower=0, upper=2> optimize_possible;", block = "data")
}

rate_hmm_tdata_code <- '
  int<lower=0,upper=1> is_state_noisy[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=0,upper=N_noisy_states> noisy_state_id[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=1,upper=N_states_hidden> N_possible_states[N_states_observed] = rep_array(0, N_states_observed);
  int<lower=0,upper=N_states_hidden> possible_states[N_states_observed, N_states_hidden] = rep_array(0, N_states_observed, N_states_hidden);

  for(s_index in 1:N_noisy_states) {
    is_state_noisy[noisy_states[s_index]] = 1;
    noisy_state_id[noisy_states[s_index]] = s_index;
  }

  //Compute possible hidden states for each observation
  //This let\'s us optimize some computation for the case where not all hidden states are possible
  {
    int is_state_possible[N_states_observed, N_states_hidden] = rep_array(0, N_states_observed, N_states_hidden);
    for(s_hidden in 1:N_states_hidden) {
      int corresponding_obs = corresponding_observation[s_hidden];
      is_state_possible[corresponding_obs, s_hidden] = 1;
      if(is_state_noisy[corresponding_obs]) {
        int noisy_id = noisy_state_id[corresponding_obs];
        for(other_index in 1:N_other_observations) {
          int s_observed = noisy_states_other_obs[noisy_id, other_index];
          is_state_possible[s_observed, s_hidden] = 1;
        }
      }
    }
    for(s_observed in 1:N_states_observed) {
      int next_possible_state = 1;
      N_possible_states[s_observed] = sum(is_state_possible[s_observed, ]);
      for(s_hidden in 1:N_states_hidden) {
        if(is_state_possible[s_observed, s_hidden]) {
          possible_states[s_observed, next_possible_state] = s_hidden;
          next_possible_state += 1;
        }
      }
    }
  }


  //Check data validity
  for(p in 1:N_series) {
    for(t in 1:(serie_max_time[p] - 1)) {
      if(predictor_sets_rect[p, t] == 0) {
        reject("serie ", p, " has missing predictors for time ", t);
      }
    }
  }
'

##' @export rate_hmm_parameters_code
rate_hmm_parameters_code <- "
  vector<lower=sensitivity_low_bound, upper=1>[N_noisy_states] sensitivity;
  simplex[N_other_observations] other_observations_probs[N_noisy_states];
"

stan_log_lik.rate_hmm <- function( ... ) {
"
  matrix[N_states_hidden, N_states_observed] observation_probs = compute_observation_matrix(
   N_states_hidden, N_states_observed, corresponding_observation,
    is_state_noisy, noisy_state_id, sensitivity, N_other_observations, noisy_states_other_obs,
    other_observations_probs);
  matrix[N_states_hidden, N_states_hidden] transition_matrices[N_predictor_sets];
  matrix[N_states_hidden, N_states_hidden] transition_matrices_t[N_predictor_sets];
  vector[N] exp_mu = exp(mu);

  //Compute all transition matrices
  for(ps in 1:N_predictor_sets) {
    vector[N_rates] rate_values = exp_mu[rate_predictors[ps]];
    transition_matrices[ps] =
      compute_transition_matrix(N_states_hidden,  N_rates, rates_from, rates_to, rate_values);
    transition_matrices_t[ps] =
      transpose(transition_matrices[ps]);
  }



  for(s in 1:N_series) {
    matrix[N_states_hidden, serie_max_time[s]] alpha = rep_matrix(0, N_states_hidden, serie_max_time[s]);
    vector[serie_max_time[s]] alpha_log_norms;
    alpha[initial_states[s], 1] = observation_probs[initial_states[s], obs_states_rect[s, 1]];


    {
      real norm = max(alpha[, 1]);
      alpha[, 1] /= norm;
      alpha_log_norms[1] = log(norm);
    }
    for(t in 2:serie_max_time[s]) {
      real col_norm;
      int ps = predictor_sets_rect[s, t - 1];
      int obs = obs_states_rect[s,t];
      if(obs != 0) {
        //Oberved something.
        int N_possible = N_possible_states[obs];
        if(!optimize_possible || N_possible == N_states_hidden) {
          vector[N_states_hidden] transition_probs =
            transition_matrices_t[ps] * alpha[, t - 1];
          alpha[, t] = observation_probs[, obs] .* transition_probs;
        } else {
          //Optimization: only update states that are possible
          int obs_prev = obs_states_rect[s,t - 1];
          int N_possible_prev = obs_prev == 0 ? N_states_hidden : N_possible_states[obs_prev];
          if(optimize_possible == 1 || N_possible_prev == N_states_hidden) {
            for(p_state_id in 1:N_possible) {
               int p_state = possible_states[obs, p_state_id];
               alpha[p_state, t] = observation_probs[p_state, obs] * dot_product(alpha[, t-1], transition_matrices[ps, , p_state]);
            }
          } else {
            for(p_state_id in 1:N_possible) {
               int p_state = possible_states[obs, p_state_id];
               vector[N_possible_prev] transition_probs;
               for(p_state_id_prev in 1:N_possible_prev) {
                 int p_state_prev = possible_states[obs_prev, p_state_id_prev];
                 transition_probs[p_state_id_prev] = alpha[ p_state_prev, t - 1] * transition_matrices[ps, p_state_prev , p_state];
               }
               alpha[p_state, t] = observation_probs[p_state, obs] * sum(transition_probs);
            }
          }
        }
      } else {
        //No observation
        vector[N_states_hidden] transition_probs = (transition_matrices_t[ps] * alpha[, t - 1]);
        alpha[, t] = transition_probs;
      }
      col_norm = max(alpha[, t]);
      alpha[, t] /= col_norm;
      alpha_log_norms[t] = log(col_norm) + alpha_log_norms[t - 1];
    }
    target += log(sum(alpha[,serie_max_time[s]])) + alpha_log_norms[serie_max_time[s]];
  }

"
}


rate_hmm_genquant_code <- '
    matrix[N_states_hidden, N_states_observed] observation_probs = compute_observation_matrix(
     N_states_hidden, N_states_observed, corresponding_observation,
      is_state_noisy, noisy_state_id, sensitivity, N_other_observations, noisy_states_other_obs,
      other_observations_probs);
'

