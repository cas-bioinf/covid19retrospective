test_that("known_epred", {
  ### Fully deterministic transitions
  deterministic_transition_matrix <- matrix(0, nrow = 5, ncol = 5)
  for(i in 1:4) {
    deterministic_transition_matrix[i, i + 1] <- 1
  }
  deterministic_transition_matrix[5,5] <- 1
  transition_matrices <- array(NA_real_, c(5,5,1,1))
  transition_matrices[,,1,1] <- deterministic_transition_matrix
  N_series <- 10

  deterministic_args <- list(nsamples = 1, max_times = rep(6, N_series), initial_states = rep(1, N_series),
                             transition_matrices = transition_matrices,
                             predictor_sets_rect = array(1L, c(N_series, 6)))
  results_simulate <- do.call(posterior_epred_rect_simulate, deterministic_args)

  results_state_prob <- do.call(posterior_epred_state_prob, deterministic_args)

  expect_true(all(results_state_prob >= 0 & results_state_prob <= 1))
  for(s in 1:N_series) {
    for(t in 1:6) {
      expect_equal(sum(results_state_prob[, t, 1, s]), 1)
    }
  }

  for(i in 1:5) {
    expect_equal(results_simulate[i,,1], rep(i, N_series))
    expect_equal(results_state_prob[i, i, 1, ], rep(1, N_series))
  }
  expect_equal(results_simulate[6,,1], rep(5, N_series))
  expect_equal(results_state_prob[5, 6, 1, ], rep(1, N_series))


  ### Analytically solvable but stochastic matrix

  N_samples <- 1000
  N_series <- 4
  geometric_transition_matrix_a <- matrix(0, nrow = 3, ncol = 3)
  geometric_transition_matrix_a[1,1] <- 0.5
  geometric_transition_matrix_a[1,2] <- 0.5
  geometric_transition_matrix_a[2,2] <- 0.3
  geometric_transition_matrix_a[2,3] <- 0.7
  geometric_transition_matrix_a[3,3] <- 1

  geometric_transition_matrix_b <- matrix(0, nrow = 3, ncol = 3)
  geometric_transition_matrix_b[1,1] <- 0.2
  geometric_transition_matrix_b[1,2] <- 0.8
  geometric_transition_matrix_b[2,2] <- 0.1
  geometric_transition_matrix_b[2,3] <- 0.9
  geometric_transition_matrix_b[3,3] <- 1

  transition_matrices <- array(NA_real_, c(3,3,2,N_samples))
  for(i in 1:N_samples){
    transition_matrices[,,1,i] <- geometric_transition_matrix_a
    transition_matrices[,,2,i] <- geometric_transition_matrix_b
  }
  predictor_sets_rect <- array(1L, c(N_series, 6))
  predictor_sets_rect[c(2,4),] <- 2

  geometric_args <- list(nsamples = N_samples, max_times = rep(6, N_series),
                         initial_states = c(1,1,2,2),
                         transition_matrices = transition_matrices,
                         predictor_sets_rect = predictor_sets_rect)
  results_simulate <- do.call(posterior_epred_rect_simulate, geometric_args)
  results_state_prob <- do.call(posterior_epred_state_prob, geometric_args)

  tolerance_state_probs <- 1e-6
  expect_true(all(results_state_prob >= 0 & results_state_prob <= 1))
  for(s in 1:N_series) {
    for(t in 1:6) {
        expect_true(all(abs(colSums(results_state_prob[, t, , s]) - 1) < tolerance_state_probs))
    }
  }


  state_probs_simulate <- array(NA_real_, c(N_series, 6, 3))
  for(time in 1:6) {
    for(serie in 1:4) {
      for(state in 1:3) {
        state_probs_simulate[serie, time, state] <- mean(results_simulate[time, serie, ] == state)
      }
    }
  }
  expect_true(all(state_probs_simulate[,1,3] == 0))
  expect_true(all(state_probs_simulate[c(3,4),,1] == 0))
  expect_true(all(results_state_prob[3, 1, ,] == 0))
  expect_true(all(results_state_prob[1,,,c(3,4)] == 0))

  tolerance <- 0.05
  for(time in 2:6) {
    expect_true(abs(state_probs_simulate[1, time, 1] - 0.5^(time - 1)) < tolerance)
    expect_true(all(abs(results_state_prob[1, time, , 1] - 0.5^(time - 1)) < tolerance_state_probs))
    expect_true(abs(state_probs_simulate[2, time, 1] - 0.2^(time - 1)) < tolerance)
    expect_true(all(abs(results_state_prob[1, time, , 2] - 0.2^(time - 1)) < tolerance_state_probs))
    expect_true(abs(state_probs_simulate[3, time, 2] - 0.3^(time - 1)) < tolerance)
    expect_true(all(abs(results_state_prob[2, time, , 3] - 0.3^(time - 1)) < tolerance_state_probs))
    expect_true(abs(state_probs_simulate[4, time, 2] - 0.1^(time - 1)) < tolerance)
    expect_true(all(abs(results_state_prob[2, time, , 4] - 0.1^(time - 1)) < tolerance_state_probs))
  }

  expect_true(abs(state_probs_simulate[1, 3, 3] - 0.5 * 0.7) < tolerance)
  expect_true(all(results_state_prob[3, 3, , 1] == 0.5 * 0.7))
  expect_true(abs(state_probs_simulate[2, 3, 3] - 0.8* 0.9) < tolerance)
  expect_true(all(results_state_prob[3, 3, , 2] == 0.8 * 0.9))

})
