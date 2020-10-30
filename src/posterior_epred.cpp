// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_outerproduct(const Eigen::VectorXd & x) {
  Eigen::MatrixXd m = x * x.transpose();
  return m;
}


// [[Rcpp::export]]
Eigen::VectorXd posterior_epred_state_prob_internal(int n_samples, const Eigen::VectorXi &max_times,
                                                    const Eigen::VectorXi &initial_states,
                                                    int n_states,
                                                    int n_time,
                                                    int n_predictor_sets,
                                                    const Eigen::VectorXd &transition_matrices,
                                                    const Eigen::MatrixXi &predictor_sets_rect) {

  using namespace Eigen;

  int n_series = max_times.size();
  if(initial_states.size() != n_series) {
    throw Rcpp::exception("n_series");
  }
  if(predictor_sets_rect.rows() != n_series){
    throw Rcpp::exception("n_series predictor_sets_rect");
  }

  if(predictor_sets_rect.cols() != n_time) {
    throw Rcpp::exception("n_time predictor_sets_rect");
  };
  if(transition_matrices.size() != n_states * n_states * n_samples * n_predictor_sets) {
    throw Rcpp::exception("Transition matrices dim");
  }

  std::vector<std::vector<MatrixXd>> transition_matrices_converted(n_samples);
  for(size_t s = 0; s < n_samples; ++s) {
    transition_matrices_converted[s].reserve(n_predictor_sets);
    for(size_t ps = 0; ps < n_predictor_sets; ++ps) {
      //MatrixXd tmp(n_states, n_states);
      size_t start_index = s * n_states * n_states * n_predictor_sets + ps * n_states * n_states;

      // for(size_t i = 0; i < n_states * n_states; ++i) {
      //   tmp << transition_matrices(seqN(start_index, n_states * n_states));
      // }
      Map<const MatrixXd> tmp(&transition_matrices(start_index), n_states, n_states);
      transition_matrices_converted[s].push_back(tmp.transpose());
    }
  }


  auto res_index = [n_states, n_samples, n_series, n_time](size_t state, size_t time, size_t sample, size_t serie) {
    return serie * (n_samples * n_states * n_time) + sample * (n_states * n_time) + time * n_states + state;
  };

  VectorXd res = VectorXd::Constant(n_states * n_samples * n_series * n_time, R_NaN);

  for(size_t serie = 0; serie < n_series; ++serie) {
    for(size_t sample = 0; sample < n_samples; ++ sample) {
      VectorXd state_probs = VectorXd::Zero(n_states);
      state_probs[initial_states(serie) - 1] = 1.0;

      //Can't get this to work using a loop instead
      //res(seqN(res_index(0, 0, sample, serie), n_states)) << state_probs;
      {
        size_t start_index = res_index(0, 0, sample, serie);
        for(size_t i = 0; i < n_states; ++i) {
          res(start_index + i) = state_probs[i];
        }
      }
      for(size_t time = 1; time < max_times[serie]; ++time) {
        int predictor_set = predictor_sets_rect(serie, time) - 1;
        //const MatrixXd &t_matrix = transition_matrices_converted[sample][predictor_set];
        if(sample >= transition_matrices_converted.size() || predictor_set >=
           transition_matrices_converted[sample].size()) {
          Rcpp::Rcout << "R1" << std::endl;
          throw Rcpp::exception("Range");
        }
        state_probs = transition_matrices_converted[sample][predictor_set] * state_probs;

        //res(seqN(res_index(0, time, sample, serie), n_states)) << state_probs;
        {
          size_t start_index = res_index(0, time, sample, serie);
          for(size_t i = 0; i < n_states; ++i) {
            if(start_index + i >= res.size()) {
              Rcpp::Rcout << std::endl << "R2: " << start_index << " " << res.size() << " "  << i << std::endl;
              throw Rcpp::exception("Range2");
            }
            res(start_index + i) = state_probs[i];
          }
        }
      }
    }
  }

  return res;
}
