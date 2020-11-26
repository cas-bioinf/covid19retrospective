hypotheses <-
  list(
       hcq_death = list(caption = "HCQ associated with risk of death", group = "hcq", type = "death"),
       hcq_hospital = list(caption = "HCQ associated with time in hospital", group = "hcq", type = "hospital"),
       az_death = list(caption = "Azithromycin associated with risk of death", group = "az", type = "death"),
       az_hospital = list(caption = "Azithromycin associated with time in hospital", group = "az", type = "hospital"),
       favipiravir_death = list(caption = "Favipiravir associated with risk of death", group = "favipiravir", type = "death"),
       favipiravir_hospital = list(caption = "Favipiravir associated with time in hospital", group = "favipiravir", type = "hospital"),
       convalescent_plasma_death = list(caption = "Conv. plasma associated with risk of death", group = "convalescent_plasma", type = "death"),
       convalescent_plasma_hospital = list(caption = "Conv. plasma associated with time in hospital", group = "convalescent_plasma", type = "hospital"),
       d_dimer_death = list(caption = "High D-dimer associated with risk of death", group = "markers", type = "death"),
       IL_6_death = list(caption = "High IL-6 associated with risk of death", group = "markers", type = "death")
  ) %>% purrr::imap(
         function(def, name) {
           def$name <- name
           def
         })

hypotheses_df <- purrr::map_dfr(hypotheses, as.data.frame) %>%
  mutate(name = factor(name, levels = name))

frequentist_hypothesis_res_from_coxph <- function(
  hypothesis, coxph_fit, coefficient_name, transition, adjusted, model_check = "OK", multiplier = 1, model_name = "competing_coxph") {

  summ <- summary(coxph_fit)
  coefficient_id <- summ$cmap[coefficient_name, transition]
  data.frame(hypothesis = hypothesis$name,
             model = model_name,
             adjusted = adjusted,
             estimand = "log(HR)",
             p_value = summ$coefficients[coefficient_id, "Pr(>|z|)"],
             point_estimate = multiplier * summ$coefficients[coefficient_id, "coef"],
             ci_low = multiplier * log(summ$conf.int[coefficient_id, "lower .95"]),
             ci_high = multiplier * log(summ$conf.int[coefficient_id, "upper .95"]),
             model_check = model_check,
             data_version = get_data_version()
  )
}

# single outcome
frequentist_hypothesis_res_from_coxph1 <- function(
  hypothesis, adjusted, point_estimate, test_stat, df, p_value, ci_low, ci_high, model_check = "OK") {

  data.frame(hypothesis = hypothesis$name,
             model = "coxph_single",
             adjusted = adjusted,
             estimand = "log(HR)",
             point_estimate = point_estimate,
             p_value = p_value,
             ci_low = ci_low,
             ci_high = ci_high,
             model_check = model_check,
             data_version = get_data_version()
  )
}


bayesian_hypothesis_res_from_jm <- function(
  hypothesis, jm_fit, coefficient_name, adjusted, model_check = "OK") {
  bayesian_hypothesis_res_from_draws(
    draws = tidybayes::tidy_draws(jm_fit)[[coefficient_name]],
    model = "jm",
    estimand = "log(HR)",
    hypothesis = hypothesis,
    adjusted = adjusted,
    model_check = model_check)
}

bayesian_hypothesis_res_from_draws <- function(
  draws, model, estimand,
  hypothesis, adjusted, model_check = "OK") {

  reference_point <- 0
  reference_location <- ecdf(draws)(reference_point)
  data.frame(hypothesis = hypothesis$name,
             model = model,
             adjusted = adjusted,
             estimand = estimand,
             widest_CI_excl_reference = abs(reference_location - 0.5) * 2,
             point_estimate = mean(draws),
             ci_low = quantile(draws, 0.025),
             ci_high = quantile(draws, 0.975),
             model_check = model_check,
             data_version = get_data_version(),
             row.names = NULL
  )
}
