hypotheses <-
  list(
       hcq_reduces_death = list(caption = "HCQ associated with risk of death", group = "hcq"),
       hcq_increases_discharged = list(caption = "HCQ associated with time in hospital", group = "hcq"),
       az_reduces_death = list(caption = "Azithromycin associated with risk of death", group = "az"),
       az_increases_discharged = list(caption = "Azithromycin associated with time in hospital", group = "az"),
       d_dimer_increases_death = list(caption = "High D-dimer associated with risk of death", group = "markers"),
       IL_6_increases_death = list(caption = "High IL-6 associated with risk of death", group = "markers")
  ) %>% purrr::imap(
         function(def, name) {
           def$name <- name
           def
         })

hypotheses_df <- purrr::map_dfr(hypotheses, as.data.frame) %>%
  mutate(name = factor(name, levels = name))

frequentist_hypothesis_res_from_coxph <- function(
  hypothesis, coxph_fit, coefficient_name, transition, adjusted) {

  summ <- summary(coxph_fit)
  coefficient_id <- summ$cmap[coefficient_name, transition]
  data.frame(hypothesis = hypothesis$name,
             model = "coxph",
             adjusted = adjusted,
             estimand = "log(HR)",
             p_value = summ$coefficients[coefficient_id, "Pr(>|z|)"],
             point_estimate = summ$coefficients[coefficient_id, "coef"],
             ci_low = log(summ$conf.int[coefficient_id, "lower .95"]),
             ci_high = log(summ$conf.int[coefficient_id, "upper .95"])
  )
}


bayesian_hypothesis_res_from_jm <- function(
  hypothesis, jm_fit, coefficient_name, adjusted) {
  bayesian_hypothesis_res_from_draws(
    draws = tidybayes::tidy_draws(jm_fit)[[coefficient_name]],
    model = "jm",
    estimand = "log(HR)",
    hypothesis = hypothesis,
    adjusted = adjusted)
}

bayesian_hypothesis_res_from_draws <- function(
  draws, model, estimand,
  hypothesis, adjusted) {

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
             row.names = NULL
  )
}
