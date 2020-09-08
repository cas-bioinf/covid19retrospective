hypotheses <-
  list(hcq_reduces_death = list(caption = "HCQ reduces risk of death"),
       hcq_increases_discharged = list(caption = "HCQ reduces time in hospital")
       ) %>% purrr::imap(
         function(def, name) {
           def$name <- name
           def
         })


frequentist_hypothesis_res_from_coxph <- function(
  hypothesis, coxph_fit, coefficient_name, transition, adjusted) {

  summ <- summary(coxph_fit)
  coefficient_id <- summ$cmap[coefficient_name, transition]
  data.frame(hypothesis = hypothesis$name,
             model = "coxph",
             adjusted = "age, sex",
             p_value = summ$coefficients[coefficient_id, "Pr(>|z|)"],
             ci_low = summ$conf.int[coefficient_id, "lower .95"],
             ci_high = summ$conf.int[coefficient_id, "upper .95"]
  )
}
