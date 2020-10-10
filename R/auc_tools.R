outcome_caption = c("alive_12" = "12-day mortality", "alive_30" = "30-day mortality",
                    "alive_all" = "Non-censored mortality",
                    "alive_imputed" = "Alive at data collection",
                    "not_ventilated_all" = "Not ventilated (uncensored only)",
                    "not_oxygen_all" = "Not on oxygen (uncensored only)",
                    "not_ventilated_admission" = "Not ventilated on admission",
                    "not_oxygen_admission" = "Not on oxygen on admission")
score_caption = c("ACP_grade" = "Li et al. (ACP)",
                  "chen_liu" = "Chen & Liu",
                  "chen_liu_incl" = "Chen & Liu",
                  "shi" = "Shi et al.",
                  "caramelo" = "Caramelo et al.",
                  "bello_chavolla", "Bello-Chavolla et al.")

outcome_note_caption = c("incl" = " (inclusion crit.)")

my_auc <- function(f, d, outcome_note = "", subgroup = "") {
  roc = pROC::roc(f, d, direction = "<", levels = c(TRUE, FALSE))
  ci = pROC::ci.auc(roc)

  outcome = as.character(formula.tools::lhs(f))
  score = as.character(formula.tools::rhs(f))


  list(roc = roc,
       summary = tibble::tibble(
         auc = ci[2],
         auc_low = ci[1],
         auc_high = ci[3],
         outcome = outcome,
         score = score,
         outcome_note = outcome_note,
         subgroup = subgroup
       )
  )
}

plot_my_auc <- function(auc_obj) {
  if(!(auc_obj$summary$score %in% names(score_caption))) {
    stop("Score does not have caption")
  }
  if(!(auc_obj$summary$outcome %in% names(outcome_caption))) {
    stop("Outcome does not have caption")
  }

  if(auc_obj$summary$outcome_note != "") {
    if(!(auc_obj$summary$outcome_note %in% names(outcome_note_caption))) {
      stop("Outcome note does not have caption")
    }
    outcome_note_caption_value <- outcome_note_caption[auc_obj$summary$outcome_note]
  } else {
    outcome_note_caption_value <- ""
  }


  caption = paste0(score_caption[auc_obj$summary$score], ", ", outcome_caption[auc_obj$summary$outcome],
                   outcome_note_caption_value)
  plot(auc_obj$roc, auc.polygon=TRUE, max.auc.polygon=TRUE, print.auc = TRUE, main = caption)
}

plot_all_auc <- function(my_auc_list) {
  my_auc_list %>% walk(plot_my_auc)
}

my_auc_from_estimate <- function(auc_estimate, outcome, score, outcome_note, subgroup = "") {
  list(roc = NULL,
       summary = tibble::tibble(
         auc = auc_estimate,
         auc_low = NA_real_,
         auc_high = NA_real_,
         outcome = outcome,
         score = score,
         outcome_note = outcome_note,
         subgroup = subgroup
       )
  )

}

compare_auc_to_orig <- function(orig_auc_list, our_auc_list, show_outcome = TRUE, ...) {
  get_summaries <- function(x) {
    map_df(x, ~ .x$summary)
  }

  source_factor <- function(x) {
    factor(x, levels = c("our","orig"), labels = c("Present study", "Original study"))
  }
  all_auc <- our_auc_list %>% map_df(get_summaries) %>%
    mutate(source = source_factor("our"), score = paste(score_caption[score], subgroup))
  all_orig <- orig_auc_list %>% get_summaries() %>% mutate(source = source_factor("orig"), score = paste(score_caption[score], subgroup))
  data_for_plot <- rbind(all_auc, all_orig) %>%
    mutate(
      id = interaction(outcome_note, outcome, source),
           size = if_else(source == "orig", 2, 1))

  res <- data_for_plot %>% ggplot(aes(x = auc, xmin = auc_low, xmax = auc_high, y = id, color = source, shape = source)) +
    geom_vline(aes(xintercept = auc),data = all_orig, color = "gray", linetype = "dashed") +
    geom_errorbarh(height = 0, na.rm = TRUE) + geom_point(aes(size = size)) +
    scale_size(range = c(2,4), guide = FALSE) +
    scale_x_continuous("C-statistic") +
    scale_shape_discrete("Source") +
    scale_color_discrete("Source") +
    facet_wrap(~score, scales = "free", ...)

  if(!show_outcome) {
    res <- res + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
  }

  res
}
