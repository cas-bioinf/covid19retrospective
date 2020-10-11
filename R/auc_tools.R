outcome_caption = c("alive_12" = "12-day mortality", "alive_30" = "30-day mortality",
                    "alive_all" = "Non-censored mortality",
                    "alive_imputed" = "Alive at data collection",
                    "not_ventilated_all" = "Not ventilated (uncensored only)",
                    "not_oxygen_all" = "Not on oxygen (uncensored only)",
                    "not_ventilated_admission" = "Not ventilated on admission",
                    "not_oxygen_admission" = "Not on oxygen on admission")
score_caption = c("ACP_grade" = "Li et al. (ACP)",
                  "chen_liu" = "Chen & Liu",
                  "shi" = "Shi et al.",
                  "caramelo" = "Caramelo et al.",
                  "bc_base" = "Bello-Chavolla et al.",
                  "bc0" = "Bello-Chavolla et al.",
                  "bc1" = "Bello-Chavolla et al.",
                  "caramelo_base" = "Caramelo et al.",
                  "caramelo1" = "Caramelo et al.",
                  "caramelo2" = "Caramelo et al.",
                  "age" = "Age only"
)

outcome_note_caption = c("incl" = " (inclusion crit.)")

my_auc_summary <- function(f, d, outcome_note = "", subgroup = "") {
  roc = pROC::roc(f, d, direction = "<", levels = c(TRUE, FALSE))
  ci = pROC::ci.auc(roc)

  outcome = as.character(formula.tools::lhs(f))
  score = as.character(formula.tools::rhs(f))


  tibble::tibble(
         auc = ci[2],
         auc_low = ci[1],
         auc_high = ci[3],
         outcome = outcome,
         score = score,
         outcome_note = outcome_note,
         subgroup = subgroup
   )
}

my_auc_summary_from_estimate <- function(auc_estimate, outcome, score, outcome_note = "", subgroup = "") {
  tibble::tibble(
    auc = auc_estimate,
    auc_low = NA_real_,
    auc_high = NA_real_,
    outcome = outcome,
    score = score,
    outcome_note = outcome_note,
    subgroup = subgroup
  )
}


plot_my_auc <- function(f, d, outcome_note = "") {

  outcome = as.character(formula.tools::lhs(f))
  score = as.character(formula.tools::rhs(f))

  if(!(score %in% names(score_caption))) {
    stop("Score does not have caption")
  }
  if(!(outcome %in% names(outcome_caption))) {
    stop("Outcome does not have caption")
  }

  if(outcome_note != "") {
    if(!(outcome_note %in% names(outcome_note_caption))) {
      stop("Outcome note does not have caption")
    }
    outcome_note_caption_value <- outcome_note_caption[outcome_note]
  } else {
    outcome_note_caption_value <- ""
  }


  caption = paste0(score_caption[score], ", ", outcome_caption[outcome],
                   outcome_note_caption_value)

  roc = pROC::roc(f, d, direction = "<", levels = c(TRUE, FALSE))

  plot(roc, auc.polygon=TRUE, max.auc.polygon=TRUE, print.auc = TRUE, main = caption)
}

plot_all_auc <- function(my_auc_list) {
  my_auc_list %>% walk(plot_my_auc)
}


compare_auc_to_orig <- function(orig_auc, our_auc, show_outcome = TRUE, ...) {

  source_factor <- function(x) {
    factor(x, levels = c("our","orig"), labels = c("Present study", "Original study"))
  }
  all_auc <- our_auc %>%
    mutate(source = source_factor("our"), score = paste(score_caption[score], subgroup))

  if(nrow(orig_auc) > 0) {
    all_orig <- orig_auc %>% mutate(source = source_factor("orig"), score = paste(score_caption[score], subgroup))
    data_for_plot_base <- rbind(all_auc, all_orig)
  } else {
    data_for_plot_base <- all_auc
    all_orig <- all_auc %>% head(0) # Just to enforce correct columns
  }
  data_for_plot <- data_for_plot_base %>%
    mutate(
      id = interaction(outcome_note, outcome, source),
           size = if_else(source == "orig", 2, 1))

  res <- data_for_plot %>% ggplot(aes(y = auc, ymin = auc_low, ymax = auc_high, x = id, color = source, shape = source)) +
    geom_hline(yintercept = 0.5, color = "darkblue", size = 2) +
    geom_hline(aes(yintercept = auc),data = all_orig, color = "gray", linetype = "dashed") +
    geom_errorbar(width = 0, na.rm = TRUE) + geom_point(aes(size = size)) +
    scale_size(range = c(2,4), guide = FALSE) +
    scale_y_continuous("AUC") +
    scale_shape_discrete("Source") +
    scale_color_discrete("Source") +
    facet_wrap(~score, scales = "free_x", ...) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

  if(!show_outcome) {
    res <- res + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  res
}
