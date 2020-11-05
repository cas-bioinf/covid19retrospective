plot_multiverse <- function(hypo_res, adjustments_to_hide = c(), x_limits = NULL) {
  all_adjustments_dict <- c(
                       "hospital" = "Site",
                       "admitted" = "Admitted",
                       "supportive" = "Supportive",
                       "age" = "Age",
                       "age_norm" = "Age",
                       "sex" = "Sex",
                       "smoking" = "Smoking",
                       "obesity" = "Obesity",
                       "comorbidities \\(sum\\)" = "Comorb.",
                       "comorbidities" = "Comorb.",
                       "hcq" = "HCQ",
                       "az" = "AZ",
                       "favipiravir" = "FPV",
                       "all treatments" = "All treatments",
                       "convalescent_plasma" = "C. plasma")
  all_adjustments = names(all_adjustments_dict)

  model_labels <- c("brms_categorical" = "Categorical (all)",
                    "brms_categorical_7" = "Categorical 7",
                    "brms_categorical_28" = "Categorical 28",
                    "brms_cox" = "Bayes Cox",
                    "coxph1" = "Cox (single)",
                    "coxph" = "Cox (competing)",
                    "HMMsimple_group" = "HMM A",
                    "HMMsimple_rate" = "HMM B",
                    "HMMcomplex_group" = "HMM C",
                    "jm" = "JM"
  )

  single_adjustment_regex <- paste(all_adjustments, collapse = "|")
  adjustment_regex <- paste0("^((",single_adjustment_regex,")(, |$))+")
  recognized_adjustments <- grepl(adjustment_regex, hypo_res$adjusted)
  if(!all(recognized_adjustments)) {
    print(hypo_res$adjusted[!recognized_adjustments])
    stop("Unexpected adjustment")
  }

  spacers <- hypo_res %>% group_by(model, hypothesis, estimand, caption, group, type) %>%
    summarise(is_spacer = TRUE, point_estimate = NA, ci_low = NA, ci_high = NA, adjusted = NA,
              model_check = NA, data_version = NA)

  for(adj in all_adjustments) {
    hypo_res[[paste0("adj.", adj)]] <- grepl(paste0("(, |^)", adj,"(,|$)"), hypo_res$adjusted)
    spacers[[paste0("adj.", adj)]] <- NA
  }




  my_facet <- facet_wrap(~caption, ncol = 1, scales = "free_y")

  # Clean up adjustments
  hypo_res <- hypo_res %>%
    mutate(is_spacer = FALSE) %>%
    rbind(spacers) %>%
    mutate(model_label = factor(model_labels[model], levels = model_labels)) %>%
    mutate(id_tmp = 1:n()) %>%
    arrange(model_label, id_tmp) %>%
    select(-id_tmp) %>%
    mutate(
       #id = factor(paste0(model, "__", adjusted), levels = rev(unique(paste0(model, "__", adjusted)))),
       id = factor(1:n(), levels = rev(1:n())),
       adj.age = adj.age | adj.age_norm,
       adj.comorbidities = adj.comorbidities | `adj.comorbidities \\(sum\\)`,
       adj.hcq = adj.hcq | `adj.all treatments`,
       adj.az = adj.az | `adj.all treatments`,
       adj.favipiravir = adj.favipiravir | `adj.all treatments`,
       adj.convalescent_plasma = adj.convalescent_plasma | `adj.all treatments`
       #adj.other_treatments = `adj.all treatments`
    ) %>%
    rename(full_adjusted = adjusted) %>%
             select(-adj.age_norm, -`adj.comorbidities \\(sum\\)`,
                    -`adj.all treatments`)


  my_labeller <- function(breaks) {
    hypo_res_for_labels <- hypo_res %>%
      group_by(model, hypothesis) %>%
      mutate(model = if_else(n() == 1 | as.integer(id) == max(as.integer(id)), as.character(model_label), "")) %>%
      ungroup()
    hypo_res_for_labels$model[as.integer(breaks)]
  }

  spacers_geom <- geom_hline(data = hypo_res %>% filter(is_spacer), color = "gray", linetype = "dashed",
                             aes(yintercept = id))

  shared_theme <- theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())

  if(!is.null(x_limits)) {
    limits_addition <- expand_limits(x = x_limits)
  } else {
    limits_addition <- NULL
  }

  plot_adj <- hypo_res %>%
    pivot_longer(cols = starts_with("adj."), names_prefix = "adj.", names_to = "adjustment",
                 values_to = "is_adjusted") %>%
    filter(!(adjustment %in% adjustments_to_hide)) %>%
    mutate(adjustment = factor(adjustment, levels = all_adjustments, labels = all_adjustments_dict)) %>%
    ggplot(aes(x = adjustment, y = id, color = is_adjusted)) +
      #geom_tile() + scale_fill_discrete(guide = FALSE, na.color = "white") +
      geom_point(na.rm = TRUE, shape = 15) +
      spacers_geom +
      scale_color_manual(values = c("white", "black"), guide = FALSE) +
      scale_y_discrete("Model type", labels = my_labeller) +
      scale_x_discrete("Adjustment") +
      my_facet +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.2, size = 9),
            axis.text.y = element_text(size = 8),
            strip.text = element_text(color = "white"),
            strip.background = element_blank(),
            plot.margin = margin(r = 2)) +
    shared_theme

  plot_cis <- hypo_res %>%
   ggplot(aes(x = point_estimate, xmin = ci_low, xmax = ci_high, y = id)) +
    geom_vline(xintercept = 0, color = "blue") +
    geom_vline(xintercept = log(1.5), color = "green") +
    geom_vline(xintercept = -log(1.5), color = "green") +
    geom_point(na.rm = TRUE) +
    geom_errorbarh(na.rm = TRUE) +
    spacers_geom +
    scale_x_continuous("Model coefficient estimate") +
    scale_y_discrete("") +
    my_facet +
    limits_addition +
    theme(axis.text.y = element_blank(),
          plot.margin = margin(l = 0)) +
    shared_theme

  (plot_adj | plot_cis ) + plot_layout(widths = c(1.8,4)) & theme(axis.title = element_text(size = 11, face = "bold"))
}
