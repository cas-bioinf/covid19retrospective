model_color <- "#3993dd"
data_color <- "#91171f"
decoration_color <- "#ca5310"

posterior_state_plot <- function(predictions_df, series_id = unique(predictions_df$.serie),
                                 observed_data = NULL,
                                 predicted_labels = NULL) {

  if(is.factor(predictions_df$.predicted) && !is.null(predicted_labels)) {
    warning("predicted_labels ignored as .predicted is already a factor")
  }

  if(!is.factor(predictions_df$.predicted)) {
    if(!is.null(observed_data) && is.null(predicted_labels)) {
      predicted_labels <- levels(observed_data$serie_data$.observed)
    }

    if(!is.null(predicted_labels)) {
      predicted_levels <- 1:length(predicted_labels)
    } else {
      predicted_levels <- 1:max(predictions_df$.predicted)
      predicted_labels <- as.character(predicted_levels)
    }
  }

  predictions_grouped <- predictions_df %>%
    select(-`.observed`) %>%
    filter(.serie %in% series_id) %>%
    group_by_at(setdiff(names(predictions_df), c(".sample", ".observed"))) %>%
    summarise(n_samples = n()) %>%
    group_by_at(setdiff(names(predictions_df), c(".observed", ".predicted", ".sample"))) %>%
    mutate(state_prob = n_samples / sum(n_samples),
           .predicted = factor(.predicted, levels = predicted_levels, labels = predicted_labels)) %>%
    ungroup()

  if(!is.null(observed_data)) {
    observed_data <- validate_brmshmmdata(observed_data)
    observed_data_plot <- observed_data$serie_data %>% filter(.serie %in% series_id, !is.na(.observed))
    observed_geom <- geom_line(data = observed_data_plot, aes(x = .time, y = .observed, group = .serie), color = data_color, inherit.aes = FALSE, size = 2)
  } else {
    observed_geom <- NULL
  }

  predictions_grouped %>%
          ggplot(aes(x = .time, y = .predicted, fill = state_prob)) +
          geom_raster() +
          observed_geom +
          scale_fill_gradient(low = "white", high = model_color, limits = c(0,1)) +
          facet_wrap(~.serie)
}

pp_check_last_state <- function(fit, predicted_rect) {
  last_state_time <- function(x) {
    max(which(!is.na(x) & x > 0))
  }

  last_state_times <- fit$data_processed$standata$obs_states_rect %>% apply(MARGIN = 1, last_state_time)
  select_matrix <- array(NA, c(length(last_state_times), 2))
  select_matrix[,1] <- 1:length(last_state_times)
  select_matrix[,2] <- last_state_times

  last_states_observed <- fit$data_processed$standata$obs_states_rect[select_matrix]

  last_states <- apply(predicted_rect, MARGIN = 3, FUN = function(x) { x[select_matrix[, 2:1]] })

  bayesplot::ppc_bars(last_states_observed, t(last_states))
}


pp_check_state_at <- function(fit, predicted_rect, day, newdata = NULL) {
  last_state_time <- function(x) {
    max(which(!is.na(x) & x > 0))
  }

  if(is.null(newdata)) {
    pred_rawdata <- fit$data
  } else {
    pred_rawdata <- newdata
  }

  data_hmm <- make_data_hmm(pred_rawdata)

  last_state_times <- data_hmm$standata$obs_states_rect %>% apply(MARGIN = 1, last_state_time)

  use_rows <- last_state_times >= day

  states_at_day_observed <- data_hmm$standata$obs_states_rect[use_rows, day]

  states_at <- predicted_rect[day, use_rows, ]

  bayesplot::ppc_bars(states_at_day_observed, t(states_at))
}


pp_check_transitions <- function(fit, predicted_rect, states_from = NULL, states_to = NULL,
                                 binwidth = NULL, scale = "prob", ...) {
  N_states_observed <- fit$data_processed$standata$N_states_observed
  if(is.null(states_from)) {
    states_from <- 1:N_states_observed
  }
  if(is.null(states_to)) {
    states_to <- 1:N_states_observed
  }

  obs_states_rect_t <- t(fit$data_processed$standata$obs_states_rect)
  include <- obs_states_rect_t != 0

  n_transitions_func <- function(rect, from, to) {
    n_time <- dim(rect)[1]
    from_match <- include[1:(n_time - 1),] & rect[1:(n_time - 1),] == from
    to_match <- include[2:n_time,] & rect[2:n_time,] == to
    sum(from_match & to_match)
  }
  prob_transition_func <- function(rect, from, to) {
    n_time <- dim(rect)[1]
    from_match <- include[1:(n_time - 1),] & rect[1:(n_time - 1),] == from
    to_match <- include[2:n_time,] & rect[2:n_time,] == to
    sum(from_match & to_match) / sum(from_match & include[2:n_time,])
  }

  if(scale == "prob") {
    summary_func <- prob_transition_func
    extra <- expand_limits(x = c(0,1))
    if(is.null(binwidth)) {
      binwidth = 0.02
    }
  } else if(scale == "count") {
    summary_func <- n_transition_func
    extra <- NULL
    if(is.null(binwidth)) {
      binwidth = 1
    }
  } else {
    stop("Unrecognized scale")
  }

  N_items <- length(states_from) * length(states_to)
  y <- array(NA_integer_, N_items)
  yrep <- array(NA_integer_, c(dim(predicted_rect)[3], N_items))
  group <- array(NA_character_, N_items)
  index <- 1
  for(from in states_from) {
    for(to in states_to) {

      y[index] <- summary_func(obs_states_rect_t, from, to)

      yrep[, index] <- apply(predicted_rect, MARGIN = 3, FUN = summary_func, from = from, to = to)

      group[index] <- paste0(from, " -> ", to)

      index <- index + 1
    }
  }

  group <- factor(group)

  bayesplot::ppc_stat_grouped(y, yrep, group, binwidth = binwidth) + extra + theme(axis.text.x = element_blank(), strip.text = element_text(size = 7))
}


pp_check_transitions_direction <- function(fit, predicted_rect, states_from = NULL, binwidth = 0.02) {
  N_states_observed <- fit$data_processed$standata$N_states_observed

  obs_states_rect_t <- t(fit$data_processed$standata$obs_states_rect)
  indices <- obs_states_rect_t == 0

  compute_directions <- function(rect) {
    rect[indices] <- NA
    if(!is.null(states_from)) {
      rect[! (rect %in% states_from) ] <- NA
    }
    n_time <- dim(rect)[1]
    diffs <- apply(rect, MARGIN = 1, FUN = diff)
    diffs <- as.vector(diffs[!is.na(diffs)])
    res <- c(mean(diffs < 0), mean(diffs == 0), mean(diffs > 0))
  }

  y <- compute_directions(obs_states_rect_t)
  yrep <- apply(predicted_rect, MARGIN = 3, FUN = compute_directions)

  bayesplot::ppc_stat_grouped(y, t(yrep), group = c("< 0", "0", "> 0"), binwidth = binwidth)
}

do_common_pp_checks <- function(fit) {
  epred_rect <- posterior_epred_rect(fit)
  predicted_rect <- posterior_epred_to_predicted(fit, epred_rect)
  predicted <- posterior_rect_to_long(fit, predicted_rect)

  predicted_df <- posterior_long_to_df(fit$data, predicted)

  all_ids <- fit$data$serie_data %>% select(hospital_id, .serie) %>%
    distinct() %>% arrange(hospital_id) %>% pull(.serie)
  step_size <- 6
  for(step in 1:ceiling(length(all_ids) / step_size)) {
    series_id <- ((step - 1) * step_size + 1) : (step* step_size)

    posterior_state_plot(predicted_df, all_ids[series_id], fit$data) %>% print()
  }

  pp_check_transitions_direction(fit, predicted_rect) %>% print()
  pp_check_transitions_direction(fit, predicted_rect, states_from = 2:4)  %>% print()


}
