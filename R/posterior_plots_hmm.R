model_color <- "#3993dd"
data_color <- "#91171f"
decoration_color <- "#ca5310"

posterior_state_plot <- function(predictions_wide, series_id = 1:max(predictions_wide$.serie), observed_data = NULL) {
  predictions_grouped <- predictions_wide %>%
    filter(.serie %in% series_id) %>%
    group_by_at(setdiff(names(predictions_wide), ".sample")) %>%
    summarise(n_samples = n()) %>%
    group_by_at(setdiff(names(predictions_wide), c(".predicted", ".sample"))) %>%
    mutate(state_prob = n_samples / sum(n_samples)) %>%
    ungroup()

  if(!is.null(observed_data)) {
    observed_data <- validate_brmshmmdata(observed_data)
    observed_data <- observed_data$serie_data %>% filter(.serie %in% series_id)
    observed_geom <- geom_line(data = observed_data, aes(x = .time, y = .observed), color = data_color, inherit.aes = FALSE, size = 2)
  } else {
    observed_geom <- NULL
  }

  print(predictions_grouped %>%
          ggplot(aes(x = .time, y = .predicted, fill = state_prob)) + geom_raster() +
          observed_geom +
          scale_fill_gradient(low = "white", high = model_color, limits = c(0,1)) +
          facet_wrap(~.serie)
  )
}
