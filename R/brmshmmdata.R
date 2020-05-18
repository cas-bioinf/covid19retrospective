brmshmmdata <- function(formula,
                        serie_data, rate_data,
                        hidden_state_data, initial_states,
                        observed_state_data = NULL, prior = NULL,
                        sensitivity_low_bound = 0.5) {

  validate_brmshmmdata(
    structure(loo::nlist(
      formula = make_brms_formula_hmm(formula),
      serie_data,
      rate_data,
      hidden_state_data,
      initial_states,
      observed_state_data,
      prior,
      sensitivity_low_bound),
      class = "brmshmmdata"
      ))
}

validate_brmshmmdata <- function(d) {
  if(!inherits(d, "brmshmmdata")) {
    stop("Must be of class brmshmmdata")
  }

  if(is.null(d$rate_data$.rate_id)) {
    d$rate_data <- d$rate_data %>% mutate(.rate_id = factor(1:n()))
  } else {
    if(!is.factor(d$rate_data$.rate_id)) {
      stop("rate_data$.rate_id must be a factor")
    }
  }

  #TODO: test for duplicate rates
  #TODO: test for states with no rates

  d$hidden_state_data <- d$hidden_state_data %>% arrange(id)
  if(!identical(d$hidden_state_data$id, 1:nrow(d$hidden_state_data))) {
    stop("Hiden state ids need to be consecutive and not duplicated and span the same range as rate_data")
  }

  all_corresponding_states <- d$hidden_state_data %>% pull(corresponding_obs) %>%
    unique() %>% sort()


  if(is.null(d$observed_state_data)) {

    if(!identical(all_corresponding_states, 1:length(all_corresponding_states))) {
      stop("If observed_state_data is NULL, corresponding_obs cannot have gaps")
    }

    d$observed_state_data <- data.frame(id = all_corresponding_states, is_noisy = FALSE)
  }

  d$observed_state_data <- d$observed_state_data %>% arrange(id)
  if(!identical(d$observed_state_data$id, 1:nrow(d$observed_state_data))) {
    stop("Observed state ids need to be consecutive.")
  }

  if(!all(all_corresponding_states %in% d$observed_state_data$id)) {
    stop("Some corresponding observations are not in observed_state_data")
  }

  d
}
