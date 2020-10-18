#' @param  optimize_possible Let's the algorithm optimize for models where some
#' hidden states are impossible, given an observation. The optimization
#' some overhead and so can be turned of. Takes values 0 - no such optimalization),
#'  1 - avoid redundant rows in matrix to multiply and 2 - also avoid some multiplications.
brmshmmdata <- function(formula,
                        serie_data, rate_data,
                        hidden_state_data, initial_states,
                        observed_state_data = NULL, prior = NULL,
                        sensitivity_low_bound = 0.5, optimize_possible = 2) {

  validate_brmshmmdata(
    structure(loo::nlist(
      formula = make_brms_formula_hmm(formula),
      serie_data,
      rate_data,
      hidden_state_data,
      initial_states,
      observed_state_data,
      prior,
      sensitivity_low_bound,
      optimize_possible),
      class = "brmshmmdata"
      ))
}

validate_id <- function(id, name_for_message,
                        force_unique = FALSE, force_no_gaps = FALSE, force_no_na = TRUE,
                        reference_levels = NULL, reference_levels_for_message = NULL) {
  if(force_no_na && any(is.na(id))) {
    stop(paste0(name_for_message, " contains NA values"))
  }
  if(!is.factor(id)) {
    if(is.double(id)) {
      if(!identical(round(id), id)) {
        stop(paste0(name_for_message, " includes non-integer numbers"))
      }
      id <- as.integer(id)
    }
    if(is.integer(id) || is.character(id)) {
      if(is.null(reference_levels)) {
        id <- factor(id, levels = sort(unique(id)))
      } else {
        na_orig <- is.na(id)
        id <- factor(id, levels = reference_levels)
        if(any(is.na(id) & !na_orig)) {
          stop(paste0("Some components of ", name_for_message, " do not fit into reference levels (", reference_levels_for_message, " )."))
        }
      }
    } else {
      stop(paste0(name_for_message, " must be a factor, integer or character"))
    }
  } else if(!is.null(reference_levels) && !all(levels(id) == reference_levels)) {
    stop(paste0("Levels of ", name_for_message, " are different than reference levels (", reference_levels_for_message, " )."))
  }

  if(force_unique && !identical(sort(unique(as.integer(id))), 1:length(id))) {
    stop(paste0(name_for_message, " has to be unique and without gaps"))
  } else if(force_no_gaps && !identical(sort(unique(as.integer(id))), 1:max(as.integer(id)))) {
    stop(paste0(name_for_message, " has to be without gaps"))
  }

  id
}

validate_brmshmmdata <- function(d) {
  if(!inherits(d, "brmshmmdata")) {
    stop("Must be of class brmshmmdata")
  }

  d$rate_data <- dplyr::ungroup(d$rate_data)
  d$serie_data <- dplyr::ungroup(d$serie_data)
  d$hidden_state_data <- dplyr::ungroup(d$hidden_state_data)
  d$observed_state_data <- dplyr::ungroup(d$observed_state_data)

  if(is.null(d$rate_data$.rate_id)) {
    d$rate_data <- d$rate_data %>% mutate(.rate_id = factor(1:n()))
  } else {
    d$rate_data$.rate_id <- validate_id(d$rate_data$.rate_id, "rate_data$.rate_id", force_no_gaps = TRUE)
  }

  #TODO: test for duplicate (rate_from, rate_to)
  #TODO: test for states with no rates


  d$hidden_state_data$id <- validate_id(d$hidden_state_data$id, "hidden_state_data$id", force_unique = TRUE, force_no_gaps = TRUE)
  d$hidden_state_data <- d$hidden_state_data %>% arrange(as.integer(id))

  d$rate_data$.from <- validate_id(d$rate_data$.from, "rate_data$.from",
                                   reference_levels = levels(d$hidden_state_data$id),
                                   reference_levels_for_message = "hidden_state_data$id")

  d$rate_data$.to <- validate_id(d$rate_data$.to, "rate_data$.to",
                                 reference_levels = levels(d$hidden_state_data$id),
                                 reference_levels_for_message = "hidden_state_data$id")


  if(is.null(d$observed_state_data)) {

    d$hidden_state_data$corresponding_obs <-
      validate_id(d$hidden_state_data$corresponding_obs, "corresponding_obs (without observed data)", force_no_gaps = TRUE)


    d$observed_state_data <- data.frame(id = levels(d$hidden_state_data$corresponding_obs), is_noisy = FALSE)
  } else {
    d$hidden_state_data$corresponding_obs <-
      validate_id(d$hidden_state_data$corresponding_obs, "corresponding_obs", force_no_gaps = TRUE,
                  reference_levels = levels(d$observed_state_data$id), reference_levels_for_message = "observed_state_data$id")
  }

  d$observed_state_data$id <- validate_id(d$observed_state_data$id, force_unique = TRUE, force_no_gaps = TRUE)
  d$observed_state_data <- d$observed_state_data %>% arrange(as.integer(id))

  other_obs_cols <- names(d$observed_state_data)[grepl("^other_obs_", names(d$observed_state_data))]
  for(oc in other_obs_cols) {
    d$observed_state_data[[oc]] <- validate_id(d$observed_state_data[[oc]], oc,
                                               force_no_na = FALSE,
                                               reference_levels = levels(d$observed_state_data$id),
                                               reference_levels_for_message = "observed_state_data$id")
  }


  d$initial_states <- validate_id(d$initial_states, "initial_states", force_unique = FALSE,
                                  force_no_gaps = FALSE, reference_levels = levels(d$hidden_state_data$id),
                                  reference_levels_for_message = "hidden_state_data$id")

  d$serie_data$.observed <- validate_id(d$serie_data$.observed, "serie_data$.observed",
                                        force_no_na = FALSE,
                                        reference_levels = levels(d$observed_state_data$id),
                                        reference_levels_for_message = "observed_state_data$id")

  d$serie_data$.serie <- validate_id(d$serie_data$.serie, "serie_data$.serie",
                                     force_no_gaps = TRUE)

  #TODO make the binding between initial states and series explicit
  if(length(d$initial_state) != length(unique(d$serie_data$.serie))) {
    stop("Incorrect number of initial states")
  }

  d
}
