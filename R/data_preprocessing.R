read_data_for_analysis <- function() {
  patient_data <- read_csv(here::here("private_data", "patient_data.csv"),
                           col_types = cols(
                             hospital_id = col_character(),
                             patient_id = col_character(),
                             age = col_double(),
                             sex = col_character(),
                             days_from_symptom_onset = col_integer(),
                             admitted_for_covid = col_logical(),
                             ischemic_heart_disease = col_logical(),
                             has_hypertension_drugs = col_logical(),
                             n_hypertension_drugs = col_integer(),
                             heart_failure = col_logical(),
                             COPD = col_logical(),
                             asthma = col_logical(),
                             other_lung_disease = col_logical(),
                             diabetes = col_logical(),
                             renal_disease = col_logical(),
                             liver_disease = col_logical(),
                             smoking = col_logical(),
                             BMI = col_double(),
                             obesity = col_character(),
                             NYHA = col_factor(levels = c(1,2,3,4), ordered = TRUE),
                             creatinin = col_double(),
                             pt_inr = col_double(),
                             albumin = col_double(),
                             discontinued_medication = col_logical(),
                             best_supportive_care_from = col_integer(),
                             last_record = col_integer(),
                             first_wave = col_logical(),
                             outcome = col_factor()
                           )) %>%
    mutate(age_norm = (age - 65) / 20,
           high_creatinin = if_else(sex == "M", creatinin > 115, creatinin > 97),
           high_pt_inr = pt_inr > 1.2,
           low_albumin = albumin < 36,
           obesity = factor(case_when(is.na(BMI) ~ obesity,
                               BMI <= 25 ~ "No",
                               BMI <= 30 ~ "Overweight",
                               TRUE ~ "Obese")
                            ,levels = c("No", "Overweight", "Obese"), ordered = TRUE)
           )

  collapse_breathing <- function(x) {
    fct_collapse(x, Ventilated = c("NIPPV", "MV", "ECMO"))
  }
  breathing_data <- read_csv(here::here("private_data", "breathing_data.csv"),
                           col_types = cols(
                            patient_id = col_character(),
                            hospital_id = col_character(),
                            day = col_integer(),
                            breathing = col_factor(levels = disease_levels, ordered = TRUE),
                            breathing_low = col_factor(levels = disease_levels, ordered = TRUE),
                            breathing_high = col_factor(levels = disease_levels, ordered = TRUE)
                           )) %>%
    mutate(breathing_s = collapse_breathing(breathing),
           breathing_s_low = collapse_breathing(breathing_low),
           breathing_s_high = collapse_breathing(breathing_low))

  marker_data <- read_csv(here::here("private_data", "marker_data.csv"),
                             col_types = cols(
                              patient_id = col_character(),
                              marker = col_character(),
                              unit = col_character(),
                              hospital_id = col_character(),
                              day = col_integer(),
                              value = col_double(),
                              censored = col_character()
                            )) %>%
    # One more bug quick fix (fixed in raw data for further imports)
    mutate(marker = if_else(patient_id == "809N695Ho6" & marker == "SpO2_native" & day >= 4, "Horowitz_index", marker)) %>%
    # Filtering harmless NAs (due to error in data)
    filter(!(patient_id == "d2TgUar00l" & is.na(value)))




  multiple_units <- marker_data %>%
    group_by(marker) %>%
    summarise(n_units = length(unique(unit))) %>%
    filter(n_units > 1)
  if(nrow(multiple_units) > 1) {
    print(multiple_units)
    stop("Multiple units for some markers")
  }


  data <- loo::nlist(patient_data, breathing_data, marker_data)
  data$marker_data_wide <- prepare_marker_data_wide(data)

  data <- compute_derived_quantities_patients(data)

  # TODO Reconstruct long from wide (because wide has more derived quantities and covers all days)
  # Would require keeping units around in marker_data_wide


  # Remove unsed xxx_censored fields from wide
  censored_columns <- names(data$marker_data_wide)[grepl("_censored$", names(data$marker_data_wide))]
  for(censored_col in censored_columns) {
    if(all(is.na(data$marker_data_wide[[censored_col]]) | data$marker_data_wide[[censored_col]] == "none")) {
      data$marker_data_wide[[censored_col]] <- NULL
    }
  }

  data
}

#' @importFrom tidyr replace_na
compute_derived_quantities_patients <- function(data) {
  derived_from_wide <- data$marker_data_wide %>%
    group_by(patient_id) %>%
    summarise(ever_hcq = any(took_hcq),
              ever_az = any(took_az),
              ever_convalescent_plasma = any(took_convalescent_plasma),
              ever_favipiravir = any(took_favipiravir),
              ever_dexamethasone = any(took_dexamethasone),
              ever_remdesivir = any(took_remdesivir),
              ever_antibiotics = any(took_antibiotics),
              any_d_dimer = any(!is.na(d_dimer)),
              any_IL_6 = any(!is.na(IL_6)),
              worst_breathing = max(breathing, na.rm = TRUE),
              worst_breathing_s = max(breathing_s, na.rm = TRUE),
              had_invasive = any(breathing >= "MV" & breathing < "Death"),
              first_day_invasive = if_else(had_invasive, min(c(day[breathing >= "MV"], 1e6)), NA_real_),
              last_day_invasive = if_else(had_invasive, max(c(day[breathing >= "MV"], 0)), NA_real_),
              .groups = "drop"
    )

  data$patient_data <- data$patient_data %>%
    mutate(
      is_obese = obesity == "Obese",
      heart_problems = NYHA > 1 | ischemic_heart_disease | heart_failure,
      lung_problems = COPD | asthma | other_lung_disease,
      markers_out_of_range = high_creatinin | high_pt_inr | low_albumin,
      comorbidities_sum =
        replace_na(heart_problems, FALSE) +
        replace_na(has_hypertension_drugs, FALSE) +
        replace_na(lung_problems, FALSE) +
        replace_na(diabetes, FALSE) +
        replace_na(renal_disease, FALSE) +
        replace_na(liver_disease, FALSE) +
        replace_na(smoking, FALSE) +
        replace_na(markers_out_of_range, FALSE) +
        replace_na(is_obese, FALSE),
      comorbidities_sum_na = 2 * (
        replace_na(heart_problems, 0.5) +
          replace_na(has_hypertension_drugs, 0.5) +
          replace_na(lung_problems, 0.5) +
          replace_na(diabetes, 0.5) +
          replace_na(renal_disease, 0.5) +
          replace_na(liver_disease, 0.5) +
          replace_na(smoking, 0.5) +
          replace_na(markers_out_of_range, 0.5) +
          replace_na(is_obese, 0.5))
      )

  data$patient_data <- data$patient_data %>%
    left_join(derived_from_wide, by = "patient_id")


  data
}

prepare_marker_data_wide <- function(data) {
  patient_ranges <- data$breathing_data %>%
    select(hospital_id, patient_id, day) %>%
    rbind(data$marker_data %>% select(hospital_id, patient_id, day)) %>%
    group_by(hospital_id, patient_id) %>%
    summarise(max_day = max(day), min_day = min(day), .groups = "drop")
  patient_ranges2 <- data$marker_data %>%
    group_by(hospital_id, patient_id) %>%
    summarise(max_day = max(day), min_day = min(day), .groups = "drop")


  complete_max_day <- max(patient_ranges$max_day)
  complete_min_day <- min(patient_ranges$min_day)

  wider_spec_censored <-
    tibble(
      marker = data$marker_data %>% pull(marker) %>% unique()
    ) %>%
    mutate(.name = paste0(marker, "_censored"), .value = "censored")

  wider_spec <-
    build_wider_spec(data$marker_data, names_from = "marker", values_from = "value") %>%
    rbind(wider_spec_censored)
  marker_wide <- data$marker_data %>%
    select(-unit) %>% pivot_wider_spec(wider_spec, id_cols = one_of("patient_id", "hospital_id", "day"), values_fn = function(x) { if(length(x) > 1) { "ERROR" } else { x } })

  which((marker_wide %>% select(oxygen_flow:az_censored) %>% as.matrix()) > 1, arr.ind = TRUE)

  all_values <- marker_wide %>% select(oxygen_flow:az_censored) %>% unlist()
  if(any(grepl("ERROR", all_values))) {
    stop("Error in pivot")
  }

  all_patient_days <- patient_ranges %>%
    crossing(day = complete_min_day:complete_max_day) %>%
    filter(day <= max_day, day >= min_day) %>%
    select(-max_day, -min_day)

  markers_wide <- all_patient_days %>%
    left_join(data$breathing_data, by = c("patient_id", "hospital_id", "day")) %>%
    left_join(marker_wide, by = c("patient_id", "hospital_id", "day"))


  compute_derived_markers_wide(markers_wide, data)
}


compute_derived_markers_wide <- function(markers_wide, data) {
  nrow_before <- nrow(markers_wide)

  any_antibiotics <-  data$marker_data %>%
    filter(marker %in% setdiff(all_antibiotics, antibiotics_not_for_pneumonia)) %>%
    group_by(patient_id, day) %>%
    summarise(antibiotics = any(!is.na(value) & value > 0))

  any_macrolides <-  data$marker_data %>%
    filter(marker %in% all_macrolides) %>%
    group_by(patient_id, day) %>%
    summarise(macrolides = any(!is.na(value) & value > 0))

  markers_wide <- markers_wide %>%
    left_join(any_antibiotics, by = c("patient_id", "day")) %>%
    left_join(any_macrolides, by = c("patient_id", "day"))

  treatment_markers <- c(quo(hcq), quo(az),
                         quo(dexamethasone), quo(remdesivir),
                         quo(convalescent_plasma),
                         quo(favipiravir),
                         quo(antibiotics), quo(macrolides))
  for(treatment_mark in treatment_markers) {
    first_treatment <- markers_wide %>% filter(!is.na(!!treatment_mark) & !!treatment_mark > 0) %>%
      group_by(hospital_id, patient_id) %>%
      summarise(first_day = min(day))

    new_col_name <- paste0("took_", rlang::as_name(treatment_mark))
    markers_wide <- markers_wide %>%
      left_join(first_treatment, by = c("hospital_id", "patient_id")) %>%
      mutate(!!new_col_name := !is.na(first_day) & first_day <= day) %>%
      select(-first_day)

    if(any(is.na(markers_wide[[new_col_name]]))) {
      stop(paste0("NA in ", new_col_name))
    }
  }

  if(nrow(markers_wide) != nrow_before) {
    stop("Failed processing")
  }

  # Markers from patient data
  markers_wide <- markers_wide %>%
    inner_join(data$patient_data %>% select(hospital_id, patient_id, best_supportive_care_from),
               by = c("hospital_id", "patient_id")) %>%
    mutate(best_supportive_care = !is.na(best_supportive_care_from) & day >= best_supportive_care_from) %>%
    select(-best_supportive_care_from)


  if(nrow(markers_wide) != nrow_before) {
    stop("Failed processing")
  }
  markers_wide
}

impute_marker_constant <- function(markers_wide, column_name, initial_value = NA) {
  markers_wide <- markers_wide %>% arrange(patient_id, day)

  last_patient_id <- ""
  last_value <- NA
  for(i in 1:nrow(markers_wide)) {
    if(markers_wide$patient_id[i] != last_patient_id) {
      last_value = initial_value
    }
    if(is.na(markers_wide[[column_name]][i])) {
      markers_wide[[column_name]][i] <- last_value
    } else {
      last_value <- markers_wide[[column_name]][i]
    }
  }

  markers_wide
}

compute_marker_peak <- function(markers_wide, column_name, new_column_name, initial_value = NA) {
  markers_wide <- markers_wide %>% arrange(patient_id, day)

  last_patient_id <- ""
  current_max <- NA
  peaks <- numeric(nrow(markers_wide))
  for(i in 1:nrow(markers_wide)) {
    if(markers_wide$patient_id[i] != last_patient_id) {
      current_max = initial_value
    }
    if(!is.na(markers_wide[[column_name]][i])) {
      current_max = max(current_max, markers_wide[[column_name]][i])
    }
    peaks[i] <- current_max
  }

  markers_wide[[new_column_name]] <- peaks

  markers_wide
}

get_data_version <- function() {
  readLines(here::here("private_data", "data_revision.txt"))
}
