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
                             NYHA = col_factor(levels = c(1,2,3,4), ordered = TRUE),
                             creatinin = col_double(),
                             pt_inr = col_double(),
                             albumin = col_double(),
                             discontinued_medication = col_logical(),
                             best_supportive_care_from = col_integer(),
                             last_record = col_integer(),
                             outcome = col_factor()
                           ))

  breathing_data <- read_csv(here::here("private_data", "breathing_data.csv"),
                           col_types = cols(
                            patient_id = col_character(),
                            hospital_id = col_character(),
                            day = col_integer(),
                            breathing = col_factor(levels = breathing_levels, ordered = TRUE)
                          ))

  marker_data <- read_csv(here::here("private_data", "marker_data.csv"),
                             col_types = cols(
                              patient_id = col_character(),
                              marker = col_character(),
                              unit = col_character(),
                              hospital_id = col_character(),
                              day = col_integer(),
                              value = col_double(),
                              censored = col_character()
                            ))

  loo::nlist(patient_data, breathing_data, marker_data)
}
