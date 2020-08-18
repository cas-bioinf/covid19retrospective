#Throughout we use global variables from data_mapping.R

read_raw_data <- function(file, hospital_id, lang, file_version, remove_examples = TRUE) {
  file_version <- update_and_check_version(file, lang, file_version)
  patient_data <- read_patient_data(file, hospital_id, lang, file_version, remove_examples)

  progression_data <- read_progression_data(file, hospital_id, lang, file_version, patient_data)

  # TODO code adverse events
  loo::nlist(patient_data,
             breathing_data = progression_data$breathing_data,
             marker_data = progression_data$marker_data,
             adverse_events_data = progression_data$adverse_events_data)
}


update_and_check_version <- function(file, lang, file_version) {
  if(file_version == "From Sheet") {
    file_version <- readxl::read_excel(file, sheet = version_sheet(lang), range = "C8:C8", col_types = "text", col_names = "version")$version
  }

  if(!grepl("^(Old|1\\.[2-4]\\.[0-9])$", file_version)) {
    print(file_version)
    stop("Invalid version")
  }

  file_version
}

read_patient_data <- function(file, hospital_id, lang, file_version, remove_examples) {
  patient_data_raw <- read_patient_data_raw(file, hospital_id, lang, file_version)

  patient_data_raw_usable <- patient_data_raw %>% filter(!is.na(Age), !is.na(Sex))

  if(remove_examples) {
    patient_data_raw_usable <- patient_data_raw_usable %>%
      filter(!startsWith(`Patient ID`, "Example"), !grepl("^P.*klad", `Patient ID`))
  }


  logical_column <- function(x) {
    if(any( !is.na(x) && !(tolower(x) %in% c(yes_variants, no_variants))) ) {
      stop("Invalid logical column contents")
    }
    tolower(x) %in% yes_variants
  }

  outcome_labels <- c("Discharged","Hospitalized","Transferred","Death")

  # Special handling of badly filled outcomes
  if(hospital_id == "QKuFp") {
    patient_data_raw_usable <- patient_data_raw_usable %>%
      mutate(`Outcome (Discharged/Hospitalized/Transferred/Death)` =
               if_else(grepl("^439[0-9]{2}$", `Outcome (Discharged/Hospitalized/Transferred/Death)`),
                       "Hospitalized",
                       `Outcome (Discharged/Hospitalized/Transferred/Death)`))
  }


  patient_data_typed <- patient_data_raw_usable %>%
    transmute(
      hospital_id = hospital_id,
      patient_id = `Patient ID`,
      age = Age,
      sex = factor(Sex, levels = c("M","F")),
      symptom_onset = lubridate::as_date(`Date of symptom onset`),
      admission = lubridate::as_date(`Date of admission`),
      transferred_from = `Transferred from`,
      admitted_for_covid = logical_column(`Admitted for Covid`),
      ischemic_heart_disease = logical_column(`Ischemic Heart Disease`),
      has_hypertension_drugs = !(tolower(`Hypertension drugs`) %in% no_variants),
      n_hypertension_drugs = if_else(has_hypertension_drugs,
                                     as.integer(if_else(has_hypertension_drugs, `Hypertension drugs`, NA_character_)), #Inner if_else to avoid warning
                                     0L),
      heart_failure = logical_column(`Heart failure`),
      COPD = logical_column(COPD),
      asthma = logical_column(`Asthma Bronchiale`),
      other_lung_disease = logical_column(`Other lung disease`),
      diabetes = logical_column(Diabetes),
      renal_disease = logical_column(`Renal Disease`),
      liver_disease = logical_column(`Liver Disease`),
      smoking = logical_column(Smoking),
      BMI = BMI,
      NYHA = NYHA,
      creatinin = `Creatinin [μmol/L]`,
      pt_inr = `PT INR (Quick)`,
      albumin = `Albumin [g/L]`,
      discontinued_medication = !(tolower(`Discontinued medication`) %in% no_variants),
      discontinued_medication_reason = if_else(discontinued_medication, `Discontinued medication`, NA_character_),
      best_supportive_care_from = if_else(tolower(`Best supportive care from:`) %in% no_variants, lubridate::as_date(NA),
                                          janitor::excel_numeric_to_date(as.numeric(
                                            if_else(tolower(`Best supportive care from:`) %in% no_variants, NA_character_, `Best supportive care from:`) # Inner if_else to avoid warning
                                          ))),
      last_record = lubridate::as_date(`Date of last record`),
      outcome = factor(tolower(`Outcome (Discharged/Hospitalized/Transferred/Death)`), levels = tolower(outcome_labels), labels = outcome_labels)
    )


  invalid_outcome <- is.na(patient_data_raw_usable$`Outcome (Discharged/Hospitalized/Transferred/Death)`) != is.na(patient_data_typed$outcome)
  if(any(invalid_outcome)) {
    print(patient_data_raw_usable %>% filter(invalid_outcome))
    stop("Invalid outcome")
  }

  patient_data_typed
}

read_patient_data_raw <- function(file, hospital_id, lang, file_version) {
  if(file_version == "Old") {
    n_middle_text_columns <- 12
    last_column <- "Z"
  } else {
    n_middle_text_columns <- 13
    last_column <- "AA"
  }

  patient_col_types <- c("text","numeric","text","date","date", #ID - Date of admission
                         rep("text", n_middle_text_columns), #Transferred from - Smoking
                         rep("numeric", 5), # BMI - Albumin
                         "text", "text", "date", "text") #Discontinued medication - Outcome

  patients_range <- paste0("A2:",last_column, "100")


  patient_data_raw <- readxl::read_excel(file, sheet = patients_sheet(lang), range = patients_range, na = "NA", col_types = patient_col_types) %>%
    translate_patient_columns(lang)

  patient_data_raw <- patient_data_raw %>%
    mutate(`Outcome (Discharged/Hospitalized/Transferred/Death)` = translate_outcomes(`Outcome (Discharged/Hospitalized/Transferred/Death)`, lang))

  if(file_version == "Old") {
    patient_data_raw <- patient_data_raw %>% mutate(`Covid Medications before admission` = NA_character_)
  }

  patient_data_raw
}


read_progression_data <- function(file, hospital_id, lang, file_version, patient_data_typed) {
  progression_data_raw <- readxl::read_excel(file, sheet = progression_sheet(lang), range = "A1:CU2539", na = c("NA",""), col_types = c(c("text","text", "skip", "text", "text"), rep("text", 94)), col_names = TRUE)

  progression_data <- progression_data_raw %>%
    translate_progression_columns(lang) %>%
    rename(patient_id = `Patient ID`, indicator = Indicator, first_day = First_Day, unit = Unit) %>%
    filter(!is.na(indicator)) %>%
    mutate(indicator = translate_markers(indicator, lang),
           hospital_id = hospital_id)

  # Remove the "Don't edit examples" comment
  progression_data$patient_id[2:9] <- NA_character_

  # Spread Patient ID and Start date
  last_patient_id <- progression_data$patient_id[1]
  last_first_day <- progression_data$first_day[1]
  for(i in 2:nrow(progression_data)) {
    if(!is.na(progression_data$patient_id[i])) {
      last_patient_id <- progression_data$patient_id[i]
      last_first_day <- progression_data$first_day[i]
    } else {
      progression_data$patient_id[i] <- last_patient_id
      progression_data$first_day[i] <- last_first_day
    }
  }

  progression_data <- progression_data %>% filter(patient_id %in% patient_data_typed$patient_id)

  progression_data <- progression_data %>%
    pivot_longer(starts_with("Day_"), names_to = "day", values_to = "value", names_prefix = "Day_") %>%
    mutate(day = as.integer(day))

  progression_data <- progression_data %>%
    filter(!is.na(indicator))


  # Check dates
  dates_data <- progression_data %>%
    filter(indicator == "Date")
  test_dates <- dates_data %>%
    mutate(value = as.integer(value), first_day = as.integer(first_day)) %>%
    filter(value != first_day + day - 1 )

  if(nrow(test_dates) > 0) {
    stop("Dates check failed")
  }


  breathing_data <- progression_data %>% filter(indicator == "Breathing") %>%
    mutate(
      # Special cases
      value = case_when(value == "Oxygen/NIPPV" ~ "NIPPV",
                        TRUE ~ value),
      # Translate to factor
      breathing = factor(tolower(value), levels = tolower(c("AA","Oxygen", "NIPPV","MV","ECMO")),labels = c("AA","Oxygen", "NIPPV","MV","ECMO"), ordered = TRUE))

  breathing_unrecognized <- breathing_data %>% filter(is.na(breathing) != is.na(value))
  if(nrow(breathing_unrecognized) > 0) {
    print(unique(breathing_unrecognized$value))
    stop("Breathing NAs check failed")
  }

  breathing_data <- breathing_data %>% select(-value, -indicator, -unit)

  correct_first_day_admission <- function(patient_data, progression_data) {
    nrow_before <- nrow(progression_data)
    progression_data <- progression_data %>%
      inner_join(patient_data %>% select(patient_id, admission), by = "patient_id") %>%
      mutate(first_day = janitor::excel_numeric_to_date(as.numeric(first_day)),
             correction = as.integer(first_day - admission),
             day = day + correction - 1)

    if(nrow(progression_data) != nrow_before) {
      stop("Invalid join")
    }

    progression_data %>%
      select(-admission, -first_day, -correction)
  }

  breathing_data <- correct_first_day_admission(patient_data_typed, breathing_data) %>% filter(!is.na(breathing))

  adverse_events_data <- progression_data %>%
    filter(indicator == "Adverse events")


  marker_data <- progression_data %>% filter(!(indicator %in% c("Date", "Breathing", "Adverse events", dummy_markers))) %>% rename(marker = indicator)
  marker_data <- correct_first_day_admission(patient_data_typed, marker_data)

  # Correct Excel error
  marker_data <- marker_data %>%
    mutate(marker = case_when(
      grepl("^Supp\\. O[0-9]*$", marker) ~ "Supp. O2",
      grepl("^SpO[0-9]*$", marker) ~ "SpO2",
      grepl("^IL-[0-9]*$", marker) ~ "IL-6",
      TRUE ~ marker))

  if(any(is.na(marker_data$marker))) {
    stop("NA in markers")
  }


  #Check NAs (TODO)
  unrecgonized_markers <- marker_data %>% filter(!(tolower(marker) %in% tolower(names(marker_map))))

  if(nrow(unrecgonized_markers) > 0) {
    print(unique(unrecgonized_markers$marker))
    print(unrecgonized_markers)
    stop("Unrecognized markers")
  }

  names(marker_map) <- tolower(names(marker_map))

  marker_data <- marker_data %>% mutate(marker = marker_map[tolower(marker)])

  # Special handling for PCR
  pcr_data <- marker_data %>% filter(marker == "pcr_value")


  wrong_pcr <- pcr_data %>% filter(!is.na(value),!(tolower(value) %in% pcr_values_positive), !tolower(value) %in% pcr_value_negative, !grepl("^[0-9]+$",value))
  if(nrow(wrong_pcr) > 0) {
    print(wrong_pcr)
    stop("Wrong PCR")
  }

  pcr_data_long <- pcr_data %>% mutate(pcr_positive = as.numeric(tolower(value) %in% pcr_values_positive),
                                       pcr_value = as.numeric(if_else(grepl("^[0-9]+$",value), value, NA_character_))) %>%
    select(-value, -marker) %>%
    pivot_longer(c("pcr_value","pcr_positive"), names_to = "marker", values_to = "value") %>%
    mutate(censored = "none")

  # Handle NT for SpO2_native
  marker_data <- marker_data %>%
    mutate(value = if_else(marker == "SpO2_native" & tolower(value) == "nt", "< 92", value) )

  # Some special cases
  # TODO ask at hospital
  if(hospital_id == "QKuFp") {
    marker_data <- marker_data %>% mutate(value = if_else(marker == "SpO2" & value == "NEG", NA_character_, value))
  } else if(hospital_id == "YqNbe") {
    marker_data <- marker_data %>% mutate(value = if_else(marker %in% c("hcq","az") & value == "ano", ">100", value))
  }

  # Handle censoring on everything except PCR
  marker_data <- marker_data %>% filter(marker != "pcr_value") %>%
    mutate(censored = case_when(startsWith(value, "<") ~ "right",
                                startsWith(value, ">") ~ "left",
                                TRUE ~ "none"),
           value_parsed = as.numeric(gsub("^[<>]? *", "", value)))


  bad_parse_value <- marker_data %>% filter(is.na(value_parsed) & !is.na(value))
  if(nrow(bad_parse_value) > 0) {
    print(bad_parse_value)
    stop("Bad parse value")
  }


  # Join with PCR data
  marker_data <- marker_data %>%
    mutate(value = value_parsed) %>%
    select(-value_parsed) %>%
    rbind(pcr_data_long)

  # Some special handling
  if(hospital_id %in% names(unit_overwrites)) {
    for(overwrite in unit_overwrites[[hospital_id]]) {
      marker_data <- marker_data %>%
        mutate(unit = if_else(marker == overwrite$marker & unit == overwrite$old_unit, overwrite$new_unit, unit))
    }
  }

  loo::nlist(marker_data, breathing_data, adverse_events_data)
}

anonymize_for_analysis <- function(data) {
  # TODO

  # Create new more anonymized IDs
  patient_data_typed <- patient_data_typed %>%
    mutate(patient_id_rnd = paste0(hospital_id, sample(1:nrow(patient_data_typed), replace = FALSE)))

  if(length(unique(patient_data_typed$patient_id_rnd)) != nrow(patient_data_typed)) {
    stop("Bad Random IDs")
  }

}