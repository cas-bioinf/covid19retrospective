#Throughout we use global variables from data_mapping.R

read_raw_data <- function(file, hospital_id, lang, file_version, remove_examples = TRUE) {
  file_version <- update_and_check_version(file, lang, file_version)
  patient_data <- read_patient_data(file, hospital_id, lang, file_version, remove_examples)

  progression_data <- read_progression_data(file, hospital_id, lang, file_version, patient_data)

  if(file_version == "Special" && any(progression_data$marker_data$marker == "weight")) {
    #Use Height/Weight from lab markers to calculate BMI
    BMI_lab_data <-  progression_data$marker_data %>%
      select(patient_id, marker, value, day) %>%
      filter(marker %in% c("height", "weight")) %>%
      pivot_wider(names_from = "marker", values_from = "value") %>%
      group_by(patient_id) %>%
      filter(day == min(day)) %>%
      ungroup() %>%
      mutate(BMI_lab = weight / (( height /100) ^ 2)) %>%
      select(patient_id, BMI_lab)

    nrow_before <- nrow(patient_data)

    patient_data <- patient_data %>%
      left_join(BMI_lab_data, by = "patient_id") %>%
      mutate(BMI = if_else(!is.na(BMI_lab), BMI_lab, BMI)) %>%
      select(-BMI_lab)

    if(nrow(patient_data) != nrow_before) {
      stop("Bad join")
    }

    progression_data$marker_data <-
      progression_data$marker_data %>% filter(!(marker %in% c("height", "weight")))

    # Use markers from first two days to impute markers at admission
    markers_lab_data <- progression_data$marker_data %>%
      select(patient_id, marker, value, day) %>%
      filter(day >= 0, day <= 1, marker %in% c("pt_inr", "albumin", "creatinin")) %>%
      pivot_wider(names_from = "marker", values_from = "value") %>%
      group_by(patient_id) %>%
      filter(day == min(day)) %>%
      select(-day)

    # Make sure all columns exist
    for(m in c("pt_inr", "albumin", "creatinin")) {
      if(is.null(markers_lab_data[[m]])) {
        markers_lab_data[[m]] <- array(NA_real_, nrow(markers_lab_data))
      }
    }

    patient_data <- patient_data %>%
      select(-pt_inr, -albumin, -creatinin) %>% # I now those were not filled
      left_join(markers_lab_data, by = "patient_id")

    if(nrow(patient_data) != nrow_before) {
      stop("Bad join")
    }
  }


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
  if(file_version == "From Sheet C7") {
    file_version <- readxl::read_excel(file, sheet = version_sheet(lang), range = "C7:C7", col_types = "text", col_names = "version")$version
  }

  if(!grepl("^(Old|Special|1\\.[2-4]\\.[0-9])$", file_version)) {
    print(file_version)
    stop("Invalid version")
  }

  file_version
}

read_patient_data <- function(file, hospital_id, lang, file_version, remove_examples) {
  patient_data_raw <- read_patient_data_raw(file, hospital_id, lang, file_version)

  patient_data_raw_usable <- patient_data_raw %>% filter(!is.na(`Patient ID`), !is.na(Age), !is.na(Sex))

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

  if(hospital_id == "ZDYlY") {
    patient_data_raw_usable <- patient_data_raw_usable %>%
      mutate(`Admitted for Covid` = if_else(`Admitted for Covid` == "NA/ne?", "No", `Admitted for Covid`),
             Smoking = if_else(Smoking == "Exkurak", "no", Smoking),
             `Best supportive care from:` = if_else(grepl("ano.*hospitalizace, datum NA", patient_data_raw_usable$`Best supportive care from:`), "43920" ,`Best supportive care from:`)
             )
  }

  if(file_version == "Special") {
    BMI_raw <- patient_data_raw_usable$BMI
    BMI_new <- array(NA_real_, length(BMI_raw))

    height_weight_match <- grepl("[0-9]*cm/[0-9]*kg", BMI_raw)
    BMI_match <- BMI_raw[height_weight_match]
    BMI_split <- strsplit(BMI_match, "cm/", fixed = TRUE)
    BMI_new[height_weight_match] <- BMI_split %>%
      map_dbl( ~ as.numeric(gsub("kg", "", .x[2])) / ((as.numeric(.x[1]) / 100) ^ 2))

    special_match <- BMI_raw == "105kg"
    BMI_new[special_match] <- NA_real_

    patient_data_raw_usable$Obesity[special_match] <- "Overweight"
    patient_data_raw_usable$BMI <- BMI_new
    patient_data_raw_usable$Obesity[is.na(patient_data_raw_usable$Obesity)] <- "No"
  }

  patient_data_typed <- patient_data_raw_usable %>%
    transmute(
      hospital_id = hospital_id,
      patient_id = `Patient ID`,
      age = Age,
      sex = factor(Sex, levels = c("M","F")),
      symptom_onset = lubridate::as_date(`Date of symptom onset`),
      admission = lubridate::as_date(`Date of admission`),
      days_from_symptom_onset = as.double(admission - symptom_onset, unit = "days"),
      transferred_from = `Transferred from`,
      transferred_to = `Transferred to`,
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
      obesity = Obesity,
      NYHA = NYHA,
      creatinin = `Creatinin [μmol/L]`,
      pt_inr = `PT INR (Quick)`,
      albumin = `Albumin [g/L]`,
      discontinued_medication = !(tolower(`Discontinued medication`) %in% no_variants),
      discontinued_medication_reason = if_else(discontinued_medication, `Discontinued medication`, NA_character_),
      best_supportive_care_from_date = if_else(tolower(`Best supportive care from:`) %in% no_variants, lubridate::as_date(NA),
                                          janitor::excel_numeric_to_date(as.numeric(
                                            if_else(tolower(`Best supportive care from:`) %in% no_variants, NA_character_, `Best supportive care from:`) # Inner if_else to avoid warning
                                          ))),
      best_supportive_care_from = as.double(best_supportive_care_from_date - admission, unit = "days"),
      last_record_date = lubridate::as_date(`Date of last record`),
      last_record = as.double(last_record_date - admission, unit = "days"),
      outcome = factor(tolower(`Outcome (Discharged/Hospitalized/Transferred/Death)`), levels = tolower(outcome_labels), labels = outcome_labels)
    )


  invalid_outcome <- is.na(patient_data_raw_usable$`Outcome (Discharged/Hospitalized/Transferred/Death)`) != is.na(patient_data_typed$outcome)
  if(any(invalid_outcome)) {
    print(patient_data_raw_usable %>% filter(invalid_outcome))
    stop("Invalid outcome")
  }

  should_not_be_NA <- c("age", "sex", "outcome","patient_id", "last_record")
  for(no_NAs in should_not_be_NA) {
    if(any(is.na(patient_data_typed[[no_NAs]]))) {
      stop(paste0("NAs in ", no_NAs))
    }
  }

  invalid_dates <- patient_data_typed %>% filter(admission > last_record_date)
  if(nrow(invalid_dates) > 0) {
    print(invalid_dates)
    stop("Invalid dates")
  }

  patient_data_typed
}

read_patient_data_raw <- function(file, hospital_id, lang, file_version) {
  if(file_version == "Special") {
    patient_col_types <- c("text","numeric","text","date","date", #ID - Date of admission
                           rep("text", 15), #Transferred from - Smoking, BMI, Obesity
                           "numeric", #NYHA
                           "text", "text", "date", "text", "text") #Discontinued medication - Transferred to
    last_column <- "Z"
  } else {
    if(file_version == "Old") {
      n_middle_text_columns <- 12
      last_column <- "AA"
    } else {
      n_middle_text_columns <- 13
      last_column <- "AB"
    }
    patient_col_types <- c("text","numeric","text","date","date", #ID - Date of admission
                           rep("text", n_middle_text_columns), #Transferred from - Smoking
                           rep("numeric", 5), # BMI - Albumin
                           "text", "text", "date", "text", "text") #Discontinued medication - Transferred to
  }

  patients_range <- paste0("A2:",last_column, "100")


  patient_data_raw <- readxl::read_excel(file, sheet = patients_sheet(lang), range = patients_range, na = c("NA","N/A", "na"), col_types = patient_col_types) %>%
    translate_patient_columns(lang)

  patient_data_raw <- patient_data_raw %>%
    mutate(`Outcome (Discharged/Hospitalized/Transferred/Death)` = translate_outcomes(`Outcome (Discharged/Hospitalized/Transferred/Death)`, lang))

  if(file_version != "Special") {
    patient_data_raw <- patient_data_raw %>% mutate(`Obesity` = NA_character_)
  } else {
    # Using 2019 as anchor because for the "Special" format most data is before June, so presumably the
    # patients are more likely born later in the year
    # Lab markers at admission will be extracted from the lab data
    patient_data_raw <- patient_data_raw %>% mutate(Age = 2019 - `Year of birth`,
                                                    `Creatinin [μmol/L]` = NA_real_,
                                                    `PT INR (Quick)`  = NA_real_,
                                                    `Albumin [g/L]` = NA_real_
                                                    ) %>%
      select(-`Year of birth`)
  }

  if(file_version == "Old") {
    patient_data_raw <- patient_data_raw %>% mutate(`Covid Medications before admission` = NA_character_)
  }

  patient_data_raw
}


read_progression_data <- function(file, hospital_id, lang, file_version, patient_data_typed) {
  max_day <- 172
  progression_data_raw <- readxl::read_excel(file, sheet = progression_sheet(lang), range = "A2:FU2539", na = c("NA",""),
                                             col_types = c(c("text","text", "skip", "text", "text"), rep("text", max_day)),
                                             col_names = c("patient_id", "first_day","indicator", "unit", paste0("Day_", 1:max_day)))

  if(!all(is.na(progression_data_raw[[paste0("Day_", max_day)]]))) {
    stop("Potentially too few days read")
  }

  progression_data <- progression_data_raw %>%
    #translate_progression_columns(lang) %>%
    #rename(patient_id = `Patient ID`, indicator = Indicator, first_day = First_Day, unit = Unit) %>%
    filter(!is.na(indicator) | !is.na(unit)) %>%
    mutate(indicator = translate_markers(indicator, lang),
           hospital_id = hospital_id)

  # Remove the "Don't edit examples" comment
  progression_data$patient_id[2:9] <- NA_character_

  if(hospital_id == "QKuFp") {
    progression_data <- progression_data %>% mutate(indicator = if_else(!is.na(unit) & unit == "IgG", "IgG", indicator))
  }
  if(hospital_id == "ZDYlY") {
    progression_data <- progression_data %>%
      mutate(Day_12 = if_else(
        indicator == "Breathing" & grepl("z oxygenoterapie na AA$", Day_12),
        NA_character_,
        Day_12),
        indicator = gsub(" [0-9,]* ?m?g( a [0-9]*h)?$", "", indicator)) # Remove stuff like "800mg a 4h"
  }

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
    filter(!is.na(indicator), !is.na(value))


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
                        value == "SPONT" ~ "AA",
                        value == "HFNO" ~ "Oxygen",
                        value == "UPV" ~ "MV",
                        TRUE ~ value)
    ) %>%
    select(-indicator, -unit)

  #When there is "Supp. O2" but no Breathing, we can infer "Oxygen"
  suppO2_unknown_breathing <- progression_data %>% filter(indicator == "Supp. O2") %>%
    anti_join(breathing_data, by = c("patient_id","day"))

  if(nrow(suppO2_unknown_breathing) > 0) {
    cat(nrow(suppO2_unknown_breathing), " breathing Oxygen imputed from Supp. O2\n")
    breathing_data <- rbind(breathing_data,
                            suppO2_unknown_breathing %>%
                              mutate(value = "Oxygen") %>%
                              select(-indicator, -unit)
                            )
  }



  duplicate_days <- breathing_data %>% group_by(patient_id, day) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)

  if(nrow(duplicate_days) > 0) {
    print(duplicate_days)
    stop("Duplicate days")
  }


  if(nrow(breathing_data) < 5) {
    stop("Too few breathing data")
  }





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

  breathing_data <- correct_first_day_admission(patient_data_typed, breathing_data)


  # Add final outcomes from patient data
  max_breathing_data <- breathing_data %>% group_by(hospital_id, patient_id) %>%
    summarise(max_day = max(day))

  final_outcomes <- patient_data_typed %>% filter(outcome %in% c("Discharged", "Death")) %>%
    inner_join(max_breathing_data, by = c("hospital_id", "patient_id")) %>%
    transmute(hospital_id = hospital_id,
              patient_id = patient_id,
              value = outcome,
              day = if_else(outcome == "Death", # Death can overwrite last breathing (in the spirit of taking the "worse" case
                            pmax(max_day, last_record),
                            pmax(max_day + 1, last_record)))

  # Join, Potentially overwrite breathing with "Death"
  breathing_data <- breathing_data %>%
    anti_join(final_outcomes %>% filter(value == "Death"), by = c("hospital_id", "patient_id", "day")) %>%
    rbind(final_outcomes)

  # Translate to factor
  to_breathing_factor <- function(value) {
    factor(tolower(value), levels = tolower(disease_levels), labels = disease_levels, ordered = TRUE)
  }
  breathing_data <- breathing_data %>% mutate(
    breathing = to_breathing_factor(
      case_when(
        value == "ARO" ~ "MV",
        value %in% c("Oxygen or better", "Oxygen (or better)") ~ "Oxygen",
        value == "Likely AA" ~ "AA",
        TRUE ~ value
      )),
    breathing_low = to_breathing_factor(
      case_when(
        value == "ARO" ~ "Oxygen",
        value %in% c("Oxygen or better", "Oxygen (or better)", "Likely AA") ~ "AA",
        TRUE ~ value
      )),
    breathing_high = to_breathing_factor(
      case_when(
        value == "ARO" ~ "MV",
        value %in% c("Oxygen or better", "Oxygen (or better)", "Likely AA") ~ "Oxygen",
        TRUE ~ value
      ))
  )

  breathing_unrecognized <- breathing_data %>%
    filter(is.na(breathing) | is.na(breathing_low) | is.na(breathing_high))
  if(nrow(breathing_unrecognized) > 0) {
    print(unique(breathing_unrecognized$value))
    stop("Breathing NAs check failed")
  }

  breathing_data <- breathing_data %>% select(-value)




  unmatched_patients <- patient_data_typed %>% anti_join(breathing_data, by = c("patient_id", "hospital_id"))
  if(nrow(unmatched_patients) > 0) {
    cat("Some patients have no breathing data\n")
    print(unmatched_patients)
  }

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
    mutate(censored = "none", unit = if_else(marker == "pcr_value", "Ct", ""))

  # Handle NT for SpO2_native
  marker_data <- marker_data %>%
    mutate(value = if_else(marker == "SpO2_native" & tolower(value) == "nt", "< 92", value) )

  # Some special cases
  # TODO ask at hospital
  if(hospital_id == "QKuFp") {
    marker_data <- marker_data %>% mutate(value = if_else(marker == "SpO2" & value == "NEG", NA_character_, value))

    # The patient got HCQ before admission, but time unknown (should not matter for most models)
    marker_data <- rbind(marker_data,
                         data.frame(hospital_id = "QKuFp", patient_id = "12",
                                    day = -5, marker = "hcq",
                                    value = ">100", unit = "mg/day"))
  } else if(hospital_id == "YqNbe") {
    marker_data <- marker_data %>% mutate(value = if_else(marker %in% c("hcq","az") & value == "ano", ">100", value))
  } else if(hospital_id == "JXnnM") {
    marker_data <- marker_data %>% mutate(value = gsub("mg", "", value))
  }


  # Handle censoring on everything except PCR
  marker_data <- marker_data %>% filter(marker != "pcr_value") %>%
    mutate(censored = case_when(startsWith(value, "<") ~ "right",
                                startsWith(value, ">") ~ "left",
                                TRUE ~ "none"),
           value_parsed =  as.numeric(if_else(tolower(value) == "no", "0",
                                              gsub("^[<>]? *", "", gsub(",", ".", value, fixed = TRUE)))))


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

  if(file_version == "Special") {
    marker_data_lab <- read_special_lab_data_all(dir = paste0(dirname(file), "/Laborky"),
                                                 hospital_id = hospital_id,
                                                 patient_ids = patient_data_typed$patient_id,
                                                 first_days = patient_data_typed$admission)


    duplicates <- rbind(marker_data, marker_data_lab) %>%
      group_by(marker, patient_id, day) %>%
      summarise(count = n(), n_unique_vals = length(unique(value))) %>%
      filter(count > 1)

    duplicates_problematic <- duplicates %>%
      filter(marker != "pcr_positive" | n_unique_vals > 1)

    if(nrow(duplicates_problematic) > 0) {
      print(duplicates_problematic)
      stop("Duplicates in special lab data")
    }

    # Remove non-problematic duplicates (basically only PCR)
    marker_data_lab <- marker_data_lab %>%
      anti_join(marker_data, by = c("patient_id", "marker", "day"))

    marker_data <- rbind(marker_data, marker_data_lab)
  }

  duplicate_markers <- marker_data %>% group_by(hospital_id,patient_id, marker, day) %>%
    summarise(count = n(), values = paste(value, collapse = ", "), .groups = "drop") %>%
    filter(count > 1)

  if(nrow(duplicate_markers) > 0)  {
    print(duplicate_markers)
    stop("Duplicate markers")
  }


  # Some special handling
  all_overwrites <- unit_overwrites$all
  if(hospital_id %in% names(unit_overwrites)) {
    all_overwrites <- c(all_overwrites, unit_overwrites[[hospital_id]])
  }

  for(overwrite in all_overwrites) {
    marker_data <- marker_data %>%
      mutate(unit = if_else(marker == overwrite$marker & (unit == overwrite$old_unit | (overwrite$old_unit == "" & is.na(unit))), overwrite$new_unit, unit))
  }



  # Default unit translations
  target_rows <- marker_data$unit %in% names(unit_map)
  marker_data$unit[target_rows] <- unit_map[marker_data$unit[target_rows]]

  #Default unit conversions
  for(conversion in unit_conversions) {
      marker_data <- marker_data %>%
        mutate(value = if_else(marker %in% conversion$markers & unit == conversion$old_unit, conversion$mult * value, value),
               unit =  if_else(marker %in% conversion$markers & unit == conversion$old_unit, conversion$new_unit, unit))
  }

  # Check that marker and breathing data are within admission
  markers_allowed_before_admission <- c("pcr_value", "pcr_positive", "hcq", "az")
  markers_allowed_after_discharge <- c("pcr_value", "pcr_positive")
  ranges_markers <- marker_data %>%
    group_by(patient_id) %>%
    summarise(min_day = min(c(0, day[!(marker %in% markers_allowed_before_admission)])),
              max_day = max(c(0, day[!(marker %in% markers_allowed_after_discharge)])),
              markers_min = paste(marker[day == min_day], collapse = ","),
              markers_max = paste(marker[day == max_day], collapse = ","),
              .groups = "drop")

  markers_out_of_range <- patient_data_typed %>%
    select(patient_id, last_record) %>%
    inner_join(ranges_markers, by = "patient_id") %>%
    filter(min_day < -1 | max_day > last_record)

  if(hospital_id == "bAbUl") {
    # Checked manually
    markers_out_of_range <- markers_out_of_range %>% filter(min_day < - 2)
  }

  if(nrow(markers_out_of_range) > 0) {
    print(markers_out_of_range)
    stop("Markers out of range")
  }

  ranges_breathing <- breathing_data %>%
    group_by(patient_id) %>%
    summarise(min_day = min(day), max_day = max(day),
              breathing_min = breathing[day == min_day],
              breathing_max = breathing[day == max_day],
              .groups = "drop")

  breathing_out_of_range <- patient_data_typed %>%
    select(patient_id, last_record) %>%
    inner_join(ranges_breathing, by = "patient_id") %>%
    filter(min_day < -1 | max_day > last_record + 1)

  if(nrow(breathing_out_of_range) > 0) {
    print(breathing_out_of_range)
    stop("Breathing out of range")
  }

  #TODO check SuppO2 only with Oxygen

  loo::nlist(marker_data, breathing_data, adverse_events_data)
}


join_hospital_data <- function(hospital_data) {
  res <- list()
  components <- c("patient_data", "breathing_data", "marker_data", "adverse_events_data")
  for(component in components) {
    component_list <- list()
    for(i in 1:length(hospital_data)) {
      component_list[[i]] <- hospital_data[[i]][[component]]
    }
    res[[component]] <- do.call(rbind, component_list)
  }
  res
}

check_patient_overlap <- function(hospital_data_1, hospital_data_2, patients_to_merge, hospitals_to_merge) {
  pts1 <- hospital_data_1$patient_data
  pts2 <- hospital_data_2$patient_data

  hospital_id1 <- unique(pts1$hospital_id)
  hospital_id2 <- unique(pts2$hospital_id)
  if(length(hospital_id1) > 1 || length(hospital_id2) > 1) {
    stop("Non-unique hospital_id in hospital")
  }

  should_merge <- nrow(
    hospitals_to_merge %>% filter(hospital_id1 == !!hospital_id1 & hospital_id2 == !!hospital_id2)
  ) > 0

  # if(should_merge) {
  #  cat("Should merge", hospital_id1, hospital_id2)
  # }

  merged_ids <- patients_to_merge %>%
    mutate(merged = paste0(hospital_id1, "_", patient_id1, "::", hospital_id2, "_", patient_id2)) %>%
    pull(merged)

  pts1 %>% inner_join(pts2, by = "sex", suffix = c("1","2")) %>%
    filter(
      should_merge | !is.na(transferred_from1) | !is.na(transferred_from2) |
        !is.na(transferred_to1) | !is.na(transferred_to2),

      abs(age1 - age2) <= 1,
      is.na(BMI1) | is.na(BMI2) | abs(BMI1 - BMI2) < 5,
      is.na(outcome1) | is.na(outcome2) |  !(outcome1 == "Death" & outcome2 == "Discharged"),
      is.na(outcome1) | is.na(outcome2) |  !(outcome2 == "Death" & outcome1 == "Discharged"),
      outcome1 != "Discharged" | outcome2 != "Discharged",
      is.na(symptom_onset1) | is.na(symptom_onset2) |
        abs(symptom_onset2 - symptom_onset1) <  lubridate::make_difftime(day = 5),
      # # Filter those already merged
      !(paste0(hospital_id1, "_", patient_id1, "::", hospital_id2, "_", patient_id2) %in% merged_ids),
      !(paste0(hospital_id2, "_", patient_id2, "::", hospital_id1, "_", patient_id1) %in% merged_ids)

    )

  # Matches by breathing - unreliable
  # breathe1 <- hospital_data_1$breathing_data %>%
  #   filter(!is.na(breathing)) %>%
  #   inner_join(pts1 %>% select(patient_id, admission), by = "patient_id") %>%
  #   mutate(date = admission + lubridate::make_difftime(day = day))
  #
  # breathe2 <- hospital_data_2$breathing_data %>%
  #   filter(!is.na(breathing)) %>%
  #   inner_join(pts2 %>% select(patient_id, admission), by = "patient_id") %>%
  #   mutate(date = admission + lubridate::make_difftime(day = day))
  #
  # breathe1 %>%
  #   inner_join(breathe2, by = c("date"), suffix = c("1","2")) %>%
  #   group_by(hospital_id1,hospital_id2, patient_id1, patient_id2) %>%
  #   summarise(n_matches = sum(breathing1 == breathing2), n_overlap = n(), .groups = "drop") %>%
  #   #filter(n_matches > 5)
  #   filter(n_overlap > 5, n_matches > pmax(0.8 * n_overlap, n_overlap - 2))

}

merge_patients <- function(complete_data, patients_to_merge) {
  actually_merge <- patients_to_merge %>% filter(resolution != "ignore")
  for(i in 1:nrow(actually_merge)) {
    row <- actually_merge[i, ]
    if(row$resolution == "use2") {
      complete_data <- complete_data %>%
        filter_complete_data(!(hospital_id  == row$hospital_id1 & patient_id == row$patient_id1))
    } else if(row$resolution == "merge_to_1") {
      patient_row_2 <- complete_data$patient_data %>%
        filter((hospital_id  == row$hospital_id2 & patient_id == row$patient_id2))
      if(nrow(patient_row_2) != 1) {
        stop("Problems in merge")
      }

      # Remove patient row 2
      complete_data$patient_data <- complete_data$patient_data %>%
        filter(!(hospital_id  == row$hospital_id2 & patient_id == row$patient_id2))
      # Merge marker rows
      complete_data$marker_data <- complete_data$marker_data %>%
        mutate(
          should_change = hospital_id  == row$hospital_id2 & patient_id == row$patient_id2,
          hospital_id = if_else(should_change, row$hospital_id1, hospital_id),
          patient_id = if_else(should_change, row$patient_id1, patient_id),
          day = if_else(should_change, day + row$shift, day)
        ) %>%
        #Avoid duplicates from the transfer day. Currently sensible for all transfers to just
        #Drop the data from hospital 2
        filter(!(should_change & day == row$shift)) %>%
        select(-should_change)

      # Merge breathing rows
      complete_data$breathing_data <- complete_data$breathing_data %>%
        mutate(
          should_change = hospital_id  == row$hospital_id2 & patient_id == row$patient_id2,
          hospital_id = if_else(should_change, row$hospital_id1, hospital_id),
          patient_id = if_else(should_change, row$patient_id1, patient_id),
          day = if_else(should_change, day + row$shift, day)
        ) %>%
        #Avoid duplicates from the transfer day. Currently sensible for all transfers to just
        #Drop the data from hospital 2
        filter(!(should_change & day == row$shift)) %>%
        select(-should_change)

      # Update outcome and last record
      complete_data$patient_data <- complete_data$patient_data %>%
        mutate(should_change = hospital_id == row$hospital_id1 & patient_id == row$patient_id1,
               last_record = if_else(should_change, last_record + patient_row_2$last_record, last_record),
               outcome = if_else(should_change, patient_row_2$outcome, outcome)) %>%
        select(-should_change)
    } else {
      stop(paste0("Unknown resolution: ", row$resolution))
    }
  }

  duplicate_breathing <- complete_data$breathing_data %>%
    group_by(hospital_id, patient_id, day) %>%
    summarise(count = n()) %>% filter(count > 1)

  duplicate_markers <- complete_data$marker_data %>%
    group_by(hospital_id, patient_id, marker, day) %>%
    summarise(count = n()) %>% filter(count > 1)

  if(nrow(duplicate_breathing) > 0 || nrow(duplicate_markers) > 0) {
    print(duplicate_breathing)
    print(duplicate_markers)
    stop("Duplicates")
  }

  complete_data
}

merge_hospitals <- function(complete_data, hospitals_to_merge) {
  for(i in 1:nrow(hospitals_to_merge)) {
    row <- hospitals_to_merge[i, ]
    if(row$unique_ids != "no") {
      ids1 <- c(
        complete_data$patient_data %>% filter(hospital_id == row$hospital_id1) %>% pull(patient_id),
        complete_data$marker_data %>% filter(hospital_id == row$hospital_id1) %>% pull(patient_id),
        complete_data$breathing_data %>% filter(hospital_id == row$hospital_id1) %>% pull(patient_id)
      )
      ids2 <- c(
        complete_data$patient_data %>% filter(hospital_id == row$hospital_id2) %>% pull(patient_id),
        complete_data$marker_data %>% filter(hospital_id == row$hospital_id2) %>% pull(patient_id),
        complete_data$breathing_data %>% filter(hospital_id == row$hospital_id2) %>% pull(patient_id)
      )
      if(length(intersect(ids1, ids2) > 0)) {
        print(row)
        print(intersect(ids1, ids2))
        stop("IDs not unique")
      }
    }
    complete_data <- complete_data %>%
      mutate_complete_data(patient_id = if_else(hospital_id == row$hospital_id2,
                                                paste0(row$hospital_id2,"_", patient_id),
                                                patient_id),
                           hospital_id = if_else(hospital_id == row$hospital_id2,
                                                 row$hospital_id1,
                                                 hospital_id)
      )
  }
  complete_data
}

filter_complete_data <- function(complete_data, ...) {
  for(n in names(complete_data)) {
    complete_data[[n]] <- complete_data[[n]] %>% filter(...)
  }

  complete_data
}

mutate_complete_data <- function(complete_data, ...) {
  for(n in names(complete_data)) {
    complete_data[[n]] <- complete_data[[n]] %>% mutate(...)
  }

  complete_data
}


anonymize_for_analysis <- function(complete_data) {
  complete_data <- complete_data %>% mutate_complete_data(patient_id_full = paste0(hospital_id,"_", patient_id))

  n_pts <- nrow(complete_data$patient_data)
  repeat {
    new_ids <- stringi::stri_rand_strings(n_pts, 10)
    if(length(unique(new_ids)) == n_pts) {
      break
    }
  }


  names(new_ids) <- complete_data$patient_data$patient_id_full

  write_csv(data.frame(old = names(new_ids), new = new_ids), path = "H:/raw/patient_map.csv")

  anonymized_data <- complete_data %>% mutate_complete_data(patient_id = new_ids[patient_id_full])


  # TODO - once coded, adverse events can be included in anonymized
  anonymized_data$adverse_events_data <- NULL

  for(n in names(anonymized_data)) {
    anonymized_data[[n]]$patient_id_full <- NULL
  }

  cols_to_remove <- c("symptom_onset", "admission", "transferred_from", "transferred_to",
                      "last_record_date", "best_supportive_care_from_date", "discontinued_medication_reason")
  for(c in cols_to_remove) {
    anonymized_data$patient_data[[c]] <- NULL
  }

  # Remove old ordering
  anonymized_data$patient_data <- anonymized_data$patient_data %>% arrange(hospital_id, patient_id)
  anonymized_data$breathing_data <- anonymized_data$breathing_data %>% arrange(hospital_id, patient_id, day)
  anonymized_data$marker_data <- anonymized_data$marker_data %>% arrange(hospital_id, patient_id, day)

  anonymized_data
}

