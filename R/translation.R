translations_file <-  here::here("public_data", "Translations.xlsx")
sheet_name_translations <- readxl::read_excel(translations_file, "SheetNames")
patient_columns_translations <- readxl::read_excel(translations_file, "PatientColumns")
outcome_translations <- readxl::read_excel(translations_file, "Outcomes")
progression_columns_translations <- readxl::read_excel(translations_file, "ProgressionColumns")
marker_translations <- readxl::read_excel(translations_file, "Markers")


yes_variants <- c("yes","y", "ano", "a")
no_variants <- c("no","n", "ne")


version_sheet <- function(lang) {
  if_else(lang == "En", "I need help",
          sheet_name_translations %>% filter(Sheet == "version") %>% pull(Cz)
  )
}

patients_sheet <- function(lang) {
  if_else(lang == "En", "Patients", "Pacienti")
}

progression_sheet <- function(lang) {
  if_else(lang == "En", "Disease Progression",
          sheet_name_translations %>% filter(Sheet == "progression") %>% pull(Cz)
  )
}

translate_patient_columns <- function(patient_data_raw, lang) {
  if(lang == "En") {
    patient_data_raw
  } else {
    translation <- patient_columns_translations$Column
    names(translation) <- patient_columns_translations$Cz
    names(patient_data_raw) <- translation[names(patient_data_raw)]

    patient_data_raw
  }
}

translate_outcomes <- function(outcomes, lang) {
  if(lang == "En") {
    outcomes
  } else {
    translation <- outcome_translations$Outcome
    names(translation) <- outcome_translations$Cz

    translation[outcomes]
  }

}

translate_progression_columns <- function(progression_data_raw, lang) {
  if(lang == "En") {
    progression_data_raw
  } else {
    translation <- progression_columns_translations$Column
    names(translation) <- progression_columns_translations$Cz
    names_to_change <- names(progression_data_raw) %in% names(translation)
    names(progression_data_raw)[names_to_change] <- translation[names(progression_data_raw)[names_to_change]]

    progression_data_raw
  }
}

translate_markers <- function(markers, lang) {
  if(lang == "En") {
    markers
  } else {
    translation <- marker_translations$Marker
    names(translation) <- marker_translations$Cz
    rows_to_change <- markers %in% names(translation)
    markers[rows_to_change] <- translation[markers[rows_to_change]]

    markers
  }
}
