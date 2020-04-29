#' Generates fake dataset starting from breathing support histories taken from
#' Gren et al.: https://www.nejm.org/doi/full/10.1056/NEJMoa2007016
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr "%>%"
fake_data_grein <- function() {
  breathing_data_orig <- readxl::read_excel(here::here("public_data", "Grein_et_al.xlsx"), na = c("")) %>%
    select(-Remdesivir_Last_Day) %>%
    pivot_longer(starts_with("Day_"), names_prefix = "Day_", names_to = "day", values_to = "breathing") %>%
    filter(!is.na(breathing)) %>%
    transmute(
      patient = Patient,
      breathing = if_else(breathing %in% c("LoO2", "HiO2"), "Oxygen", breathing),
      day = as.integer(day)) %>%
    group_by(patient) %>%
    mutate(breathing = if_else(breathing == "ND", breathing[which.max(day[breathing %in% c("AA", "Oxygen", "NIPPV", "MV", "ECMO")])], breathing)) %>%
    ungroup() %>%
    mutate(
      breathing = factor(breathing, levels = c("Discharge", "AA", "Oxygen", "NIPPV", "MV", "ECMO", "Death"), ordered = TRUE)
    )

  #Generate "Lead-in" data
  breathing_data_list <- list()
  for(p in unique(breathing_data_orig$patient)) {
    breathing_p <- breathing_data_orig %>% filter(patient == p)
    start_value <- breathing_p %>% filter(day == 0) %>% pull(breathing)
    if(start_value >= "NIPPV") {
      n_days_before <- 1 + rnbinom(1, mu = (as.integer(start_value) - 2) * 2.5, size = 1)
      n_steps <- as.integer(start_value) - 3
      step_times <- sort(purrr::rdunif(n_steps, a = 1, b = n_days_before))
      breathing_new <- character(n_days_before)
      start <- 1
      for(s in 1:n_steps) {
        if(start <= step_times[s]) {
          breathing_new[start:(step_times[s])] <- levels(breathing_data_orig$breathing)[s + 1]
          start <- step_times[s] + 1
        }
      }
      breathing_new <- factor(breathing_new, levels = levels(breathing_data_orig$breathing), ordered = TRUE)
      if(start <= n_days_before) {
        breathing_new[start:n_days_before] <- start_value
      }

      if(any(is.na(breathing_new))) {
        stop("Error")
      }

      breathing_data_list[[p]] <- rbind(
        data.frame(patient = p, breathing = breathing_new, day = 0: (n_days_before - 1)),
        breathing_p %>% mutate(day = day + n_days_before)
      )
    } else {
      breathing_data_list[[p]] <- breathing_p
    }
  }

  breathing_data <- do.call(rbind, breathing_data_list)


  if(any(is.na(breathing_data$breathing)) || any(is.na(breathing_data$day))){
    stop("Error")
  }

  breathing_summary <- breathing_data %>% group_by(patient) %>%
    summarise(worst_condition = max(breathing),
              first_day_invasive = as.integer(min(c(1e4, day[breathing >= "MV" & breathing < "Death"]))),
              last_day_invasive = as.integer(max(c(-1, day[breathing >= "MV" & breathing < "Death"]))),
              days_hospitalized = max(day) + 1
    ) %>%
    mutate(first_day_invasive = if_else(first_day_invasive > 1e3, NA_integer_, first_day_invasive),
           last_day_invasive = if_else(last_day_invasive < 0, NA_integer_, last_day_invasive))

  patient_data <- breathing_data %>% group_by(patient) %>%
    top_n(1, day) %>%
    inner_join(breathing_summary, by = c("patient" = "patient")) %>%
    group_by(patient) %>%
    mutate(outcome = case_when(breathing == "Discharge" ~ "Discharged",
                               breathing == "Death" ~ "Death",
                               runif(n()) < 0.01 ~ "Transferred",
                               TRUE ~ "Hospitalized"),
           site = sample(c("A","B","C"), size = n(), replace = TRUE, prob = c(0.7,0.2,0.1)),
           sex = if_else(runif(n()) < 0.5, "M", "F" ),
           age = round(rnorm(n(), c(-1, 30, 38, 48, 55, 58, 65)[worst_condition], sd = c(-1, 5, 6, 8, 10, 11, 11)[worst_condition])),
           days_from_symptom_onset = rnbinom(n(), mu = 4, size = 1),
           myocardial_infarction = rpois(n(), 0.05 + (site != "A") * 0.08),
           hypertension_drugs = rnbinom(n(), mu = 0.1 + (site == "B") * 0.6 , size = 0.3),
           smoking = runif(n()) < 0.4 + (site == "C") * 0.15,
           bmi = rnorm(n(), 28, 4) + case_when(site == "A" ~ 0, site == "B" ~ 2.7, site == "C" ~ -0.5) ,
           took_hcq = runif(n()) < 0.6,
           took_az = took_hcq & runif(n()) < 0.8,
           took_kaletra = !took_hcq & site == "B" & runif(n()) < 0.2,
           first_day_antivirals = if_else(took_hcq | took_kaletra, as.integer(purrr::rdunif(n(), a = 0, b = 3)), NA_integer_),
           took_tocilizumab = worst_condition >= "MV" & site == "A" & runif(n()) < 0.1,
           first_day_tocilizumab = first_day_invasive + as.integer(purrr::rdunif(n(), a = 0, b = 1))
    ) %>%
    ungroup() %>%
    select(-breathing, -day)

  breathing_data <- breathing_data %>% filter(!(breathing %in% c("Discharge", "Death"))) %>%
    mutate(breathing = droplevels(breathing))


  marker_data_list <- list()
  next_id <- 1
  for(p in patient_data$patient) {
    breathing_p <- breathing_data %>% filter(patient == p)
    patient_p <- patient_data %>% filter(patient == p)

    #Fit a spline to breathing status
    sp <- smooth.spline(breathing_p$day, as.integer(breathing_p$breathing), nknots = 6)

    oxygen_days <- breathing_p %>% filter(breathing == "Oxygen") %>% pull(day)
    if(length(oxygen_days) > 0) {
      oxygen_flow <- round( (runif(1, 0, 3) + runif(length(oxygen_days)) + predict(sp, oxygen_days)$y - 3) * 2 ) / 2 + 2
      marker_data_list[[next_id]] <- data.frame(patient = p, day = oxygen_days, marker = "oxygen_flow", value = oxygen_flow, censored = "none")
      next_id <- next_id + 1
    }

    if(patient_p$took_hcq) {
      hcq_days <- patient_p$first_day_antivirals:min(patient_p$first_day_antivirals + 5, max(breathing_p$day))
      if(patient_p$site == "B") {
        hcq_dose <- rep(400, length(hcq_days))
        hcq_dose[1] <- 800
      } else {
        hcq_dose <- rep(200, length(hcq_days))
        hcq_dose[1] <- 400
      }
      marker_data_list[[next_id]] <- data.frame(patient = p, day = hcq_days, marker = "hcq", value = hcq_dose,  censored = "none")
      next_id <- next_id + 1


      if(patient_p$took_az) {
        az_dose <- rep(250, length(hcq_days))
        az_dose[1] <- 500
        marker_data_list[[next_id]] <- data.frame(patient = p, day = hcq_days, marker = "az", value = az_dose,  censored = "none")
        next_id <- next_id + 1
      }
    }
    if(patient_p$took_kaletra) {
      kaletra_max_days <- purrr::rdunif(1, a = 4, b = 8)
      kaletra_days <- patient_p$first_day_antivirals:min(patient_p$first_day_antivirals + kaletra_max_days, max(breathing_p$day))
      if(patient_p$took_kaletra) {
        marker_data_list[[next_id]] <- data.frame(patient = p, day = kaletra_days, marker = "kaletra", value = 400,  censored = "none")
        next_id <- next_id + 1
      }
    }
    if(patient_p$took_tocilizumab) {
      marker_data_list[[next_id]] <- data.frame(patient = p, day = patient_p$first_day_tocilizumab, marker = "tocilizumab", value = if_else(patient_p$site == "C", 6, 7),  censored = "none")
      next_id <- next_id + 1
    }

    pcr_days <- c(rbinom(1, size = 2, prob = 0.3), which(runif(max(breathing_p$day) - 5) < 0.1) + 2)
    pcr_positive <- c(TRUE, runif(length(pcr_days) - 1)) < 0.93
    if(patient_p$site == "A" && sum(pcr_positive) > 0) {
      pcr_value <- 35 - (((predict(sp, pcr_days[pcr_positive] - 2)$y / 5) * 21) + 12 + rnorm(sum(pcr_positive), sd = 9))
      pcr_value[pcr_value < 8] <- runif(sum(pcr_value < 8), 8, 10)
      pcr_value[pcr_value > 35] <- runif(sum(pcr_value > 35), 33, 35)

      marker_data_list[[next_id]] <- data.frame(patient = p, day = pcr_days[pcr_positive], marker = "pcr_value", value = pcr_value,  censored = "none")
      next_id <- next_id + 1

    }

    if(patient_p$outcome == "Discharged") {
      first_negative <- purrr::rdunif(1, max(pcr_days) + 1, max(breathing_p$day) - 2)
      pcr_days <- c(pcr_days, first_negative, first_negative + 1 + (runif(1) < 0.4))
      pcr_positive <- c(pcr_positive, FALSE, FALSE)
    }

    marker_data_list[[next_id]] <- data.frame(patient = p, day = pcr_days, marker = "pcr_positive", value = pcr_positive,  censored = "none")
    next_id <- next_id + 1

    lab_days <- which(runif(max(breathing_p$day) + 1) < 0.8) - 1
    spline_vals_norm <- predict(sp, lab_days + 0.3)$y / 5
    spline_vals_norm_crp <- predict(sp, lab_days + 1.5)$y / 5

    crp <- round(exp(spline_vals_norm_crp * log(90) + rnorm(1, sd = 0.3) + rnorm(length(lab_days), sd = log(2))))
    censored_crp <- crp < 1
    crp[censored_crp] <- 1

    d_dimer <- round( exp(spline_vals_norm - 0.8 + rnorm(1, sd = 0.1) + rnorm(length(lab_days), sd = log(1.1))), 2)
    d_dimer[d_dimer < 0.6] <- round(runif(sum(d_dimer < 0.6), 0.5, 0.8))

    ferritin <- round(450 + (1800 - 450 + rnorm(1, sd = 600)) * spline_vals_norm + rnorm(length(lab_days), sd = 100))
    ferritin[ferritin < 150] <- round(runif(sum(ferritin < 150), 150, 300))
    censored_ferritin <- ferritin > 2000
    ferritin[censored_ferritin] <- 2000


    marker_data_list[[next_id]] <- data.frame(patient = p, day = lab_days, marker = "crp", value = crp,  censored = if_else(censored_crp, "left", "none"))
    next_id <- next_id + 1
    marker_data_list[[next_id]] <- data.frame(patient = p, day = lab_days, marker = "d_dimer", value = d_dimer,  censored = "none")
    next_id <- next_id + 1

    if(patient_p$site == "A") {
      marker_data_list[[next_id]] <- data.frame(patient = p, day = lab_days, marker = "ferritin", value = ferritin,  censored = if_else(censored_ferritin, "right", "none"))
    }
    next_id <- next_id + 1

  }

  marker_data <- do.call(rbind, marker_data_list)

  base_spec <- build_wider_spec(marker_data, names_from = "marker", values_from = "value")
  censored_spec <- base_spec %>% mutate(.name = paste0(.name, "_censored"), .value = "censored")

  censored_cols <- c("crp", "ferritin")
  non_censored_cols <- setdiff(unique(marker_data$marker), censored_cols)

  marker_data_wide <- marker_data %>% pivot_wider_spec(spec = rbind(base_spec, censored_spec)) %>%
    select(-one_of(paste0(non_censored_cols, "_censored")))

  list(
    patient_data = patient_data,
    breathing_data = breathing_data,
    marker_data = marker_data,
    marker_data_wide = marker_data_wide
  )
}
