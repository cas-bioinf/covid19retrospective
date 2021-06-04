special_lab_data_ignored_patterns <- c(
  "^.adatel:",
  "^Pl.tce ?:",
  "^Protokol",
  "^-*$",
  "^Odd.leni klinick. mikrobiologie$",
  "^DIF - STANOVEN. ANALYZ.TOREM$",
  "^Nen. zn.mo, zda pacient ",
  "^Hodn\\.:",
  "^[0-9]{5} -"
)

special_lab_data_sections <- c(
  "Plazma",
  "KREVN\U00CD OBRAZ",
  "DIF ANALYZ\U00C1TOR",
  "KOAGULACE",
  "Biochemie",
  "Mo\U010D chemicky a sedim",
  "Mo\U00E8 chemicky a sedim",
  "Mo\U010D",
  "Astrup",
  "DIF MIKROSKOPICKY",
  "Imunologie",
  "RETIKULOCYTY",
  "Transfuzn\U00ED",
  "Punkt\U00E1t",
  "L\U00E9\U010Diva",
  "Markery hepatitid",
  "Toxikologie"

)

special_lab_data_markers_to_collect <-
  c("INR" = "pt_inr",
    "Ddim" = "d_dimer",
    "Ly#" = "lymphocyte_count",
    "LYa" = "lymphocyte_count",
    "Ne#" = "neutrophile_count",
    "NEa" = "neutrophile_count",
    "Eo#" = "eosinophile_count",
    "EOa" = "eosinophile_count",
    "PLT" = "platelet_count",
    "Kre" = "creatinin",
    "CRP" = "CRP",
    "PCT" = "procalcitonin",
    "Alb" = "albumin",
    "Ferr" = "ferritin",
    "IL-6" = "IL_6",
    "IL6" = "IL_6",
    "AST" = "AST",
    "LD" = "lactate_dehydrogenase",
    "Vyska" = "height",
    "Vaha" = "weight"
  )

special_lab_data_markers_to_collect_units <-
  c("INR" = "INR",
    "Ddim" = "ng/ml DDU",
    "Ly#" = "10^9/l",
    "LYa" = "10^9/l",
    "Ne#" = "10^9/l",
    "NEa" = "10^9/l",
    "Eo#" = "10^9/l",
    "EOa" = "10^9/l",
    "PLT" = "10^9/l",
    "Kre" = "\U03BCmol/l",
    "CRP" = "mg/l",
    "PCT" = "\U03BCg/l",
    "Alb" = "g/l",
    "Ferr" = "\U03BCg/l",
    "IL-6" = "ng/l",
    "IL6" = "ng/l",
    "AST" = "\U03BCkat/l",
    "LD" = "\U03BCkat/l",
    "Vyska" = "cm",
    "Vaha" = "kg"
  )

if(!identical(names(special_lab_data_markers_to_collect), names(special_lab_data_markers_to_collect_units))){
  stop("Bad names")
}

read_special_lab_data <- function(input_file, patient_id, hospital_id, first_day) {
  lines <- readLines(input_file)
  #close(in_con)
  Encoding(lines) <- "UTF-8"

  res <- list()
  current_day <- NULL
  in_section <- NULL
  covid_test_state <- "no"
  in_comments <- FALSE

  line_numbers <- 1:length(lines)

  for(ig in special_lab_data_ignored_patterns) {
    to_include <- !grepl(ig, trimws(lines))
    lines <- lines[to_include]
    line_numbers <- line_numbers[to_include]
  }


  for(i in 1:length(lines)) {
    s <- trimws(lines[i], which = "right")


    if(grepl("^V.sledky z [0-9][0-9]/[0-9][0-9]/[0-9][0-9]:$", s)) {
      date_time_str <- substr(s, 12, 19)
      current_date <- lubridate::dmy(date_time_str)
      current_day <- as.integer(current_date - first_day)
      in_section <- NULL
      in_comments <- FALSE
      if(covid_test_state == "started") {
        print(paste0("Line ", line_numbers[i]))
        stop("Covid test not stopped before next day")
      }
      covid_test_state <- "no"
    } else {
      matched <- FALSE

      for(sec in special_lab_data_sections) {
        if(s == paste0("  ",sec,":")) {
          in_section <- sec
          in_comments <- FALSE
          if(covid_test_state == "started") {
            print(paste0("Line ", line_numbers[i]))
            stop("Covid test not stopped before section")
          }
          covid_test_state <- "no"
          matched <- TRUE
        }
      }

      if(!matched) {
        if(grepl("^  N.lezy:$", s)) {
          covid_test_state <- "in_nalezy"
          matched <- TRUE
          in_section <- NULL
        } else if(grepl("^  (Koment..e a v.po.ty|Nedefinovan. t..da):$", s)) {
          in_comments <- TRUE
          in_section <- NULL
          if(covid_test_state == "started") {
            print(paste0("Line ", line_numbers[i]))
            stop("Covid test not stopped before comments")
          }
          covid_test_state <- "no"
          matched <- TRUE
        }
      }

      if(matched) {
        next
      }

      if(!is.null(in_section)) {
        if(s == "") {
          in_section <- NULL
        } else if(grepl("^      [A-Za-z0-9][^:]*: [^ ]", s)) {
          marker_value_part <- strsplit(trimws(s), ": ", fixed = TRUE)[[1]]
          if(length(marker_value_part) != 2) {
            print(paste0("Line ", line_numbers[i]))
            print(s)
            stop("Unexpected marker format")
          }
          marker_raw <- marker_value_part[1]

          if(marker_raw %in% names(special_lab_data_markers_to_collect)) {
            vals <- strsplit(marker_value_part[2], "; ", fixed = TRUE)[[1]]

            #Ignored values
            vals <- vals[!grepl("^(Txt\\+Hisviz koment..|krev sra.en.|hemol.za|chyl.zn.|nelze)$", vals)]
            if(length(vals) == 0) {
              next
            }
            value_pattern <- "^(<?>?-?[0-9]+,?[0-9]*|(Txt\\+His)?v.ce ne. [0-9]+ .*)$"
            if(!all(grepl(value_pattern, vals))) {
              print(paste0("Line ", line_numbers[i]))
              print(s)
              stop("Bad marker value")
            }

            values <- as.numeric(
              gsub("^(<|>)","",
                gsub("(Txt\\+His)?v.ce ne. ([0-9]+) [^;]*", "\\2",
                                  gsub(",", ".", vals, fixed = TRUE))))

            if(any(is.na(values))) {
              print(paste0("Line ", line_numbers[i]))

              print(s)
              print(vals)
              stop("NAs in conversion")
            }
            all_censored <- case_when(
              grepl("v.ce ne.", vals) ~  "right",
              grepl("^>", vals) ~  "right",
              grepl("^<", vals) ~ "left",
              TRUE ~ "none"
            )

            #this is potentially problematic, but likely OK for all the stuff we care about
            index_to_save <- which.max(values)
            value <- values[index_to_save]
            censored <- all_censored[index_to_save]

            res[[length(res) + 1]] <- tibble(patient_id = patient_id,
                                             hospital_id = hospital_id,
                                             marker = special_lab_data_markers_to_collect[marker_raw],
                                             unit = special_lab_data_markers_to_collect_units[marker_raw],
                                             day = current_day,
                                             value = value,
                                             censored = censored,
                                             section = in_section)
          }
        } else {
          print(paste0("Line ", line_numbers[i]))
          print(paste0("Weird line: ",s))
          #stop("Invalid line in section")
        }
        next
      }

      if(is.null(current_day)) {
        stop("Results before date")
      }

      if(grepl("Typ vzorku: St.r Covid 19", s)) {
        if(covid_test_state != "in_nalezy") {
          stop("Invalid covid transition")
        }
        covid_test_state <- "started"
      } else if(grepl("^(Koment..:.*neprovedeno|N.lez: Nejasn. v.sledek)", s)) {
        covid_test_state <- "in_nalezy"
      } else if(grepl("^N.lez: (P.vodce )?Covid 19", s)) {
        if(covid_test_state != "started") {
          stop("Invalid covid transition")
        }

        if(grepl("pozitivn.$", s)) {
          value = 1
        } else if(grepl("negativn.$", s)) {
          value = 0
        } else {
          stop("Invalid Covid test res")
        }

        res[[length(res) + 1]] <- tibble(patient_id = patient_id,
                                         hospital_id = hospital_id,
                                         marker = "pcr_positive",
                                         unit = "",
                                         day = current_day,
                                         value = value,
                                         censored = "none",
                                         section = "")
        covid_test_state <- "in_nalezy"
      } else if(grepl("^Rapid test Covid", s)) {
        if(covid_test_state != "in_nalezy") {
          stop("Invalid covid transition")
        }
      } else if(grepl("^Negativn.*nevylu.uje.*COVID", s, ignore.case = TRUE)) {
        if(covid_test_state != "in_nalezy") {
          stop("Invalid covid transition")
        }
      } else if(grepl("^Koment..", s)) {
        next
      } else {
        if(in_comments) {
          next
        }
        if(grepl("covid", s, ignore.case = TRUE)) {
          print(paste0("Line ", line_numbers[i]))
          print(paste0("Unrecognized Covid row:", s))
        }

        if(covid_test_state != "no") {
          next
        }

        print(paste0("Line ", line_numbers[i]))
        print(s)
        stop("Unrecognized row")
      }
    }
  }
  res_all <- do.call(rbind, res)

  if(any(is.na(res_all$value))) {
    print(res_all %>% filter(is.na(value)))
    stop("NAs in values")
  }

  res_all
}


read_special_lab_data_all <- function(dir, hospital_id, patient_ids, first_days) {
  res <- list()
  for(i in 1:length(patient_ids)) {
    filename <- paste0(dir, "/P", patient_ids[i], ".txt")
    if(!file.exists(filename)) {
      warning(paste0(filename, " does not exist\n"))
    } else {
      print(filename)
      res[[i]] <- read_special_lab_data(filename, patient_ids[i], hospital_id, first_days[i])
    }
  }

  all <- do.call(rbind, res) %>%
    filter(!startsWith(section,"Punkt")) # Not include punction results

  # Inspect sections where markers occur
  all %>% group_by(marker, section) %>% summarise(n()) %>% print()

  # Plots to check values across sections
  # for(m in c("albumin", "creatinine", "CRP", "IL_6", "lymphocyte_count","procalcitonin")) {
  #   print(all %>% filter(marker == m) %>%
  #           ggplot(aes(x = as.numeric(value))) + geom_histogram() + facet_wrap(~section, ncol = 1)+ ggtitle(m))
  # }


  # Merging serum and plasma measurements (also microscopy vs. analyzator Lymphocyte counts)
  # This is not 100% clean, but the values seem to not be _too_ different
  # Once again, since
  merged_sections <- all %>%
    rename(value_raw = value, censored_raw = censored) %>%
    group_by(hospital_id, patient_id, marker, day, unit) %>%
    summarise(censored = censored_raw[which.max(value_raw)],
              value = max(value_raw))

  merged_sections
}

remove_names_from_special_lab_data <- function(input_file, output_file) {
  #in_con <- file(input_file, encoding = "UTF-8", open = "rt")
  lines <- readLines(input_file)
  Encoding(lines) <- "UTF-8"
  #close(in_con)
  lines <- lines[!grepl("Pacient:|[Ll].ka. *:|Schv.lil\\(a\\)|Datum hl..en.|Potvrzuj", lines)]

  out_con <- file(output_file, open = "w+", encoding = "native.enc")
  writeLines(lines, con = out_con, useBytes = TRUE)
  close(out_con)
}

remove_names_from_special_lab_data_dir <- function(input_dir, output_dir) {
  for(f in list.files(input_dir)) {
    remove_names_from_special_lab_data(paste0(input_dir, "/", f), paste0(output_dir, "/", f))
  }
}
