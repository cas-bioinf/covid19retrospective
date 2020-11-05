dummy_markers <- c("Drugs", "Other marker (optional)", "Other Covid drug", "Yet another drug",
                   "RTG S+P", "Rosuvastatin","Rtg S+P")


marker_map <- c(
      ## Markers
        "PCR" = "pcr_value",
        "Supp. O2" = "oxygen_flow",
        "SpO2" = "SpO2",
        "SpO2 native" = "SpO2_native",
        "Horowitz index" = "Horowitz_index",
        "PEEP" = "PEEP",
        "FiO2" = "FiO2",
        "Ferritin" = "ferritin",
        "D-dimer" = "d_dimer",
        "CRP" = "CRP",
        "IL-6" = "IL_6",
        "PCT" = "procalcitonin",
        "Prokalcitonin" = "procalcitonin",
        "Lymphocyte count" = "lymphocyte_count",
        "Leukocyte count" = "leukocyte_count",
        "Leukocytes" = "leukocyte_count",
        "Leukocyte" = "leukocyte_count",
        "creatinine" = "creatinine",
        "urea" = "urea",
        "IgG" = "IgG",
      ## Covid-specific drugs
        "Hydroxychloroquine" = "hcq",
        "Kaletra" = "kaletra",
        "Azithromycin" = "az",
        "Sumamed" = "az",
        "Tocilizumab" = "tocilizumab",
        "Covalescent plasma" = "convalescent_plasma",
        "Convalescent plasma" = "convalescent_plasma",
        "Dexametason" = "dexamethasone",
        "Dexamethasone" = "dexamethasone",
        "Dexamethason" = "dexamethasone",
        "Remdesivir" = "remdesivir",
        "Favipiravir" = "favipiravir",
      ## Antibiotics
        "Amoxicillin/Clavulanic" = "amoxiclav",
        "Amoxicillin/Clavulanate" = "amoxiclav",
        "amoxiclin/klavulanat" = "amoxiclav",
        "Amoxicilin/klavulan\U00e1t" = "amoxiclav",
        "Amoksiklav" = "amoxiclav",
        "Amoxiclav" = "amoxiclav",
        "Augmentin" = "amoxiclav",
        "Piperacilin/tazobactam" = "piperacilin_tazobactam",
        "Piperacillin/Tazobactam" = "piperacilin_tazobactam",
        "Piperacilin Tazobactam" = "piperacilin_tazobactam",
        "Clarithromycin" = "clarithromycin", # Also a macrolide
        "Clarithomcin" = "clarithromycin",
        "Klarithromycin" = "clarithromycin",
        "Clacid" = "clarithromycin",
        "Klacid" = "clarithromycin",
        "Biseptol" = "cotrimoxazole",
        "sulfamethoxazol" = "sulfamethoxazol", # Combines with trimethoprim to "Cotrimoxazol"
        "trimethoprim" = "trimethoprim",
        "ofloxacin" = "ofloxacin",
        "ceftriaxone" = "ceftriaxone",
        "Ceftriaxon" = "ceftriaxone",
        "Ceftriaxom" = "ceftriaxone",
        "Ceftriaxome" = "ceftriaxone",
        "cefuroxime" = "cefuroxime",
        "cefuroxim" = "cefuroxime",
        "Cefuroxim axetim" = "cefuroxime",
        "Zinnat" = "cefuroxime",
        "Sefotak" = "cefotaxime",
        "Safotak" = "cefotaxime",
        "Taximed" = "cefotaxime",
        "cefotaxim" = "cefotaxime",
        "cefotaxime" = "cefotaxime",
        "cefepim" = "cefepime",
        "cefepime" = "cefepime",
        "Vancomycin" = "vancomycin",
        "ampicillin/sulbactam" = "ampicillin_sulbactam",
        "Ampicilin+Sulbactam" = "ampicillin_sulbactam",
        "Ciprofloksacin" = "ciprofloxacin",
        "Ciprofloxacin" = "ciprofloxacin",
        "Ciphin" = "ciprofloxacin",
        "Meropenem" = "meropenem",
        "Meronem" = "meropenem",
        "Metronidazol" = "metronidazol",
        "Linezolid" = "linezolid",
        "Doxybene" = "doxycycline",
        "Penicilin" = "penicilin",
        "Penicilin-G" = "penicilin",
        "Clindamycin" = "clindamycin",
        "Furolin" = "nitrofurantoin",
        "Cefzil" = "cefprozil",
      ## Other drugs
        "Fluconazol" = "fluconazol",
        "Flukonazol" = "fluconazol",
        "zinkorot" = "zinc",
        "Anidulafungin" = "anidulafungin",
        "Amphotericin B" = "amphotericin_b",
        "Voricoazol" = "voriconazole",
        "Aciclovir" = "aciclovir",
        "Isoprinosine" = "isoprinosine",
        "Tamiflu" = "tamiflu",
        "Lexaurin" = "lexaurin",
        "Isoprinosin" = "isoprinosin"
)

all_antibiotics <- c(
  "az",
  "amoxiclav",
  "piperacilin_tazobactam",
  "clarithromycin",
  "sulfamethoxazol",
  "trimethoprim",
  "cotrimoxazole",
  "ofloxacin",
  "ceftriaxone",
  "cefuroxime",
  "cefotaxime",
  "cefepime",
  "vancomycin",
  "ampicillin_sulbactam",
  "ciprofloxacin",
  "meropenem",
  "metronidazol",
  "linezolid",
  "doxycycline",
  "penicilin",
  "nitrofurantoin",
  "cefprozil"
)

antibiotics_not_for_pneumonia <- c(
  "nitrofurantoin",
  "metronidazol"
)

all_macrolides <- c("az", "clarithromycin")

pcr_values_positive <- c("pos","neg", "neg and pos", "susp.", "susp", "poz", "nejasn\u00fd", "nejsan\u00fd")
pcr_value_negative <- "neg"

unit_overwrites <- list(
  all = list(
    list(marker = "CRP", old_unit = "ng/l", new_unit = "mg/l")
  ),
  QKuFp = list(
    list(marker = "leukocyte_count", old_unit = "unit", new_unit = "10^9/l"),
    list(marker = "vancomycin", old_unit = "", new_unit = "g/day"),
    #TODO check the IgG with hospital
    list(marker = "IgG", old_unit = "", new_unit = "TODO")
  ),
  YqNbe = list(
    #TODO: check that d_dimer is DDU, not FEU
    #list(marker = "d_dimer", old_unit = "ng/ml DDU", new_unit = "\U03BCg/l"),
    list(marker = "d_dimer", old_unit = "ug/l", new_unit = "ng/ml DDU"),
    list(marker = "procalcitonin", old_unit = "unit", new_unit = "\U03BCg/l")
  ),
  KXqam = list(
    list(marker = "CRP", old_unit = "mg/j", new_unit = "mg/l"),
    list(marker = "d_dimer", old_unit = "ug/l", new_unit = "ng/ml DDU"),
    list(marker = "procalcitonin", old_unit = "ng/l", new_unit = "\U03BCg/l")
  ),
  ZDYlY = list(
    list(marker = "ampicillin_sulbactam", old_unit = "g", new_unit = "g/day"),
    list(marker = "cefotaxime", old_unit = "g", new_unit = "g/day"),
    list(marker = "meropenem", old_unit = "g", new_unit = "g/day"),
    list(marker = "piperacilin_tazobactam", old_unit = "g", new_unit = "g/day"),
    list(marker = "anidulafungin", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "ceftriaxone", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "ciprofloxacin", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "clarithromycin", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "dexamethasone", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "fluconazol", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "hcq", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "remdesivir", old_unit = "mg", new_unit = "mg/day"),
    list(marker = "dexamethasone", old_unit = "", new_unit = "mg/day")
  ),
  tMdnA = list(
    list(marker = "d_dimer", old_unit = "mg/l DDU", new_unit = "ng/ml DDU"),
    list(marker = "clarithromycin", old_unit = "mg/den", new_unit = "mg/day"),
    list(marker = "amoxiclav", old_unit = "mg/den", new_unit = "mg/day")
  ),
  cTdij = list(
    list(marker = "convalescent_plasma", old_unit = "", new_unit = "TU/day")
  ),
  iVMlA = list(
    list(marker = "convalescent_plasma", old_unit = "", new_unit = "TU/day")
  )



)

unit_map <- c(
  "ug/l" = "\U03BCg/l",
  "TU/d" = "TU/day",
  "g/d" = "g/day",
  "% or NT" = "%"
  )

breathing_levels <- c("AA","Oxygen", "NIPPV","MV","ECMO")
disease_levels <- c("Discharged", breathing_levels, "Death")

breathing_s_levels <- c("AA","Oxygen", "Ventilated")
disease_s_levels <- c("Discharged", breathing_s_levels, "Death")


unit_conversions <- list(
  list(markers = c(
    "cefotaxime","cefuroxime", "ampicillin_sulbactam", "amoxiclav",
    "ceftriaxone", "fluconazol", "vancomycin", "clarithromycin"),
       old_unit = "g/day", new_unit = "mg/day", mult = 1000),
  list(markers = "d_dimer", old_unit = "mg/l FEU", new_unit = "ng/ml DDU", mult = 0.5)
)
