dummy_markers <- c("Drugs", "Other marker (optional)", "Other Covid drug", "Yet another drug")


marker_map <- c(
      ## Markers
        "PCR" = "pcr_value",
        "Supp. O2" = "oxygen_flow",
        "SpO2" = "SpO2",
        "SpO2 native" = "SpO2_native",
        "Horowitz index" = "Horowitz_index",
        "PEEP" = "PEEP",
        "Ferritin" = "ferritin",
        "D-dimer" = "d_dimer",
        "CRP" = "CRP",
        "IL-6" = "IL_6",
        "PCT" = "procalcitonin",
        "Lymphocyte count" = "lymphocyte_count",
        "Leukocyte count" = "leukocyte_count",
        "Leukocytes" = "leukocyte_count",
        "Leukocyte" = "leukocyte_count",
        "urea" = "urea",
        "IgG" = "IgG",
        "ARO" = "ICU", # We know the patient was on ICU, but no details on breathing
      ## Covid-specific drugs
        "Hydroxychloroquine" = "hcq",
        "Kaletra" = "kaletra",
        "Azithromycin" = "az",
        "Tocilizumab" = "tocilizumab",
        "Covalescent plasma" = "convalescent_plasma",
        "Dexametason" = "dexamethasone",
        "Dexamethasone" = "dexamethasone",
        "Remdesivir" = "remdesivir",
      ## Antibiotics
        "Amoxicillin/Clavulanic" = "amoxiclav",
        "Amoxicillin/Clavulanate" = "amoxiclav",
        "amoxiclin/klavulanat" = "amoxiclav",
        "Augmentin" = "amoxiclav",
        "Piperacilin/tazobactam" = "piperacilin_tazobactam",
        "Piperacillin/Tazobactam" = "piperacilin_tazobactam",
        "Piperacilin Tazobactam" = "piperacilin_tazobactam",
        "Clarithromycin" = "clarithromycin", # Also a macrolide
        "Clarithomcin" = "clarithromycin",
        "sulfamethoxazol" = "sulfamethoxazol", # Combines with trimethoprim to "Cotrimoxazol"
        "trimethoprim" = "trimethoprim",
        "ofloxacin" = "ofloxacin",
        "ceftriaxone" = "ceftriaxone",
        "Ceftriaxon" = "ceftriaxone",
        "cefuroxime" = "cefuroxime",
        "cefuroxim" = "cefuroxime",
        "cefotaxim" = "cefotaxime",
        "cefotaxime" = "cefotaxime",
        "cefepim" = "cefepime",
        "cefepime" = "cefepime",
        "Vancomycin" = "vancomycin",
        "ampicillin/sulbactam" = "ampicillin_sulbactam",
        "Ampicilin+Sulbactam" = "ampicillin_sulbactam",
        "Ciprofloxacin" = "ciprofloxacin",
        "Meropenem" = "meropenem",
        "Metronidazol" = "metronidazol",
        "Linezolid" = "linezolid",
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
        "Lexaurin" = "lexaurin"
)

all_antibiotics <- c(
  "az",
  "amoxiclav",
  "piperacilin_tazobactam",
  "clarithromycin",
  "sulfamethoxazol",
  "trimethoprim",
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
  "linezolid"
)

all_macrolides <- c("az", "clarithromycin")

pcr_values_positive <- c("pos","neg", "neg and pos", "susp.", "susp", "poz")
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

unit_conversions <- list(
  list(markers = c("cefotaxime","cefuroxime", "ampicillin_sulbactam", "amoxiclav", "ceftriaxone", "fluconazol"),
       old_unit = "g/day", new_unit = "mg/day", mult = 1000)
)
