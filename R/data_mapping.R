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
      ## Covid-specific drugs
        "Hydroxychloroquine" = "hcq",
        "Kaletra" = "kaletra",
        "Azithromycin" = "az",
        "Tocilizumab" = "tocilizumab",
        "Covalescent plasma" = "convalescent_plasma",
      ## Antibiotics
        "Amoxicillin/Clavulanic" = "amoxiclav",
        "Amoxicillin/Clavulanate" = "amoxiclav",
        "Augmentin" = "amoxiclav",
        "Piperacilin/tazobactam" = "piperacilin_tazobactam",
        "Piperacillin/Tazobactam" = "piperacilin_tazobactam",
        "Clarithromycin" = "clarithromycin",
        "ofloxacin" = "ofloxacin",
        "ceftriaxone" = "ceftriaxone",
        "cefuroxime" = "cefuroxime",
        "cefuroxim" = "cefuroxime",
        "cefotaxim" = "cefotaxime",
        "cefotaxime" = "cefotaxime",
        "cefepim" = "cefepime",
        "cefepime" = "cefepime",
        "Vancomycin" = "vancomycin",
        "ampicillin/sulbactam" = "ampicillin_sulbactam",
        "trimethoprim" = "trimethoprim",
        "Ciprofloxacin" = "ciprofloxacin",
        "Meropenem" = "meropenem",
        "Metronidazol" = "Metronidazol",
      ## Other drugs
        "Fluconazol" = "fluconazol",
        "zinkorot" = "zinc",
        "Anidulafungin" = "anidulafungin",
        "sulfamethoxazol" = "sulfamethoxazol",
        "Isoprinosine" = "isoprinosine",
        "Tamiflu" = "tamiflu",
        "Lexaurin" = "lexaurin"
)

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
    #list(marker = "d_dimer", old_unit = "ng/ml DDU", new_unit = "\U03BCg/l"),
    list(marker = "d_dimer", old_unit = "ug/l", new_unit = "ng/ml DDU"),
    list(marker = "procalcitonin", old_unit = "unit", new_unit = "\U03BCg/l")
  )

)

unit_map <- c(
  "ug/l" = "\U03BCg/l",
  "TU/d" = "TU/day",
  "g/d" = "g/day",
  "% or NT" = "%"
  )


unit_conversions <- list(
  list(markers = c("cefotaxime","cefuroxime", "ampicillin_sulbactam", "amoxiclav"),
       old_unit = "g/day", new_unit = "mg/day", mult = 1000)
)
