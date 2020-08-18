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
      ## Covid-specific drugs
        "Hydroxychloroquine" = "hcq",
        "Kaletra" = "kaletra",
        "Azithromycin" = "az",
        "Tocilizumab" = "tocilizumab",
        "Covalescent plasma" = "convalescent_plasma",
      ## Antibiotics
        "Amoxicillin/Clavulanate" = "amoxiclav",
        "Augmentin" = "amoxiclav",
        "Piperacilin/tazobactam" = "piperacilin_tazobactam",
        "Clarithromycin" = "clarithromycin",
        "ofloxacin" = "ofloxacin",
        "ceftriaxone" = "ceftriaxone",
        "cefuroxim" = "cefuroxim",
        "cefotaxim" = "cefotaxim",
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
  QKuFp = list(
    list(marker = "leukocyte_count", old_unit = "unit", new_unit = "10^9/l")
    ),
  YqNbe = list(
    list(marker = "d_dimer", old_unit = "ng/ml DDU", new_unit = "μg/l"),
    list(marker = "CRP", old_unit = "ng/l", new_unit = "mg/l")
   )
)
