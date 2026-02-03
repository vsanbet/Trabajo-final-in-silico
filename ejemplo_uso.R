# ejemplo_uso.R
source("run_geyer_pipeline.R")

run_geyer_pipeline(
     path = "path/to/MaxQuant_output",
     group_regex = c("Dilution", "Thrombocytes"),
     group_names = c("Dilution", "Thrombocytes"),
     analysis_title = "Platelet-rich Plasma analysis",
     mod_file = "Oxidation (M)Sites.txt"
)
