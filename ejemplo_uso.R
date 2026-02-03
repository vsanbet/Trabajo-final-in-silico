# ejemplo_uso.R
source("run_geyer_pipeline.R")

# Obtención de grupos
run_geyer_pipeline(
     path = "ruta/a/tu/dataset", 
     analysis_title = "QC check", 
     inspect_samples = TRUE
)

# Análisis de ficheros
run_geyer_pipeline(
     path = "path/to/MaxQuant_output",
     group_regex = c("Dilution", "Thrombocytes"),
     group_names = c("Dilution", "Thrombocytes"),
     analysis_title = "Platelet-rich Plasma analysis",
     mod_file = "Oxidation (M)Sites.txt"
)
