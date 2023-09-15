# Lanzar script:
# Rscript scripts/summary_to_excel.R /ALMEIDA/carpeta1/carpeta2/fichero.csv
suppressPackageStartupMessages(require(openxlsx))

tryCatch({
  
  args <- commandArgs(trailingOnly = TRUE)
  
  # Comprueba fichero de entrada
  if(length(args)==0){
    stop("¡¡ERROR, no hay fichero de entrada!!")
  }
  
  # Asigna el argumento a una variable y creamos otras de utilidad
  FILE <- args[1]
  DIROUT <- dirname(FILE)
  DATE_RES=unlist(strsplit(FILE, "_summary_"))[1]
  RUN=gsub(".csv", "", unlist(strsplit(FILE, "_summary_"))[2])
  
  # Cargamos fichero
  data <- read.table(FILE, 
                     sep = ",", 
                     dec = ".", 
                     quote = "\"",
                     header = TRUE, 
                     check.names = FALSE, 
                     stringsAsFactors = FALSE)
  
  # Comprueba si el data frame tiene 0 filas
  if(nrow(data) == 0){
    stop("El archivo está vacío.")
  }
  
  # Preparamos datos
  data$Run <- RUN
  data1 <- data[,c("Run", "Sample", "Status_Q30", "Status_SP", "Status_SNV", "Status_Contamination", "Assembly_quality")]
  data2 <- data[,c("Run", "Sample", "Sp_detected", "Cov_Q30", "Completeness", "Contamination", "contigs", "N50_(contigs)",
                   "ST", "MLST", "Antigenic_profile", "ST_patho", "cgmlst_ST", "Resistance_genes", "Antimicrobial(class)"
                   ,"Gene_mut(resistance)" )]
  
  # Creamos workbook
  wb <- createWorkbook()
  addWorksheet(wb, paste0("CONTROL_", RUN))
  addWorksheet(wb, paste0(RUN))
  
  headerStyle <- createStyle(
    fontSize = 14, fontColour = "#FFFFFF", halign = "center",
    fgFill = "#4F81BD", border = "TopBottom", borderColour = "#4F81BD"
  )
  topAlignmentStyle <- createStyle(valign="top")
  
  writeData(wb, paste0("CONTROL_", RUN), data1, startRow = 1, startCol = 1,
            headerStyle = headerStyle)
  writeData(wb, RUN, data2, startRow = 1, startCol = 1,
            headerStyle = headerStyle)
  
  setColWidths(wb, sheet = paste0("CONTROL_", RUN), cols = 1:7, widths = 22)
  setColWidths(wb, sheet = RUN, cols = 1, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 2, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 3, widths = 22)
  setColWidths(wb, sheet = RUN, cols = 4:9, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 10, widths = 55)
  setColWidths(wb, sheet = RUN, cols = 11, widths = 20)
  setColWidths(wb, sheet = RUN, cols = 12:13, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 14:16, widths = 50)
  
  # Guarda el archivo Excel
  saveWorkbook(wb, paste0(paste0(DIROUT, "/", RUN,".xlsx")), overwrite = TRUE)
  
}, error = function(e){
  message("Ocurrió un error al leer el archivo: ", e$message)
  # stop(e)  # Para detener el script
})
