setwd("/ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/scripts")

FILES = dir()
FILES = FILES[grepl("_summary_", FILES)]

library(openxlsx)


if(length(FILES)>1){
  for(i in 1:length(FILES)){
    # i=1
    DATE_RES=unlist(strsplit(FILES[i], "_summary_"))[1]
    RUN=gsub(".csv", "", unlist(strsplit(FILES[i], "_summary_"))[2])
    
    data <- read.table(FILES[i], 
                       sep = ",", 
                       dec = ".", 
                       quote = "\"",
                       header = T, 
                       check.names = F, 
                       stringsAsFactors = F)
    
  }
  data$Run <- RUN
  
  data1 <- data[,c("Run", "Sample", "Status_Q30", "Status_SP", "Status_SNV", "Status_Contamination", "Assembly_quality")]
  data2 <- data[,c("Run", "Sample", "Sp_detected", "Cov_Q30", "Completeness", "Contamination", "contigs", "N50_(contigs)",
                   "ST", "MLST", "Antigenic_profile", "ST_patho", "cgmlst_ST", "Resistance_genes", "Antimicrobial(class)"
                   ,"Gene_mut(resistance)" )]
  # Crea un nuevo archivo Excel
  if(i==1){
    wb <- createWorkbook()
  }
  # Agrega tu dataframe al archivo Excel
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
  
  # Fusiona celdas para la columna "WHO"
  
  # mergeCells(wb, SAMPLE, cols = 3, rows = 2:5)   # "VOI (Variants of Interest)"
  # mergeCells(wb, SAMPLE, cols = 3, rows = 6:8)   # "VUM (Variants under Monitoring)"
  # mergeCells(wb, SAMPLE, cols = 3, rows = 9:10)  # "Sublinajes de interés"
  # mergeCells(wb, SAMPLE, cols = 3, rows = 11:dim(varfinal)[1]+1) # "Otros"
  # addStyle(wb, SAMPLE, style = topAlignmentStyle, cols = 3, rows = 2:5)
  # addStyle(wb, SAMPLE, style = topAlignmentStyle, cols = 3, rows = 6:8)
  # addStyle(wb, SAMPLE, style = topAlignmentStyle, cols = 3, rows = 9:10) 
  # addStyle(wb, SAMPLE, style = topAlignmentStyle, cols = 3, rows = 11:dim(varfinal)[1]+1)
  
  setColWidths(wb, sheet = paste0("CONTROL_", RUN), cols = 1:7, widths = 22)
  setColWidths(wb, sheet = RUN, cols = 1, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 2, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 3, widths = 22)
  setColWidths(wb, sheet = RUN, cols = 4:9, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 10, widths = 55)
  setColWidths(wb, sheet = RUN, cols = 11, widths = 20)
  setColWidths(wb, sheet = RUN, cols = 12:13, widths = 18)
  setColWidths(wb, sheet = RUN, cols = 14:16, widths = 50)
  
  
}

# Guarda el archivo Excel
saveWorkbook(wb, paste0(RUN,".xlsx"), overwrite = TRUE)


detach("package:openxlsx", unload = TRUE)


  
  
}else{
  "NO 'SUMMARIES' FOUND!"  
}


