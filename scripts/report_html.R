# Rscript scripts/report_html.R /ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/scripts/RECOPILACIÓN_RESULTADOS/


suppressPackageStartupMessages(require(ReportingTools))

tryCatch({
  
  args <- commandArgs(trailingOnly = TRUE)
  
  # Comprueba fichero de entrada
  if(length(args)==0){
    stop("¡¡ERROR, no hay directorio de entrada!!")
  }
  INPUT_PATH=args[1]
  FILES=dir(INPUT_PATH)
  FILES=FILES[grepl("_summary_", FILES)]
  if(length(FILES)>1){
    DIROUT <- dirname(FILES[1])
    for(i in 1:length(FILES)){
      # i=1
      DATE_RES=unlist(strsplit(FILES[i], "_summary_"))[1]
      RUN=gsub(".csv", "", unlist(strsplit(FILES[i], "_summary_"))[2])
      
      data <- read.table(paste0(INPUT_PATH, FILES[i]), 
                         sep = ",", 
                         dec = ".", 
                         quote = "\"",
                         header = T, 
                         check.names = F, 
                         stringsAsFactors = F)
      data$Run <- RUN
      data$Date_results <- DATE_RES
      datar <- data[,c("Run", "Sample", "Sp_detected", "Cov_Q30", "Completeness", "Contamination", "contigs", "N50_(contigs)",
                       "ST", "MLST", "Antigenic_profile", "ST_patho", "cgmlst_ST", "Resistance_genes", "Antimicrobial(class)"
                       ,"Gene_mut(resistance)", "Status_Q30", "Status_SP", "Status_SNV", "Status_Contamination", "Assembly_quality",
                       "Date_results")]
      if(i == 1){
        sdata <- datar
      }else{
        sdata <- rbind(sdata, datar)
      }
    }
    DATE=format(Sys.time(), '%y%m%d')
    DITOUT_HTML <- paste0(DIROUT, "/reportHTML")
    
    if(!file.exists(DITOUT_HTML)) {
      dir.create(DITOUT_HTML)
    }
    htmlRep <- HTMLReport(shortName = paste0(DATE, '_DGSP.html'),
                          title = paste0('RESULTADOS ANÁLISIS DGSP EFSA  <br> (última actualización ', format(Sys.time(), '%d %B, %Y'), ")"),
                          reportDirectory = DITOUT_HTML)
    
    publish(sdata, htmlRep)
    finish(htmlRep)
  }else{
    "NO 'SUMMARIES' FOUND!"  
  }
  
}, error = function(e){
  message("Ocurrió un error al leer el archivo: ", e$message)
  # stop(e)  # Para detener el script
})










 