setwd("/ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/scripts")

FILES = dir()
FILES = FILES[grepl("_summary_", FILES)]

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

  
  
}else{
  "NO 'SUMMARIES' FOUND!"  
}


