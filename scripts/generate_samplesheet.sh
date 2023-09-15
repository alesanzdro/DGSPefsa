#!/bin/bash

# uso:
# generate-samplesheet.sh <dir/fastq>

# Comprueba si se ha proporcionado un argumento
if [ -z "$1" ]; then
  echo "Por favor, proporciona una ruta de carpeta como argumento."
  exit 1
fi

# Cambia al directorio proporcionado
cd "$1" || exit

# Crea un archivo de salida y escribe las cabeceras de las columnas
echo "sample,fastq_1,fastq_2" > samplesheet.csv

# Recorre todos los archivos R1 y busca su correspondiente archivo R2
for fq1 in *_R1_001.fastq.gz; do
  sample=$(echo "$fq1" | awk -F '_' '{print $1"_"$2"_"$3}')
  fq2=$(echo "$fq1" | sed 's/_R1_/_R2_/')
  
  # Añade la información al archivo de salida
  if [ -e "$fq2" ]; then
    echo "$sample,$fq1,$fq2" >> samplesheet.csv
  fi
done
