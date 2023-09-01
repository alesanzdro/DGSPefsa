#!/bin/bash

# Función para mostrar ayuda
show_help() {
    echo "Uso: ./download_loci.sh [DIRECTORIO]"
    echo -e "\nDonde:"
    echo -e "  DIRECTORIO es la carpeta donde se guardarán los archivos descargados. Si no se especifica, se usará el directorio actual."
    exit 1
}

# Verificar si se solicitó ayuda
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_help
fi

# Si se especifica un directorio, úsalo; de lo contrario, usa el directorio actual
target_dir="${1:-.}"

# Crea el directorio si no existe
mkdir -p "$target_dir"

# Inicialización de variables
base_url="https://pubmlst.org/bigsdb?db=pubmlst_campylobacter_seqdef&page=downloadAlleles&locus="
total_locus=1854
locus_downloaded=0
locus_failed=0
failed_list=""

# Función para descargar un locus específico
download_locus() {
    locus="$1"
    wget "${base_url}${locus}" -O "${target_dir}/${locus}.fasta"
    # Verifica si la descarga fue exitosa
    if [ $? -eq 0 ]; then
        # Incrementa el contador de descargas exitosas
        ((locus_downloaded++))
    else
        # Incrementa el contador de descargas fallidas y añade el locus a la lista de fallos
        ((locus_failed++))
        failed_list="${failed_list}\n${locus}"
    fi
}

# Bucle principal para descargar los locus
for i in $(seq -f "CAMP%04g" 1 $total_locus); do
    download_locus $i
done

# Resumen de descargas
echo -e "\n*** RESUMEN ***"
echo -e "locus_downloaded: $locus_downloaded"
echo -e "locus_failed: $locus_failed"
if [ "$locus_failed" -gt 0 ]; then
    echo -e "list_locus_failed: $failed_list"
fi

