#!/bin/bash

RESOURCES_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources"

# Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to for example your .bashrc file.
export CGE_RESFINDER_RESGENE_DB=${RESOURCES_PATH}"/db_cge/resfinder"
export CGE_RESFINDER_RESPOINT_DB=${RESOURCES_PATH}"/db_cge/pointfinder"
export CGE_DISINFINDER_DB=${RESOURCES_PATH}"/db_cge/disinfinder"
export CGE_SEROTYPEFINDER_DB=${RESOURCES_PATH}"/db_cge/serotypefinder"
export CGE_PLASMIDFINDER_DB=${RESOURCES_PATH}"/db_cge/plasmidfinder"
export CGE_RESFINDERFSA_DB=${RESOURCES_PATH}"/db_cge/kmerresistance/ResFinder"
export CGE_PMLST_DB=${RESOURCES_PATH}"/db_cge/pmlst"

export TRIMMOMATIC_ADAPTERS=${RESOURCES_PATH}"/trimmomatic/adapters"
#DATABASE_KMERFINDER="/software/resources/kmerfinder/databases/bacteria/bacteria"
DATABASE_KMERFINDER="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/kmerfinder/databases/bacteria/bacteria"


#export PERL5LIB=/software/miniconda3/envs/dgsp_amr_mlst_detection/lib/perl5/5.32


################################################################################
# Log functions
################################################################################
# Print and log messages in LOG_FILE
log_string() 
{
    echo -e "$@" | tee -a ${LOG_FILE}
}

# Run a command and log any message printed to stdout in LOG_FILE
log_command() 
{
    echo -e "$@" |  tee -a ${LOG_FILE}
    echo "$@" | bash 2>>${LOG_FILE} | tee -a ${LOG_FILE}
}




#==============================================================================
# NECESARIO MODIFICAR O TENER CONSTRUIDO
RUN="230210_LSPV002"
SRUN=$(echo $RUN | awk -F "_" '{print $2}')
INPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW"
SAMPLESHEET=${INPUT_PATH}"/"${RUN}"/"${SRUN}"/samplesheet.csv"
IR=${INPUT_PATH}"/"${RUN}"/"${SRUN}
DATE=$(date +"%y%m%d")
OUTPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/analysis/"${RUN}
CONDAPATH="/software/miniconda3/envs"
THREADS=24
SPADESMEM=72


#export LD_LIBRARY_PATH=${CONDAPATH=}/dgsp_efsa_sp/lib
# AMR MLST env



##########################################################
#
# 01 Create OUTPUT dirs
#
##########################################################

if [[ ! -e ${OUTPUT_PATH} ]]; then
    mkdir -p ${OUTPUT_PATH}"/stats/fastqc/raw"
    mkdir -p ${OUTPUT_PATH}"/stats/fastqc/trim"
    mkdir -p ${OUTPUT_PATH}"/stats/fastp"
    mkdir -p ${OUTPUT_PATH}"/stats/mash"
    mkdir -p ${OUTPUT_PATH}"/out/0_fastq"
    mkdir -p ${OUTPUT_PATH}"/logs/fastp"
    mkdir -p ${OUTPUT_PATH}"/logs/efsa"
    mkdir -p ${OUTPUT_PATH}"/logs/mash"
    mkdir -p ${OUTPUT_PATH}"/report"
    mkdir -p ${OUTPUT_PATH}"/tmp"
    elif [[ ! -d ${OUTPUT_PATH} ]]; then
    echo "${OUTPUT_PATH} already exists but is not a directory" 1>&2
fi


##########################################################
#
# ANALYSIS
#
##########################################################

# De momento generamos sample sheet
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R1_" | awk -F "_R1_" '{print $1}'
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R1_"
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R2_"

# Con el extra read, eliminamos la primera línea, en teoría los valores son reales
{
    read
    while IFS=, read -r sample fastq_1 fastq_2
    do
        echo "$sample FASTQ1: $fastq_1 FASTQ2: $fastq_2"

        ##########################################################
        #
        # 02  Assessment of the genomic sequence quality
        #
        ##########################################################

        sample="22_LM_03313"
        fastq_1="22_LM_03313_S25_R1_001.fastq.gz"
        fastq_2="22_LM_03313_S25_R2_001.fastq.gz"

        source /software/miniconda3/etc/profile.d/conda.sh
        conda activate ${CONDAPATH}/dgsp_efsa_qc

        ################################################################################
        # FastQC Quality RAW
        ################################################################################
        fastqc --quiet --threads ${THREADS} -outdir ${OUTPUT_PATH}/stats/fastqc/raw ${IR}/${fastq_1} ${IR}/${fastq_2}

        
        ################################################################################
        # FastP Quality  trimming
        ################################################################################
        fastp --thread ${THREADS} --in1  ${IR}/${fastq_1} --in2 ${IR}/${fastq_2} \
        --cut_tail --cut_window_size=10 --cut_mean_quality=20 --length_required=50 --correction \
        --json ${OUTPUT_PATH}/stats/fastp/${sample}.report.fastp.json --html ${OUTPUT_PATH}/stats/fastp/${sample}.report.fastp.html \
        --out1 ${OUTPUT_PATH}/out/0_fastq/${sample}_1.fastq.gz --out2 ${OUTPUT_PATH}/out/0_fastq/${sample}_2.fastq.gz  >> ${OUTPUT_PATH}/logs/fastp/${sample}.log 2>&1
   
        ################################################################################
        # CHECK Q30 EFSA
        ################################################################################
        # si la longitud media de lectura = (2 × 301 pb) y Q30 < 70 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 251 pb) y Q30 < 75 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 151 pb) y Q30 < 80 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 76 pb) y Q30 < 85 %, pipeline fallará debido a la baja tasa de Q30;
        # si el número total de bases es menor que el tamaño del genoma del patógeno multiplicado por el umbral mínimo de cobertura (30x), pipeline fallará debido a una cobertura insuficiente después del recorte.
        
        VAR_C1_LENGTH=301
        VAR_C1_GENOME=2880000
        # Recuperamos valor Q30

        BASESQ30=$(grep "Q30" ${OUTPUT_PATH}/logs/fastp/${sample}.log | tail -n 2 | awk -F " " '{print $3}' | awk -F "(" '{print $1}' | awk '{total += $1}END{ printf total}')
        COVQ30=$(echo 2k $BASESQ30 $VAR_C1_GENOME /p | dc)
        VALUEQ30=$(grep "Q30" ${OUTPUT_PATH}/logs/fastp/${sample}.log | tail -n 2 | awk -F " " '{print $3}' | awk -F "(" '{print $2}' | sed 's/%//g' | awk '{total += $1; count += 1}END{ printf "%4.3f\n",  total / count}')
        if((VAR_C1_LENGTH == 301)); then
            if (( $(echo "$COVQ30 < 30" | bc -l) )) || (( $(echo "$VALUEQ30 < 70" | bc -l) )); then     
                STATUSQ30="FAIL" 
            else
                STATUSQ30="PASS" 
            fi
        elif((VAR_C1_LENGTH == 251)); then
            if (( $(echo "$COVQ30 < 30" | bc -l) )) || (( $(echo "$VALUEQ30 < 75" | bc -l) )); then     
                STATUSQ30="FAIL" 
            else
                STATUSQ30="PASS" 
            fi
        elif((VAR_C1_LENGTH == 151)); then
            if (( $(echo "$COVQ30 < 30" | bc -l) )) || (( $(echo "$VALUEQ30 < 80" | bc -l) )); then     
                STATUSQ30="FAIL" 
            else
                STATUSQ30="PASS" 
            fi
        elif((VAR_C1_LENGTH == 76)); then
            if (( $(echo "$COVQ30 < 30" | bc -l) )) || (( $(echo "$VALUEQ30 < 85" | bc -l) )); then     
                STATUSQ30="FAIL" 
            else
                STATUSQ30="PASS" 
            fi
        else
            echo "Error, read lenght!"
        fi
        
        echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        echo "Q30 CHECK" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'READ_SIZE: %s\n' "$VAR_C1_LENGTH"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'GENOME_SIZE: %s\n' "$VAR_C1_GENOME"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'BASES: %s\n' "$BASESQ30"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'COV_ESTIMATED: %s\n' "$COVQ30"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'RATE Q30: %s\n' "$VALUEQ30"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
        printf 'STATUS: %s\n' "$STATUSQ30"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log

        if [ "$STATUSQ30" = "PASS" ]; then
        ################################################################################
        # FastQC Quality TRIM
        ################################################################################
        fastqc --quiet --threads ${THREADS} --outdir ${OUTPUT_PATH}/stats/fastqc/trim ${OUTPUT_PATH}/out/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/out/0_fastq/${sample}_2.fastq.gz
        conda deactivate    

        ##########################################################
        #
        # 03 Species validation
        #
        ##########################################################
        conda activate ${CONDAPATH}/dgsp_efsa_sp
        
        RefSeq88n.msh  refseq.genomes.k21s1000.msh

        #MASH_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/mash/RefSeq88n.msh"
        MASH_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/mash/refseq.genomes.k21s1000.msh"
        mash screen ${MASH_DB} ${OUTPUT_PATH}/out/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/out/0_fastq/${sample}_2.fastq.gz > ${OUTPUT_PATH}/logs/mash/${sample}_88.txt

        #nextflow run ctmrbio/BACTpipe --mashscreen_database path/to/refseq.genomes.k21s1000.msh --reads '*_R{1,2}.fastq.gz'
        nextflow run ctmrbio/BACTpipe --mashscreen_database ${MASH_DB} --reads ${OUTPUT_PATH}/out/0_fastq/${sample}_{1,2}.fastq.gz

        RefSeq88n.msh
nextflow run /ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/BACTpipe-2.7.0/bactpipe.nf --mashscreen_database ${MASH_DB} --reads '/ALMEIDA/PROJECTS/BACTERIAS/DGSP/analysis/230210_LSPV002/out/0_fastq/22_LM_03313_{1,2}.fastq.gz' --outdir p



        #https://mash.readthedocs.io/en/latest/data.html

        wget ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt 

        conda deactivate    
        else
            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            printf '%s FAILED IN Q30 CHECK\n' "$sample"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            conda deactivate       
        fi
   
        
    done
} < ${SAMPLESHEET}




