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


export PERL5LIB=/software/miniconda3/envs/dgsp_efsa_sp/lib/perl5/5.32
export PATH=/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/signalp-5.0b/bin:$PATH
export PATH="${CONDAPATH}/dgsp_efsa_contamination/INNUca:${CONDAPATH}/dgsp_efsa_contamination/ReMatCh/ReMatCh:$PATH"

MASH_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/mash/refseq.genomes.k21s1000.msh"
#MASH_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/mash/RefSeq88n.msh"


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
BACTPIPE="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/BACTpipe-2.7.0/bactpipe.nf"
THREADS=16
SPADESMEM=72


#export LD_LIBRARY_PATH=${CONDAPATH=}/dgsp_efsa_sp/lib
# AMR MLST env



##########################################################
#
# 01 Create OUTPUT dirs
#
##########################################################

if [[ ! -e ${OUTPUT_PATH} ]]; then
    # En la carpeta pipeline es donde correremos innuca, bactpipe, confindr... los controles de calidad
    mkdir -p ${OUTPUT_PATH}"/pipeline"

    mkdir -p ${OUTPUT_PATH}"/reports/fastqc/raw"
    mkdir -p ${OUTPUT_PATH}"/reports/fastqc/trim"
    mkdir -p ${OUTPUT_PATH}"/reports/fastp"

    mkdir -p ${OUTPUT_PATH}"/logs/fastp"
    mkdir -p ${OUTPUT_PATH}"/logs/efsa"
    mkdir -p ${OUTPUT_PATH}"/logs/BACTpipe"
    mkdir -p ${OUTPUT_PATH}"/logs/innuca"
    mkdir -p ${OUTPUT_PATH}"/logs/confindr"
    mkdir -p ${OUTPUT_PATH}"/logs/checkm"

    innuca
   
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


        # VARIABLES
        VAR_C1_LENGTH=301
        VAR_C1_GENOME=2880000
        VAR_C2_SPE="Listeria monocytogenes"
        VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)


        # DIRECTORIOS ASOCIADOS A MUESTRA / ESPECIE
        mkdir -p ${OUTPUT_PATH}"/pipeline/"${sample}"/fastq"
        mkdir -p ${OUTPUT_PATH}"/pipeline/"${sample}"/assembly"
        
        if [[ "$VAR_C2_SPE" == "Salmonella enterica" ]];then
            #mkdir -p ${OUTPUT_PATH}"/logs/seqsero2"
            #mkdir -p ${OUTPUT_PATH}"/pipeline/"${sample}"/seqsero2"
        elif [[ "$VAR_C2_SPE" == "Listeria monocytogenes" ]]; then
            mkdir -p ${OUTPUT_PATH}"/logs/lissero"
            #mkdir -p ${OUTPUT_PATH}"/pipeline/"${sample}"/lissero"
        fi

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
        fastqc --quiet --threads ${THREADS} -outdir ${OUTPUT_PATH}/reports/fastqc/raw ${IR}/${fastq_1} ${IR}/${fastq_2}

        
        ################################################################################
        # FastP Quality  trimming
        ################################################################################
        fastp --thread ${THREADS} --in1  ${IR}/${fastq_1} --in2 ${IR}/${fastq_2} \
        --cut_tail --cut_window_size=10 --cut_mean_quality=20 --length_required=50 --correction \
        --json ${OUTPUT_PATH}/reports/fastp/${sample}.report.fastp.json --html ${OUTPUT_PATH}/reports/fastp/${sample}.report.fastp.html \
        --out1 ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_1.fastq.gz --out2 ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_2.fastq.gz  >> ${OUTPUT_PATH}/logs/fastp/${sample}.log 2>&1
   
        ################################################################################
        # CHECK Q30 EFSA
        ################################################################################
        # si la longitud media de lectura = (2 × 301 pb) y Q30 < 70 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 251 pb) y Q30 < 75 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 151 pb) y Q30 < 80 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 76 pb) y Q30 < 85 %, pipeline fallará debido a la baja tasa de Q30;
        # si el número total de bases es menor que el tamaño del genoma del patógeno multiplicado por el umbral mínimo de cobertura (30x), pipeline fallará debido a una cobertura insuficiente después del recorte.
        


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
            fastqc --quiet --threads ${THREADS} --outdir ${OUTPUT_PATH}/reports/fastqc/trim ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_2.fastq.gz
            conda deactivate    

            ##########################################################
            #
            # 03 Species validation
            #
            ##########################################################
            conda activate ${CONDAPATH}/dgsp_efsa_sp

            #mash screen ${MASH_DB} ${OUTPUT_PATH}/out/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/out/0_fastq/${sample}_2.fastq.gz > ${OUTPUT_PATH}/logs/mash/${sample}_88.txt
       
            cd ${OUTPUT_PATH}/pipeline/${sample}
            nextflow-20.10.0-all run ${BACTPIPE} \
            --mashscreen_database ${MASH_DB} \
            --reads './fastq/*_{1,2}.fastq.gz' >> ${OUTPUT_PATH}/logs/BACTpipe/${sample}.log 2>&1

            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            echo "SPECIES CHECK" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            cat ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv >> ${OUTPUT_PATH}/logs/efsa/${sample}.log
            C2_STATUS=$(cat ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk '{print $2}') 
            C2_SPE=$(cat ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk -F "\t" '{print $5}' | sed -e "s/\[\|\]\|'//g") 

            printf 'OUR SAMPLE: %s\n' "$VAR_C2_SPE"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log


            if [[ "$C2_STATUS" = "PASS" ]]  && [[ "$C2_SPE" == "$VAR_C2_SPE" ]]; then     
                STATUSBACT="PASS"
                printf 'STATUS: %s\n' "$STATUSBACT"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log

                # CHECK SEROVARES !!!!!
                if [[ "$VAR_C2_SPE" = "Salmonella enterica" ]];then
                    # Raw reads k-mer ("-m k"), for separated paired-end raw reads ("-t 2")
                    SeqSero2_package.py -m k -t 2 -i ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_2.fastq.gz >> ${OUTPUT_PATH}/pipeline/${sample}/${sample}.seqsero2.txt 2>&1
                fi

                if [[ "$VAR_C2_SPE" = "Listeria monocytogenes" ]];then
                    (cd ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill && lissero ${sample}.contigs.fa >> ${OUTPUT_PATH}/pipeline/${sample}/${sample}.lissero.txt 2>&1)
                    var2erase=$(echo ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/)
                    #IMPORTANTE EN SED SE EMPLEA ":" PARA SEPARAS CAMPOS Y PODER ELIMINAR PATH PARCIAL
                    lissero ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa --logfile ${OUTPUT_PATH}/logs/lissero/${sample}.log | sed "s:^$var2erase::" > ${OUTPUT_PATH}/pipeline/${sample}/${sample}.lissero.txt 
                fi

                ##########################################################
                #
                # 04 Contamination
                #
                ##########################################################
                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                echo "Contamination CHECK" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                confindr.py -i ${OUTPUT_PATH}/pipeline/${sample}/fastq -o ConFindr >> ${OUTPUT_PATH}/logs/confindr/${sample}.log 2>&1
                conda deactivate

                conda activate ${CONDAPATH}/dgsp_efsa_contamination

                INNUca.py -j ${THREADS} \
                -i ${OUTPUT_PATH}/pipeline/${sample}/fastq \
                -o INNUca \
                -s "$VAR_C2_SPE" -g $VAR_C3_GENOME >> ${OUTPUT_PATH}/logs/innuca/${sample}.log 2>&1

                # Recogemos resultados contaminación
                C3_confindr=$(awk -F "," '{print $4}' ${OUTPUT_PATH}/pipeline/${sample}/ConFindr/confindr_report.csv | grep -v "ContamStatus" | uniq)
                C3_warning=$(awk -F "," '{print $3}' ${OUTPUT_PATH}/pipeline/${sample}/ConFindr/confindr_report.csv | grep -v "NumContamSNVs" | awk '{total += $1; count += 1}END{ printf "%4.3f\n",  total / count}')
                
                #Si el número de SNV contaminados supera los siguientes umbrales específicos de especies 
                #(como lo sugiere Deneke et al., 2021a), se envía un mensaje de advertencia:
                #• S. enterica: 7
                #• L. monocytogenes: 3
                #• E.coli: 4
                if [[ "$VAR_C2_SPE" = "Listeria monocytogenes" ]];then
                    if (( $(echo "$C3_warning > 3" | bc -l) ));then
                        CONTSNV="FAIL" 
                    else
                        CONTSNV="PASS" 
                    fi
                elif [[ "$VAR_C2_SPE" = "Salmonella enterica" ]];then
                    if (( $(echo "$C3_warning > 7" | bc -l) )); then
                        CONTSNV="FAIL" 
                    else
                        CONTSNV="PASS" 
                    fi
                elif [[ "$VAR_C2_SPE" = "Escherichia coli" ]];then
                    if (( $(echo "$C3_warning > 4" | bc -l) )); then 
                        CONTSNV="FAIL" 
                    else
                        CONTSNV="PASS" 
                    fi
                else
                    if (( $(echo "$C3_warning > 4" | bc -l) )); then     
                        CONTSNV="FAIL" 
                    else
                        CONTSNV="PASS" 
                    fi
                fi


                if [[ "$CONTSNV" = "FAIL" ]]; then 
                    echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    echo "========================================" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    echo "WARNING!" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    echo "ConFindr show SNV contamination" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    echo "========================================" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                else

                repinnuca=$(ls ${OUTPUT_PATH}/pipeline/${sample}/INNUca/samples_report.*)
                C3_innuca=$(awk -F "\t" '{print $3}' $repinnuca  | grep -v "samples_passQC")

                if [[ "$C3_confindr" = "False" ]]  && [[ "$C3_innuca" = "PASS" ]]; then     
                    STATUSCONT="PASS"
                    printf 'STATUS: %s\n' "$STATUSCONT"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log

                ##########################################################
                #
                # 05 Assembly
                #
                ##########################################################
                conda activate ${CONDAPATH}/dgsp_efsa_sp
                
                shovill --depth 100 --kmers 31,33,55,77,99,127 --minlen 500 --cpus ${THREADS} \
                --R1 ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_1.fastq.gz \
                --R2 ${OUTPUT_PATH}/pipeline/${sample}/fastq/${sample}_1.fastq.gz \
                --outdir ${OUTPUT_PATH}"/pipeline/"${sample}"/assembly/shovill"
                
                cp ${OUTPUT_PATH}/pipeline/${sample}/assembly/shovill/contigs.fa ${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.contigs.fa
                statswrapper.sh in=${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.contigs.fa > ${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.assembly_stats.txt

                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                echo "Assembly quality CHECK" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                mkdir -p ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs
                #Desglosamos fichero contigs en fastas individuales
                #(cd ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs && cat ../${sample}.contigs.fa | awk '{
                (cd ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs && cat ../../BACTpipe_results/shovill/${sample}.contigs.fa | awk '{
                        if (substr($0, 1, 1)==">") {filename=(substr($1,2) ".fa")}
                        print $0 >> filename
                        close(filename)
                }'
                )

                checkm lineage_wf -t ${THREADS} \
                -x fa ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs \
                ${OUTPUT_PATH}/pipeline/${sample}/checkm >> ${OUTPUT_PATH}/logs/checkm/${sample}.log 2>&1
            
                #----------------------------------------------------------------------------------------------------------------------------------------------------------
                #  genes    c__Bacilli (UID354)      515         328           182        1   314   12   1   0   0       99.45            5.22              86.67          
                #----------------------------------------------------------------------------------------------------------------------------------------------------------
             
                #revisar /software/miniconda3/envs/dgsp_efsa_sp/checkm_data
 
                cd /ALMEIDA/PROJECTS/BACTERIAS/DGSP/analysis/230210_LSPV002/pipeline/22_LM_03313/checkm_test
                checkm lineage_wf -t ${THREADS} -x fa bins/22_LM_03313/ salida
                else
                    echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    printf '%s FAILED IN CONTAMINATION CHECK\n' "$sample"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                    conda deactivate
                fi



            else
                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                printf '%s FAILED IN SPECIES CHECK\n' "$sample"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
                conda deactivate
            fi

        else
            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            printf '%s FAILED IN Q30 CHECK\n' "$sample"  | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            echo "****************************************" | tee -a ${OUTPUT_PATH}/logs/efsa/${sample}.log
            conda deactivate       
        fi
   
        
    done
} < ${SAMPLESHEET}




