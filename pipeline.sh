#!/bin/bash

# RESOURCES_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources"

# # Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to for example your .bashrc file.
# export CGE_RESFINDER_RESGENE_DB=${RESOURCES_PATH}"/db_cge/resfinder"
# export CGE_RESFINDER_RESPOINT_DB=${RESOURCES_PATH}"/db_cge/pointfinder"
# export CGE_DISINFINDER_DB=${RESOURCES_PATH}"/db_cge/disinfinder"
# export CGE_SEROTYPEFINDER_DB=${RESOURCES_PATH}"/db_cge/serotypefinder"
# export CGE_PLASMIDFINDER_DB=${RESOURCES_PATH}"/db_cge/plasmidfinder"
# export CGE_RESFINDERFSA_DB=${RESOURCES_PATH}"/db_cge/kmerresistance/ResFinder"
# export CGE_PMLST_DB=${RESOURCES_PATH}"/db_cge/pmlst"

# export TRIMMOMATIC_ADAPTERS=${RESOURCES_PATH}"/trimmomatic/adapters"
# #DATABASE_KMERFINDER="/software/resources/kmerfinder/databases/bacteria/bacteria"
# DATABASE_KMERFINDER="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/kmerfinder/databases/bacteria/bacteria"

# Creamos enlace simbólico carpeta resources
#ln -sf /ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/ /software/resources/dgsp
#unlink /software/resources/dgsp


################################################################################
# Log functions
################################################################################
# Print and log messages in LOG_FILE
# log_string() {
#     echo -e "$@" | tee -a ${LOG_FILE}
# }

# # Run a command and log any message printed to stdout in LOG_FILE
# log_command() {
#     echo -e "$@" | tee -a ${LOG_FILE}
#     echo "$@" | bash 2>>${LOG_FILE} | tee -a ${LOG_FILE}
# }


# REVISAR EJEMPLO USO DOCKER
#https://onestopdataanalysis.com/checkm-completeness-contamination/
#==============================================================================
# NECESARIO MODIFICAR O TENER CONSTRUIDO
#RUN="230331_LSPV005"
RUN=$1
SRUN=$(echo $RUN | awk -F "_" '{print $2}')

INPUT_PATH="/home/susana/DGSP/RAW"
OUTPUT_PATH="/home/susana/DGSP/analysis_efsa/"${RUN}
RESOURCES="/software/resources"
EFSA_PATH="/software/DGSPefsa"
THREADS=14

#INPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/sample_data/RAW"
#OUTPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/sample_data/ANALYSIS/"${RUN}
#RESOURCES="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/DGSPefsa/resources"
#THREADS=20

CONDAPATH="/software/miniconda3/envs"

SAMPLESHEET=${INPUT_PATH}/${RUN}/${SRUN}/samplesheet.csv
IR=${INPUT_PATH}/${RUN}/${SRUN}
DATE=$(date +"%y%m%d")
#DATE="230919"

VAR_C1_LENGTH=301
VAR_C4_COMPLETENESS=98
VAR_C4_CONTAMINATION=2

PATHARIBA=${RESOURCES}/ariba
PATHcgMLST=${RESOURCES}/cgMLST_data
PATHCONFINDR=${RESOURCES}/confindr_db
PATHMASH=${RESOURCES}/mash/refseq.genomes.k21s1000.msh
# PATHMLST=${RESOURCES}/pubmlst

export PATH=${RESOURCES}/signalp-5.0b/bin:$PATH
export PERL5LIB=/software/miniconda3/envs/dgsp_efsa_sp/lib/perl5/5.32
export PATH="${CONDAPATH}/dgsp_efsa_contamination/INNUca:${CONDAPATH}/dgsp_efsa_contamination/ReMatCh/ReMatCh:$PATH"


#export LD_LIBRARY_PATH=${CONDAPATH=}/dgsp_efsa_sp/lib
# AMR MLST env

##########################################################
#
# 01 Create OUTPUT dirs
#
##########################################################

if [[ ! -e ${OUTPUT_PATH} ]]; then

    # Carpetas principales
    mkdir -p ${OUTPUT_PATH}/0_fastq
    mkdir -p ${OUTPUT_PATH}/1_qc
    mkdir -p ${OUTPUT_PATH}/2_assembly
    mkdir -p ${OUTPUT_PATH}/3_typing
    mkdir -p ${OUTPUT_PATH}/4_results

    # La carpeta pipeline es donde correremos innuca, bactpipe, confindr... los controles de calidad
    mkdir -p ${OUTPUT_PATH}/tmp/pipeline

    # Carpeta para recoger múltiples reportes
    mkdir -p ${OUTPUT_PATH}/1_qc/fastqc_raw
    mkdir -p ${OUTPUT_PATH}/1_qc/fastqc_trim
    mkdir -p ${OUTPUT_PATH}/1_qc/fastp

    mkdir -p ${OUTPUT_PATH}/log/fastp
    mkdir -p ${OUTPUT_PATH}/log/efsa
    mkdir -p ${OUTPUT_PATH}/log/BACTpipe
    mkdir -p ${OUTPUT_PATH}/log/innuca
    mkdir -p ${OUTPUT_PATH}/log/confindr
    mkdir -p ${OUTPUT_PATH}/log/checkm
    mkdir -p ${OUTPUT_PATH}/log/seq_typing
    mkdir -p ${OUTPUT_PATH}/log/patho_typing
    mkdir -p ${OUTPUT_PATH}/log/chewbbaca
    mkdir -p ${OUTPUT_PATH}/log/resfinder
    mkdir -p ${OUTPUT_PATH}/log/mlst
    mkdir -p ${OUTPUT_PATH}/log/sistr
    mkdir -p ${OUTPUT_PATH}/log/abricate
    mkdir -p ${OUTPUT_PATH}/log/ariba
    mkdir -p ${OUTPUT_PATH}/log/lissero

elif [[ ! -d ${OUTPUT_PATH} ]]; then
    echo "${OUTPUT_PATH} already exists but is not a directory" 1>&2
fi

#echo -e "Sample,Fq1,Fq2,Expected_sp,Genome_size,Min_genome_size,Max_genome_size,Max_SNV,\
#Max_Contigs,Read_size,Bases_Q30,Cov_Q30,Rate_Q30,Status_Q30,Sp_detected,Status_SP,SNV_detected,\
#Status_SNV,Expected_completeness,Expected_contamination,Marker_lineage,Completeness,Contamination,\
#Strain_heterogeneity,Taxonomy_(contained),Genome_size_(Mbp),Gene_count,out_genome_size_(Mbp)_mean,\
#out_genome_size_(Mbp)_std,out_gene_count_mean,out_gene_count_std,contigs,N50_(contigs),\
#Status_Contamination,Assembly_quality,ST,MLST,Antigenic_profile,ST_patho,cgmlst_ST,h1,h2,o_antigen,serovar,Resistance_genes,\
#Antimicrobial(class),Gene_mut(resistance)"  > "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"

echo -e "\"Sample\",\"Fq1\",\"Fq2\",\"Expected_sp\",\"Genome_size\",\"Min_genome_size\",\"Max_genome_size\",\"Max_SNV\",\
\"Max_Contigs\",\"Read_size\",\"Bases_Q30\",\"Cov_Q30\",\"Rate_Q30\",\"Status_Q30\",\"Sp_detected\",\"Status_SP\",\"SNV_detected\",\
\"Status_SNV\",\"Expected_completeness\",\"Expected_contamination\",\"Marker_lineage\",\"Completeness\",\"Contamination\",\
\"Strain_heterogeneity\",\"Taxonomy_(contained)\",\"Genome_size_(Mbp)\",\"Gene_count\",\"out_genome_size_(Mbp)_mean\",\
\"out_genome_size_(Mbp)_std\",\"out_gene_count_mean\",\"out_gene_count_std\",\"contigs\",\"N50_(contigs)\",\
\"Status_Contamination\",\"Assembly_quality\",\"ST\",\"MLST\",\"Antigenic_profile\",\"ST_patho\",\"cgmlst_ST\",\"h1\",\"h2\",\"o_antigen\",\"serovar\",\"Resistance_genes\",\
\"Antimicrobial(class)\",\"Gene_mut(resistance)\"" > "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"

##########################################################
#
# ANALYSIS
#
##########################################################

# Con el extra read, eliminamos la primera línea, en teoría los valores son reales
{
    read -r
    while IFS=, read -r sample fastq_1 fastq_2; do
        echo "$sample FASTQ1: $fastq_1 FASTQ2: $fastq_2"

        # sample="17_STEC_10960"
        # fastq_1="17_STEC_10960_S65_R1_001.fastq.gz"
        # fastq_2="17_STEC_10960_S65_R2_001.fastq.gz"
        # sample="22_LMON_07334"
        # fastq_1="22_LMON_07334_S49_R1_001.fastq.gz"
        # fastq_2="22_LMON_07334_S49_R2_001.fastq.gz"
        # sample="22_SALM_01804"
        # fastq_1="22_SALM_01804_S67_R1_001.fastq.gz"
        # fastq_2="22_SALM_01804_S67_R2_001.fastq.gz"
        # sample="23_CAMP_01451"
        # fastq_1="23_CAMP_01451_S43_R1_001.fastq.gz"
        # fastq_2="23_CAMP_01451_S43_R1_001.fastq.gz"
        
        #sample="22_SALM_10116"
        #fastq_1="22_SALM_10116_S76_R1_001.fastq.gz"
        #fastq_2="22_SALM_10116_S76_R2_001.fastq.gz"
        
        #sample="22_SALM_10116"
        #fastq_1="22_SALM_10116_S76_R1_001.fastq.gz"
        #fastq_2="22_SALM_10116_S76_R2_001.fastq.gz"

        #22_SALM_15935	22_SALM_15935_S17_R1_001.fastq.gz	22_SALM_15935_S17_R2_001.fastq.gz



        SPE=$(echo $fastq_1 | awk -F "_" '{print $2}')

        if [[ "$SPE" == "SAL" ]] || [[ "$SPE" == "SALM" ]]; then
            VAR_C1_GENOME=5000000
            VAR_C2_SPE="Salmonella enterica"
            VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)
            VAR_C3_SNV=7
            VAR_C4_GENOME_MIN=4.3
            VAR_C4_GENOME_MAX=5.3
            VAR_C4_CONTIGS=300
        elif [[ "$SPE" == "LM" ]] || [[ "$SPE" == "LMON" ]]; then
            VAR_C1_GENOME=2880000
            VAR_C2_SPE="Listeria monocytogenes"
            VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)
            # Modificamos de 3 a 10, para hacer más laxo el filtro
	    VAR_C3_SNV=10
            VAR_C4_GENOME_MIN=2.7
            VAR_C4_GENOME_MAX=3.2
            VAR_C4_CONTIGS=200
        elif [[ "$SPE" == "STEC" ]] || [[ "$SPE" == "ECOL" ]]; then
            VAR_C1_GENOME=5000000
            VAR_C2_SPE="Escherichia coli"
            VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)
            VAR_C3_SNV=7
            VAR_C4_GENOME_MIN=4.5
            VAR_C4_GENOME_MAX=5.9
            VAR_C4_CONTIGS=500
        elif [[ "$SPE" == "YERS" ]]; then
            VAR_C1_GENOME=4650000
            VAR_C2_SPE="Yersinia enterocolitica"
            VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)
            VAR_C3_SNV=7
            VAR_C4_GENOME_MIN=4.1
            VAR_C4_GENOME_MAX=5.1
            VAR_C4_CONTIGS=500
            VAR_C4
        elif [[ "$SPE" == "CAMP" ]]; then
            VAR_C1_GENOME=1700000
            VAR_C2_SPE="Campylobacter jejuni"
            VAR_C3_GENOME=$(echo "scale=3; $VAR_C1_GENOME/1000000" | bc -l)
            VAR_C3_SNV=3
            VAR_C4_GENOME_MIN=1.5
            VAR_C4_GENOME_MAX=2
            VAR_C4_CONTIGS=150
        fi

        ##########################
        # DIRECTORIOS ASOCIADOS A MUESTRA
        mkdir -p ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq
        mkdir -p ${OUTPUT_PATH}/tmp/pipeline/${sample}/assembly

        ##########################################################
        #
        # 02  Assessment of the genomic sequence quality
        #
        ##########################################################
        source /software/miniconda3/etc/profile.d/conda.sh
        conda activate ${CONDAPATH}/dgsp_efsa_sp

        ################################################################################
        # FastQC Quality RAW
        ################################################################################
        fastqc --quiet --threads ${THREADS} -outdir ${OUTPUT_PATH}/1_qc/fastqc_raw "${IR}/${fastq_1}" "${IR}/${fastq_2}"

        ################################################################################
        # FastP Quality  trimming
        ################################################################################
        fastp --thread ${THREADS} --in1 "${IR}/${fastq_1}" --in2 "${IR}/${fastq_2}" \
            --cut_tail --cut_window_size=10 --cut_mean_quality=20 --length_required=50 --correction \
            --json ${OUTPUT_PATH}/1_qc/fastp/${sample}.report.fastp.json --html ${OUTPUT_PATH}/1_qc/fastp/${sample}.report.fastp.html \
            --out1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz --out2 ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz >> ${OUTPUT_PATH}/log/fastp/${sample}.log 2>&1

        # Creamos enlaces simbólicos para poder correr el pipeline
        ln -sf ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq/${sample}_1.fastq.gz
        ln -sf ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq/${sample}_2.fastq.gz

        ################################################################################
        # CHECK Q30 EFSA
        ################################################################################
        # si la longitud media de lectura = (2 × 301 pb) y Q30 < 70 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 251 pb) y Q30 < 75 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 151 pb) y Q30 < 80 %, pipeline fallará debido a la baja tasa de Q30;
        # si la longitud media de lectura = (2 × 76 pb) y Q30 < 85 %, pipeline fallará debido a la baja tasa de Q30;
        # si el número total de bases es menor que el tamaño del genoma del patógeno multiplicado por el umbral mínimo de cobertura (30x), pipeline fallará debido a una cobertura insuficiente después del recorte.

        # Recuperamos valor Q30

        BASESQ30=$(grep "Q30" ${OUTPUT_PATH}/log/fastp/${sample}.log | tail -n 2 | awk -F " " '{print $3}' | awk -F "(" '{print $1}' | awk '{total += $1}END{ printf total}')
        COVQ30=$(echo 2k $BASESQ30 $VAR_C1_GENOME /p | dc)
        VALUEQ30=$(grep "Q30" ${OUTPUT_PATH}/log/fastp/${sample}.log | tail -n 2 | awk -F " " '{print $3}' | awk -F "(" '{print $2}' | sed 's/%//g' | awk '{total += $1; count += 1}END{ printf "%4.3f\n",  total / count}')
        if ((VAR_C1_LENGTH == 301)); then
            if (($(echo "$COVQ30 < 25" | bc -l))) || (($(echo "$VALUEQ30 < 70" | bc -l))); then
                CONTROL_1_Q30="FAIL"
            else
                CONTROL_1_Q30="PASS"
            fi
        elif ((VAR_C1_LENGTH == 251)); then
            if (($(echo "$COVQ30 < 30" | bc -l))) || (($(echo "$VALUEQ30 < 75" | bc -l))); then
                CONTROL_1_Q30="FAIL"
            else
                CONTROL_1_Q30="PASS"
            fi
        elif ((VAR_C1_LENGTH == 151)); then
            if (($(echo "$COVQ30 < 30" | bc -l))) || (($(echo "$VALUEQ30 < 80" | bc -l))); then
                CONTROL_1_Q30="FAIL"
            else
                CONTROL_1_Q30="PASS"
            fi
        elif ((VAR_C1_LENGTH == 76)); then
            if (($(echo "$COVQ30 < 30" | bc -l))) || (($(echo "$VALUEQ30 < 85" | bc -l))); then
                CONTROL_1_Q30="FAIL"
            else
                CONTROL_1_Q30="PASS"
            fi
        else
            echo "Error, read lenght!"
        fi

        echo "****************************************" | tee ${OUTPUT_PATH}/log/efsa/${sample}.log
        echo "Q30 CHECK" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        echo "************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'READ_SIZE: %s\n' "$VAR_C1_LENGTH" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'GENOME_SIZE: %s\n' "$VAR_C1_GENOME" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'BASES: %s\n' "$BASESQ30" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'COV_ESTIMATED: %s\n' "$COVQ30" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'RATE Q30: %s\n' "$VALUEQ30" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
        printf 'STATUS: %s\n' "$CONTROL_1_Q30" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

        if [ "$CONTROL_1_Q30" = "PASS" ]; then
            ################################################################################
            # FastQC Quality TRIM
            ################################################################################
            fastqc --quiet --threads ${THREADS} --outdir ${OUTPUT_PATH}/1_qc/fastqc_trim ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz

            ##########################################################
            #
            # 03 Species validation
            #
            ##########################################################
            #conda activate ${CONDAPATH}/dgsp_efsa_sp

            #mash screen ${PATHMASH} ${OUTPUT_PATH}/out/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/out/0_fastq/${sample}_2.fastq.gz > ${OUTPUT_PATH}/log/mash/${sample}_88.txt

            if [[ ! -d "${OUTPUT_PATH}/tmp/pipeline/${sample}" ]]; then
                echo "Error: El directorio ${OUTPUT_PATH}/tmp/pipeline/${sample} no existe."
                exit 1
            fi

            cd ${OUTPUT_PATH}/tmp/pipeline/${sample}


            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            echo "SPECIES CHECK" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            echo "************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            printf 'SPE_LOOK: %s\n' "$VAR_C2_SPE" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

            
            # Intento 1
            nextflow-20.10.0-all run ${RESOURCES}/BACTpipe-2.7.0/bactpipe.nf \
            --mashscreen_database ${PATHMASH} \
            --reads './fastq/*_{1,2}.fastq.gz' >> ${OUTPUT_PATH}/log/BACTpipe/${sample}.log 2>&1

            # Verificar si el directorio existe después del intento 1
            if [ ! -d ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ ]; then
            	sleep 5
                # Intento 2
                nextflow-20.10.0-all run ${RESOURCES}/BACTpipe-2.7.0/bactpipe.nf \
                --mashscreen_database ${PATHMASH} \
                --reads './fastq/*_{1,2}.fastq.gz' >> ${OUTPUT_PATH}/log/BACTpipe/${sample}.log 2>&1
                
                # Verificar si el directorio existe después del intento 2
                if [ ! -d ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ ]; then
                    # Intento 2
                    sleep 10
                    nextflow-20.10.0-all run ${RESOURCES}/BACTpipe-2.7.0/bactpipe.nf \
                    --mashscreen_database ${PATHMASH} \
                    --reads './fastq/*_{1,2}.fastq.gz' >> ${OUTPUT_PATH}/log/BACTpipe/${sample}.log 2>&1

                fi
            fi

            if [ $? -eq 0 ]; then
                cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv >> ${OUTPUT_PATH}/log/efsa/${sample}.log
                C2_STATUS=$(cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk '{print $2}')
                C2_SPE=$(cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk -F "\t" '{print $5}' | sed -e "s/\[\|\]\|'//g")
                cp ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa ${OUTPUT_PATH}/2_assembly
                printf 'SPE_FIND: %s\n' "$C2_SPE" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            else
                echo "Error: El comando Nextflow falló."
                C2_SPE="None"
                printf 'SPE_FIND: %s\n' "$C2_SPE" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            fi


#            nextflow-20.10.0-all run ${RESOURCES}/BACTpipe-2.7.0/bactpipe.nf \
#            --mashscreen_database ${PATHMASH} \
#            --reads './fastq/*_{1,2}.fastq.gz' >> ${OUTPUT_PATH}/log/BACTpipe/${sample}.log 2>&1 || {
#                echo "Error: El comando Nextflow falló."
#                exit 1
#            }
#
#            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
#            echo "SPECIES CHECK" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
#            echo "************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
#            cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv >> ${OUTPUT_PATH}/log/efsa/${sample}.log
#            C2_STATUS=$(cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk '{print $2}')
#            C2_SPE=$(cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/mash_screen/all_samples.mash_screening_results.tsv | awk -F "\t" '{print $5}' | sed -e "s/\[\|\]\|'//g")
#
#            cp ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa ${OUTPUT_PATH}/2_assembly
#
#            printf 'SPE_LOOK: %s\n' "$VAR_C2_SPE" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
#            printf 'SPE_FIND: %s\n' "$C2_SPE" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

            if [[ "$C2_STATUS" = "PASS" ]] && [[ "$C2_SPE" == "$VAR_C2_SPE" ]]; then
                CONTROL_2_BACT="PASS"
                printf 'STATUS: %s\n' "$CONTROL_2_BACT" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                #23_SALM_04474	23_SALM_04474_S70_R1_001.fastq.gz	23_SALM_04474_S70_R2_001.fastq.gz
                #FILE_SEQSERO=/home/susana/DGSP/analysis_efsa/230519_LSPV008/3_typing/PLUS/23_SALM_04474.seqsero2.txt
                # CHECK SEROVARES !!!!!
                if [[ "$VAR_C2_SPE" = "Salmonella enterica" ]]; then
                    mkdir -p ${OUTPUT_PATH}/3_typing/PLUS
                    # Raw reads k-mer ("-m k"), for separated paired-end raw reads ("-t 2")
                    SeqSero2_package.py -m k -t 2 -i ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz >> ${OUTPUT_PATH}/3_typing/PLUS/${sample}.seqsero2.txt 2>&1
                    FILE_SEQSERO=${OUTPUT_PATH}/3_typing/PLUS/${sample}.seqsero2.txt
                    if [ -f "$FILE_SEQSERO" ]; then
                        STANTI=$(grep 'Predicted antigenic profile:' $FILE_SEQSERO | awk -F"\t" '{print $2}' | sed 's/,/;/g' | tr -d '[:space:]')
                        STPATHO=$(grep 'Predicted serotype:' $FILE_SEQSERO | awk -F: '{print $2}' | sed 's/,/;/g' | tr -d '[:space:]')
                    else
                        STANTI="-"
                        STPATHO="-"
                    fi
                fi

                if [[ "$VAR_C2_SPE" = "Listeria monocytogenes" ]]; then
                    mkdir -p ${OUTPUT_PATH}/3_typing/PLUS
                    #lissero ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa >> ${OUTPUT_PATH}/3_typing/${sample}.lissero.txt 2>&1
                    var2erase=${OUTPUT_PATH}/2_assembly/
                    #IMPORTANTE EN SED SE EMPLEA ":" PARA SEPARAS CAMPOS Y PODER ELIMINAR PATH PARCIA
                    lissero ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa --logfile ${OUTPUT_PATH}/log/lissero/${sample}.log | sed "s:^$var2erase::" > ${OUTPUT_PATH}/3_typing/PLUS/${sample}.lissero.txt

                    FILE_LISSERO=${OUTPUT_PATH}/3_typing/PLUS/${sample}.lissero.txt
                    if [ -f "$FILE_LISSERO" ]; then
                        STANTI=$(awk -F"\t" '{print $2}' ${FILE_LISSERO}| tail -1 | sed 's/, /;/g')
                    else
                        STANTI="-"
                        STPATHO="-"
                    fi
               fi


                ##########################################################
                #
                # 04 Contamination
                #
                ##########################################################
                echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                echo "CONTAMINATION CHECK" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                echo "********************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                # Updated database rMLST for Campylobacter jejuni
                # https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases
                # registro en https://pubmlst.org/bigsdb
                # Campylobacter jejuni/coli isolates (pubmlst_campylobacter_isolates)
                # Campylobacter jejuni/coli typing (pubmlst_campylobacter_seqdef)
                # Ribosomal MLST genomes (pubmlst_rmlst_isolates)
                # Ribosomal MLST typing (pubmlst_rmlst_seqdef)
                ## Descargo los MLST de varios genes: aspA.fas  glnA.fas  gltA.fas  glyA.fas  pgm.fas  tkt.fas  uncA.fas
                # cat *.fas > Campylobacter_db_cgderived.fasta
                ## Ve voy a la base de datos que está en /home/asanzc/.confindr_db y copio este fichero
                ## Hago copia del fichero refseq
                # mv refseq.msh saved_refseq.msh
                ## combino los ficheros anteriores y genero un nuevo refseq.msh incluyendo campylobacter
                # cat *fasta > combined.fasta
                #
                # mash sketch -o refseq combined.fasta
                ## elimino combined.fasta
                ## Esto no ha funcionado, por que no se ha añadido bien
                ## dejo el archivo refseq.msh antiguo
                # mv saved_refseq.msh refseq.msh
                # Creo fichero download_loci.sh para descargar cgMLST de Campylobacter
                #chmod +x download_loci.sh
                #./download_loci.sh ~/.confindr_db/database_Campylobacter_cgMLST
                # Eliminamos comentarios como CAMP0123 NO ES UN LOCUS! DEL FICHERO FASTA FINAL
                #cat ~/.confindr_db/database_Campylobacter_cgMLST/*fasta | grep -v "^CAMP" > ~/.confindr_db/Campylobacter_cgMLST.fasta
                # De momento se crea excpeción para que se coga el mlst

                if [[ "$SPE" == "CAMP" ]]; then
                    # Para que sólo coja los fastq de la muestra en concreto
                    #confindr.py -t $THREADS -i ${OUTPUT_PATH}/0_fastq/${sample}/fastq -o ConFindr --cgmlst ~/.confindr_db/Campylobacter_cgMLST.fasta >> ${OUTPUT_PATH}/log/confindr/${sample}.log 2>&1
                    confindr.py -t $THREADS -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq -o ConFindr --cgmlst ${HOME}/.confindr_db/Campylobacter_jejuni-coli.fasta >> ${OUTPUT_PATH}/log/confindr/${sample}.log 2>&1
                else
                    confindr.py -t $THREADS -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq -o ConFindr >> ${OUTPUT_PATH}/log/confindr/${sample}.log 2>&1
                fi

                conda deactivate

                conda activate ${CONDAPATH}/dgsp_efsa_contamination

                INNUca.py -j ${THREADS} \
                    -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/fastq \
                    -o INNUca \
                    -s "$VAR_C2_SPE" -g $VAR_C3_GENOME >> ${OUTPUT_PATH}/log/innuca/${sample}.log 2>&1

                # Recogemos resultados contaminación
                C3_confindr=$(awk -F "," '{print $4}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/ConFindr/confindr_report.csv | grep -v "ContamStatus" | uniq)
                C3_warning=$(awk -F "," '{print $3}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/ConFindr/confindr_report.csv | grep -v "NumContamSNVs" | awk '{total += $1; count += 1}END{ printf "%4.3f\n",  total / count}')

                conda deactivate

                #Si el número de SNV contaminados supera los siguientes umbrales específicos de especies
                #(como lo sugiere Deneke et al., 2021a), se envía un mensaje de advertencia:
                #• S. enterica: 7
                #• L. monocytogenes: 3
                #• E.coli: 4
                if (($(echo "$C3_warning > $VAR_C3_SNV" | bc -l))); then
                    echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                    echo "========================================" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                    echo "WARNING!" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                    echo "ConFindr show SNV contamination" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                    echo "========================================" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                else
                    repinnuca=$(ls -t ${OUTPUT_PATH}/tmp/pipeline/${sample}/INNUca/samples_report.* | head -n 1)
                    C3_innuca=$(awk -F "\t" '{print $15}' $repinnuca | tail -n 1)

                    if [[ "$C3_confindr" = "False" ]] && { [[ "$C3_innuca" = "PASS" ]] || [[ "$C3_innuca" = "WARNING" ]]; }; then

                        CONTROL_3_CONT="PASS"
                        printf 'STATUS: %s\n' "$CONTROL_3_CONT" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                        ##########################################################
                        #
                        # 05 Assembly
                        #
                        ##########################################################
                        conda activate ${CONDAPATH}/dgsp_efsa_sp

                        #shovill --depth 100 --kmers 31,33,55,77,99,127 --minlen 500 --cpus ${THREADS} \
                        #--R1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                        #--R2 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                        #--outdir ${OUTPUT_PATH}"/pipeline/"${sample}"/assembly/shovill"

                        #cp ${OUTPUT_PATH}/pipeline/${sample}/assembly/shovill/contigs.fa ${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.contigs.fa
                        #statswrapper.sh in=${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.contigs.fa > ${OUTPUT_PATH}/pipeline/${sample}/assembly/${sample}.assembly_stats.txt

                        echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        echo "Assembly quality CHECKM" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        echo "************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        #mkdir -p ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs
                        #Desglosamos fichero contigs en fastas individuales
                        #(cd ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs && cat ../${sample}.contigs.fa | awk '{
                        #(
                        #    cd ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs && cat ../../BACTpipe_results/shovill/${sample}.contigs.fa | awk '{
                        #    if (substr($0, 1, 1)==">") {filename=(substr($1,2) ".fa")}
                        #    print $0 >> filename
                        #    close(filename)
                        #    }'
                        #)

                        #checkm lineage_wf -t ${THREADS} \
                        #    -x fa ${OUTPUT_PATH}/pipeline/${sample}/assembly/contigs \
                        #    ${OUTPUT_PATH}/pipeline/${sample}/checkm >> ${OUTPUT_PATH}/log/checkm/${sample}.log 2>&1

                        #checkm lineage_wf -t ${THREADS} \
                        #    -x fa ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/\
                        #    ${OUTPUT_PATH}/pipeline/${sample}/checkm >> ${OUTPUT_PATH}/log/checkm/${sample}.log 2>&1

                        # En caso de no tener más de 16GB opción --reduced_tree
                        #checkm lineage_wf --reduced_tree -t ${THREADS} \
                        #    -x fa ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/\
                        #    ${OUTPUT_PATH}/pipeline/${sample}/checkm >> ${OUTPUT_PATH}/log/checkm/${sample}.log 2>&1

                        #checkm tree --reduced_tree --ali -t ${THREADS} -x fa \
                        #    ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/ \
                        #    ${OUTPUT_PATH}/pipeline/${sample}/checkm \
                        #    >> ${OUTPUT_PATH}/log/checkm/${sample}_01_tree.log 2>&1

                        #checkm tree_qa -o 2 --tab_table \
                        #-f ${OUTPUT_PATH}/pipeline/${sample}/checkm/tree_qa.tsv \
                        #${OUTPUT_PATH}/pipeline/${sample}/checkm \
                        #>> ${OUTPUT_PATH}/log/checkm/${sample}_02_tree_qa.log 2>&1

                        #http://www.mselab.cn/detail/81/
                        #checkm data setRoot /software/resources/checkm_data
                        # https://github.com/Ecogenomics/CheckM/issues/244
                        # https://github.com/Ecogenomics/CheckM/wiki/Workflows#using-custom-marker-genes

                        # --reduced_tree \
                        mkdir -p ${OUTPUT_PATH}/tmp/pipeline/${sample}/taxonomy
                        checkm lineage_wf \
                        -t ${THREADS} \
                        -f ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/lineage_wf.tsv \
                        --nt \
                        --tab_table \
                        --pplacer_threads ${THREADS} \
                        -x fa \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_01_lineage_wf.log 2>&1

                        checkm taxonomy_wf -t ${THREADS} \
                        -x fa \
                        --tab_table \
                        -f ${OUTPUT_PATH}/tmp/pipeline/${sample}/taxonomy/checkm_results \
                        domain Bacteria \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/taxonomy \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_02_taxonomy.log 2>&1

                        checkm tree_qa -o 2 --tab_table \
                        -f ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/tree_qa.tsv \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_03_tree_qa.log 2>&1

                        checkm tetra ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/tetra_profile.tsv \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_04_tetra_qa.log 2>&1


                        #https://136.159.60.117/metaerg/mockEvenCell/metaspades/metabat1500/SCG/
                        #checkm lineage_set <output folder> <marker file>
                        #checkm analyze <marker file> <bin folder> <output folder>
                        #checkm qa <marker file> <output folder>

                        #mkdir /ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/checkm
                        #cd /ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/checkm
                        #wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz
                        #wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz
                        #wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/metadata_field_desc.tsv
                        #tar -xzvf ar122_metadata_r95.tar.gz
                        #tar -xzvf bac120_metadata_r95.tar.gz
                        #cat ${OUTPUT_PATH}/pipeline/${sample}/checkm/storage/bin_stats.analyze.tsv | \
                        #awk -F "{" '{print $2'} | sed "s/}//g" | sed "s/, /\n/g" | sed "s/: /\t/" | \
                        #sed -e "s/\# //g" | sed "s/'//g" | sed "s/ /_/g" > ${OUTPUT_PATH}/pipeline/${sample}/checkm/${sample}_checkm.tsv
                        # nf-core modules install checkm/qa

                        # https://www.biostars.org/p/447744/
                        # ChackM uses single-copy genes to evaluate the completeness and contamination of a genome (or a pseudo-genome).
                        # If all the genes are found in the genome then completeness is 100% (since they are all essential proteins).
                        # If they appear more than once then it's probably contaminated (because two copies are usually lethal).
                        # So for P14_4 you can see that there are 104 markers, 39 of which appear only once, 32 appear twice and 31 three
                        # times (2 appear 4 times) so since all the genes are found the genome is probably complete, but since there are
                        # multiple copies you are probably looking at 2.3 genomes instead of one (that's 130% contamination), 91.24% of
                        # the contamination is probably from another strain of the main bacteria.

                        #----------------------------------------------------------------------------------------------------------------------------------------------------------
                        #  genes    c__Bacilli (UID354)      515         328           182        1   314   12   1   0   0       99.45            5.22              86.67
                        #----------------------------------------------------------------------------------------------------------------------------------------------------------
                        # Juntamos columnas que nos interesa y reformateamos
                        # Formateo resultado analisis
                        awk -F "\t" '{print $2"\t"$12"\t"$13"\t"$14}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/lineage_wf.tsv | \
                        awk -F "\t" '{
                                        for (f = 1; f <= NF; f++) { a[NR, f] = $f } 
                                    }
                                    NF > nf { nf = NF }
                                    END {
                                        for (f = 1; f <= nf; f++) {
                                            for (r = 1; r <= NR; r++) {
                                                printf a[r, f] (r==NR ? RS : FS)
                                            }
                                        }
                                    }' | \
                        sed "s/ /_/g" > ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv

                        # Métricas tree
                        awk -F "\t" '{print $5"\t"$8"\t"$9"\t"$15"\t"$16"\t"$17"\t"$18}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/tree_qa.tsv | \
                        awk -F "\t" '{
                                        for (f = 1; f <= NF; f++) { a[NR, f] = $f } 
                                    }
                                    NF > nf { nf = NF }
                                    END {
                                        for (f = 1; f <= nf; f++) {
                                            for (r = 1; r <= NR; r++) {
                                                printf a[r, f] (r==NR ? RS : FS)
                                            }
                                        }
                                    }' | \
                        sed "s/ /_/g" | sed "s/Lineage:_/out_/g" >> ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv

                        #Métricas ensamblaje
                        cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/storage/bin_stats.analyze.tsv | \
                        awk -F "{" '{print $2}' | sed "s/}//g" | sed "s/, /\n/g" | sed "s/: /\t/" | \
                        sed -e "s/\# //g" | sed "s/'//g" | sed "s/ /_/g" | grep "contigs" \
                        >> ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv

                        # Gráficas
                        mkdir -p ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/plots

                        checkm marker_plot ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/ \
                        -x fa ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/plots \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_05_marker_plot.log 2>&1

                        checkm len_hist \
                        -x fa ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/plots \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_06_len_hist.log 2>&1

                        checkm nx_plot \
                        -x fa ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/plots \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_07_nx_plot.log 2>&1

                        checkm tetra_plot \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/ \
                        -x fa ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/plots \
                        ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/tetra_profile.tsv \
                        95 \
                        >> ${OUTPUT_PATH}/log/checkm/${sample}_08_tetra_plot.log 2>&1


                        # ponemos resultados también en documento efsa
                        cat ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv >> ${OUTPUT_PATH}/log/efsa/${sample}.log
                        echo "________________________________________" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                        C4_completeness=$(grep "Completeness" ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | awk '{print $2}')
                        C4_contamination=$(grep "Contamination" ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | awk '{print $2}')

                        # Extraer la segunda columna y unirla con ';'
                        #vcheckm=$(awk '{print $2}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | tr '\n' ',')
                        vcheckm=$(awk '{print $2}' ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | tr '\n' ',' | awk -F "," '{print "\""$1"\"",$2,$3,$4,"\"" $5 "\"",$6,$7,$8,$9,$10,$11,$12,$13}' | tr '\ ' ',')


                        # Eliminar el último ','
                        #vcheckm=${vcheckm%?}


                        if (($(echo "$C4_completeness < $VAR_C4_COMPLETENESS" | bc -l))) || (($(echo "$C4_contamination > $VAR_C4_CONTAMINATION" | bc -l))); then
                            CONTROL4_CHECKM="FAIL"
                            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "========================================" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "WARNING!" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "CHECKM control has not passed" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "========================================" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        else
                            CONTROL4_CHECKM="PASS"
                            printf 'STATUS: %s\n' "$CONTROL4_CHECKM" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                            C4_genome_size=$(grep "out_genome_size_(Mbp)_mean" ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | awk '{print $2}')
                            C4_contigs=$(grep "contigs" ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | grep -v "N50" | awk '{print $2}')
                            C4_N50=$(grep "contigs" ${OUTPUT_PATH}/tmp/pipeline/${sample}/checkm/${sample}_checkm.tsv | grep "N50" | awk '{print $2}')


                            # Vemos si el genoma está dentro del rango
                            if (($(echo "$C4_genome_size >= $VAR_C4_GENOME_MIN"| bc -l))) && (($(echo "$C4_genome_size <= $VAR_C4_GENOME_MAX" | bc -l))); then
                                C4_genome_size_control="PASS"
                            else
                                C4_genome_size_control="FAIL"
                            fi


                            if (($(echo "$C4_N50 < 30000" | bc -l))) || (($(echo "$C4_contigs > $VAR_C4_CONTIGS" | bc -l))) || [[ "$C3_confindr" = "True" ]] || [[ "$C4_genome_size_control" = "FAIL" ]]; then
                                CONTROL4_CHECKM_quality="ASSEMBLY_QUALITY: BAD"
                                vCONTROL4_CHECKM_quality="BAD"

                            else
                                CONTROL4_CHECKM_quality="ASSEMBLY_QUALITY: GOOD"
                                vCONTROL4_CHECKM_quality="GOOD"
                            fi

                            echo ${CONTROL4_CHECKM_quality} | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log



                            ##########################################################
                            #
                            # 05 typing
                            #
                            ##########################################################
                            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "TYPING" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log

                            ##########################################################

                            #mkdir -p ${OUTPUT_PATH}"/typing/"${sample}

                            # Copiamos contigs en carpeta typing
                            #cp ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa ${OUTPUT_PATH}"/typing/"${sample}

                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            if [[ "$VAR_C2_SPE" = "Salmonella enterica" ]]; then
                                # Raw reads k-mer ("-m k"), for separated paired-end raw reads ("-t 2")
                                #SeqSero2_package.py -m k -t 2 -i ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz >> ${OUTPUT_PATH}/pipeline/${sample}/${sample}.seqsero2.txt 2>&1
                                #cp ${OUTPUT_PATH}/pipeline/${sample}/${sample}.seqsero2.txt ${OUTPUT_PATH}/pipeline/${sample}/SEROtyping_${sample}_seqsero2.txt
                                #---------------------------------------------------------
                                # 05.2 Detección de AMR
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/AMR
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} genomicepidemiology/resfinder:4.3.3 \
                                -ifa ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                -o ${OUTPUT_PATH}/3_typing/AMR/${sample} -s senterica --acquired --point \
                                >> ${OUTPUT_PATH}/log/resfinder/${sample}_resfinder.log 2>&1


                                FILE_RESISTANT="${OUTPUT_PATH}/3_typing/AMR/${sample}/pheno_table.txt"
                                if [[ -f "$FILE_RESISTANT" ]]; then
                                    RESFINDER_RESISTANT=$(awk -F'\t' '{gsub(/ /, "_", $1); gsub(/ /, "_", $2)} $3 == "Resistant" && $4 == 3 {print $1"("$2")"}' "$FILE_RESISTANT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_RESISTANT="-"
                                fi

                                FILE_RESGENE="${OUTPUT_PATH}/3_typing/AMR/${sample}/ResFinder_results_tab.txt"
                                if [[ -f "$FILE_RESGENE" ]]; then
                                    RESFINDER_GENES=$(awk -F'\t' '$4 ~ /^[0-9]+(\.[0-9]+)?$/ && $4 >= 90 {gsub(/ /, "_", $1); print $1}' "$FILE_RESGENE" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_GENES="-"
                                fi

                                FILE_RESMUT="${OUTPUT_PATH}/3_typing/AMR/${sample}/PointFinder_results.txt"
                                if [[ -f "$FILE_RESMUT" ]]; then
                                    RESFINDER_MUTS=$(awk -F'\t' 'NR>1 {gsub(/ /, "_", $1); gsub(/ /, "_", $4); gsub(/,_/, ";", $4); print $1"("$4")"}' "$FILE_RESMUT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_MUTS="-"
                                fi

                                [ -z "$RESFINDER_GENES" ] && RESFINDER_GENES="-"
                                [ -z "$RESFINDER_RESISTANT" ] && RESFINDER_RESISTANT="-"
                                [ -z "$RESFINDER_MUTS" ] && RESFINDER_MUTS="-"
                                # docker run -it staphb/ncbi-amrfinderplus /bin/bash
                                # amrfinder -p ecoli.faa --plus -o AMRFinder_complete.tsv --threads 4 --ident_min $(echo "scale=2; 90/100" | bc -l ) \
                                # --coverage_min $(echo "scale=2; 80/100" | bc -l ) --name ecoli --protein_output ecoli_args.faa --database /amrfinder/data/2023-04-17.1 
                                # awk -F '	' '{ if ($3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;

                                #---------------------------------------------------------
                                # 05.3 MLST y perfil de virulencia
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/MLST

                                # stringMLST.py \
                                # --predict \
                                # -1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                                # -2 ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                # -p \
                                # --prefix ${RESOURCES}/string_mlst_db/Escherichia_coli_1/Escherichia_coli_1 \
                                # -k 35 \
                                # -o ${OUTPUT_PATH}/3_typing/MLST/${sample}_stringmlst_ecoli1.txt

                                # stringMLST.py \
                                # --predict \
                                # -1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                                # -2 ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                # -p \
                                # --prefix ${RESOURCES}/string_mlst_db/Escherichia_coli_2/Escherichia_coli_2 \
                                # -k 35 \
                                # -o ${OUTPUT_PATH}/3_typing/MLST/${sample}_stringmlst_ecoli2.txt

                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/sistr
                                sistr \
                                --qc \
                                -vv \
                                --alleles-output ${OUTPUT_PATH}/3_typing/PLUS/sistr/${sample}_allele-results.json \
                                --novel-alleles ${OUTPUT_PATH}/3_typing/PLUS/sistr/${sample}_novel-alleles.fasta \
                                --cgmlst-profiles ${OUTPUT_PATH}/3_typing/PLUS/sistr/${sample}_cgmlst-profiles.csv \
                                -f tab \
                                -o ${OUTPUT_PATH}/3_typing/PLUS/sistr/${sample}_sistr-output.tab \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                > ${OUTPUT_PATH}/log/sistr/${sample}.log 2>&1

                                FILE_SISTR=${OUTPUT_PATH}/3_typing/PLUS/sistr/${sample}_sistr-output.tab

                                if [[ -e $FILE_SISTR ]]; then
                                    ST_CGMLST=$(awk -F'\t' '{print $1}' $FILE_SISTR | tail -1 | sed 's/,/;/g')
                                    ST_H1=$(awk -F'\t' '{print $9}' $FILE_SISTR | tail -1 | sed 's/,/;/g')
                                    ST_H2=$(awk -F'\t' '{print $10}' $FILE_SISTR | tail -1 | sed 's/,/;/g')
                                    ST_O_antigen=$(awk -F'\t' '{print $11}' $FILE_SISTR | tail -1 | sed 's/,/;/g')
                                    SEROVAR=$(awk -F'\t' '{print $15}' $FILE_SISTR | tail -1 | sed 's/,/;/g')
                                else
                                    ST_CGMLST="-"
                                    ST_H1="-"
                                    ST_H2="-"
                                    ST_O_antigen="-"
                                    SEROVAR="-"
                                fi

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/mlst:2.23.0 \
                                mlst \
                                --threads ${THREADS} \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out \
                                2> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.log

                                # Limpiamos un poco el fichero
                                var2erase=${OUTPUT_PATH}/2_assembly/
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out 

                                FILE_MLST=${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out
                                if [[ -f "$FILE_MLST" ]]; then
				    MLST=$(awk -F'\t' '{result=""; for(i=4; i<=NF; i++) result = result $i ";"; sub(/;$/, "", result); gsub(/,/, "-", result); print result}' "$FILE_MLST")
				    STVAR=$(awk '{if(NR==1) print $3}' "$FILE_MLST")
                                else
                                    MLST="-"
                                    STVAR="-"
                                fi

                                # ABRICATE
                                ##########
                                # Para ver bases de datos disponibles
                                #docker run --rm staphb/abricate:1.0.1-insaflu-220727 abricate --list
                                # Nos devuelve:
                                # DATABASE	SEQUENCES	DBTYPE	DATE
                                # argannot	2223	nucl	2022-Nov-16
                                # ncbi	5386	nucl	2022-Nov-16
                                # plasmidfinder	460	nucl	2022-Nov-16
                                # vfdb	2597	nucl	2022-Nov-16
                                # ecoli_vf	2701	nucl	2022-Nov-16
                                # card	2631	nucl	2022-Nov-16
                                # megares	6635	nucl	2022-Nov-16
                                # resfinder	3077	nucl	2022-Nov-16
                                # ecoh	597	nucl	2022-Nov-16
                                # insaflu	34	nucl	2022-Nov-16
                                mkdir -p ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                2> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.log

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                --summary ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                > ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # Limpiamos un poco el fichero
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # awk '{print $6}' ${OUTPUT_PATH}"/typing/"${sample}"/MLST_VIRULENCE-PROFILE/"${sample}"_abricate.out" | uniq | grep -v "GENE" | awk 'BEGIN { ORS = " " } { print }'
                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/ariba
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHARIBA}:${RESOURCES}/ariba staphb/ariba:latest \
                                ariba run \
                                --threads ${THREADS} --force \
                                ${PATHARIBA}/Salmonella_enterica/ref_db \
                                ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                ${OUTPUT_PATH}/3_typing/PLUS/ariba/${sample} \
			                    > ${OUTPUT_PATH}/log/ariba/${sample}.log 2>&1

                                #---------------------------------------------------------
                                # 05.4  Allele calling
                                #---------------------------------------------------------
                                #mkdir -p ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}
                                #ln -sf ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}
                                #PATHMLST=/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/cgMLST_data

                                chewBBACA_BLAST_SCORE_RATIO=0.6
                                chewBBACA_MINIMUM_LENGTH=0
                                chewBBACA_TRANSLATION_TABLE=11
                                chewBBACA_SIZE_THRESHOLD=None
                                BACT=Salmonella_enterica

                                # Ya venimos con el esquema descargado y adaptado según las especificaciones de la EFSA
                                # docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/chewbbaca:3.1.2 \
                                # chewBBACA.py CreateSchema \
                                # --cpu ${THREADS} \
                                # -i ${PATHMLST}"/"${BACT}"/ecoli_INNUENDO_wgMLST/" \
                                # --bsr ${chewBBACA_BLAST_SCORE_RATIO} \
                                # --l ${chewBBACA_MINIMUM_LENGTH} \
                                # --t ${chewBBACA_TRANSLATION_TABLE} \
                                # --st ${chewBBACA_SIZE_THRESHOLD} \
                                # -o ${OUTPUT_PATH}"/typing/"${sample}"/ALLELE_CALLING/scheme" \
                                # --ptf ${PATHMLST}"/"${BACT}"/ecoli_INNUENDO_wgMLST/"${BACT}".trn"

                                # Allele Calling
                                #################
                                # The allele call step determines the allelic profiles of the analyzed strains,
                                # identifying known and novel alleles in the analyzed genomes.
                                # Novel alleles are assigned an allele identifier and added to the schema.
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py AlleleCall \
                                --cpu ${THREADS} --force-continue \
                                -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ \
                                --bsr ${chewBBACA_BLAST_SCORE_RATIO} \
                                --l ${chewBBACA_MINIMUM_LENGTH} \
                                --t ${chewBBACA_TRANSLATION_TABLE} \
                                --st ${chewBBACA_SIZE_THRESHOLD} \
                                -g ${PATHcgMLST}/schema/${BACT} \
                                --ptf ${PATHcgMLST}/schema/${BACT}/${BACT}.trn \
                                -o ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_1_allelecall.log 2>&1

                                # Get PATH results
                                PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}/$(ls ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample})
                                #PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}
                                # Paralog detection
                                ###################
                                # The next step in the analysis is to determine if some of the loci can be considered paralogs
                                # based on the result of the wgMLST allele calling. The AlleleCall module returns a list of
                                # Paralogous genes in the paralogous_counts.tsv file that can be found on the results32_wgMLST folder.
                                # The paralogous_counts.tsv file contains a set of 10 loci that were identified as possible paralogs.
                                # These loci should be removed from the schema due to the potential uncertainty in allele assignment.
                                # To remove the set of 10 paralogous loci from the allele calling results, run the following command:
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py RemoveGenes -i ${PATHALLELE}/results_alleles.tsv -g ${PATHALLELE}/paralogous_counts.tsv -o ${PATHALLELE}/results_alleles_NoParalogs.tsv \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_2_removegenes.log 2>&1

                                # Evaluate genome quality NOT AVAILABLE IN THIS chewBBACA VERSION!
                                # docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                # chewBBACA.py TestGenomeQuality -i ${PATHALLELE}/results_alleles_NoParalogs.tsv -n 13 -t 300 -s 5 -o ${PATHALLELE}/genome_quality


                                # ExtractCgMLST - Determine the set of loci that constitute the core genome
                                ###################
                                # We can now determine the set of loci in the core genome based on the allele calling results.
                                # The set of loci in the core genome is determined based on a threshold of loci presence in the analysed genomes.
                                # We can run the ExtractCgMLST module to determine the set of loci in the core genome for the loci presence thresholds of 95%, 99% and 100%.
                                ###
                                # Requirements to define a core genome MLST (cgMLST) schema
                                # cgMLST schemas are defined as the set of loci that are present in all strains under analysis or, due to
                                # sequencing/assembly limitations, >95% of strains analyzed. In order to have a robust definition of a cgMLST
                                # schema for a given bacterial species, a set of representative strains of the diversity of a given species
                                # should be selected. Furthermore, since cgMLST schema definition is based on pre-defined thresholds, only
                                # when a sufficient number of strains have been analyzed can the cgMLST schema be considered stable.
                                # This number will always depend on the population structure and diversity of the species in question, with
                                # non-recombinant monomorphic species possibly requiring a smaller number of strais to define cgMLST schemas
                                # than panmictic highly recombinogenic species that are prone to have large numbers of accessory genes and
                                # mobile genetic elements. It is also important to refer that the same strategy described here can be used
                                # to defined lineage specific schemas for more detailed analysis within a given bacterial lineage. Also,
                                # by definition, all the loci that are not considered core genome, can be classified as being part of an accessory genome MLST (agMLST) schema.
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py ExtractCgMLST -i ${PATHALLELE}/results_alleles.tsv -o ${PATHALLELE} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_3_extractgcmlst.log 2>&1

                                # rm -r ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}

                                #echo -e "$sample,$fastq_1,$fastq_2,$VAR_C2_SPE,$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,$CONTROL_1_Q30,$C2_SPE,$CONTROL_2_BACT,$C3_warning,$CONTROL_3_CONT,$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,$CONTROL4_CHECKM,$vCONTROL4_CHECKM_quality,$STVAR,$MLST,$STANTI,$STPATHO,$ST_CGMLST,$ST_H1,$ST_H2,$ST_O_antigen,$SEROVAR,$RESFINDER_GENES,$RESFINDER_RESISTANT,$RESFINDER_MUTS" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
				echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",$C3_warning,\"$CONTROL_3_CONT\",$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,\"$CONTROL4_CHECKM\",\"$vCONTROL4_CHECKM_quality\",\"$STVAR\",\"$MLST\",\"$STANTI\",\"$STPATHO\",\"$ST_CGMLST\",\"$ST_H1\",\"$ST_H2\",\"$ST_O_antigen\",\"$SEROVAR\",\"$RESFINDER_GENES\",\"$RESFINDER_RESISTANT\",\"$RESFINDER_MUTS\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"


                            elif [[ "$VAR_C2_SPE" = "Listeria monocytogenes" ]]; then

                                #---------------------------------------------------------
                                # 05.0 Declaramos variables que no recogemos en esta especie
                                #---------------------------------------------------------
                                ST_CGMLST="-"
                                ST_H1="-"
                                ST_H2="-"
                                ST_O_antigen="-"
                                SEROVAR="-"

                                #---------------------------------------------------------
                                # 05.1 Serotipado y patotipado
                                #---------------------------------------------------------
                                # (cd ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill && lissero ${sample}.contigs.fa >> ${OUTPUT_PATH}/pipeline/${sample}/${sample}.lissero.txt 2>&1)
                                # var2erase=${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/
                                #IMPORTANTE EN SED SE EMPLEA ":" PARA SEPARAS CAMPOS Y PODER ELIMINAR PATH PARCIAL
                                # lissero ${OUTPUT_PATH}/pipeline/${sample}/BACTpipe_results/shovill/${sample}.contigs.fa --logfile ${OUTPUT_PATH}/log/lissero/${sample}.log | sed "s:^$var2erase::" > ${OUTPUT_PATH}/pipeline/${sample}/${sample}.lissero.txt


                                #---------------------------------------------------------
                                # 05.2 Detección de AMR
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/AMR
                                # docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} genomicepidemiology/resfinder:4.3.3 \
                                # -ifa ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                # -o ${OUTPUT_PATH}/3_typing/AMR/${sample} -s listeria --acquired --point \
                                # >> ${OUTPUT_PATH}/log/resfinder/${sample}_resfinder.log 2>&1


                                FILE_RESISTANT="${OUTPUT_PATH}/3_typing/AMR/${sample}/pheno_table.txt"
                                if [[ -f "$FILE_RESISTANT" ]]; then
                                    RESFINDER_RESISTANT=$(awk -F'\t' '{gsub(/ /, "_", $1); gsub(/ /, "_", $2)} $3 == "Resistant" && $4 == 3 {print $1"("$2")"}' "$FILE_RESISTANT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_RESISTANT="-"
                                fi

                                FILE_RESGENE="${OUTPUT_PATH}/3_typing/AMR/${sample}/ResFinder_results_tab.txt"
                                if [[ -f "$FILE_RESGENE" ]]; then
                                    RESFINDER_GENES=$(awk -F'\t' '$4 ~ /^[0-9]+(\.[0-9]+)?$/ && $4 >= 90 {gsub(/ /, "_", $1); print $1}' "$FILE_RESGENE" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_GENES="-"
                                fi

                                FILE_RESMUT="${OUTPUT_PATH}/3_typing/AMR/${sample}/PointFinder_results.txt"
                                if [[ -f "$FILE_RESMUT" ]]; then
                                    RESFINDER_MUTS=$(awk -F'\t' 'NR>1 {gsub(/ /, "_", $1); gsub(/ /, "_", $4); gsub(/,_/, ";", $4); print $1"("$4")"}' "$FILE_RESMUT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_MUTS="-"
                                fi

                                [ -z "$RESFINDER_GENES" ] && RESFINDER_GENES="-"
                                [ -z "$RESFINDER_RESISTANT" ] && RESFINDER_RESISTANT="-"
                                [ -z "$RESFINDER_MUTS" ] && RESFINDER_MUTS="-"
                                # docker run -it staphb/ncbi-amrfinderplus /bin/bash
                                # amrfinder -p ecoli.faa --plus -o AMRFinder_complete.tsv --threads 4 --ident_min $(echo "scale=2; 90/100" | bc -l ) \
                                # --coverage_min $(echo "scale=2; 80/100" | bc -l ) --name ecoli --protein_output ecoli_args.faa --database /amrfinder/data/2023-04-17.1 
                                # awk -F '	' '{ if ($3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;

                                #---------------------------------------------------------
                                # 05.3 MLST y perfil de virulencia
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/MLST

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/mlst:2.23.0 \
                                mlst \
                                --threads ${THREADS} \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out \
                                2> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.log

                                # Limpiamos un poco el fichero
                                var2erase=${OUTPUT_PATH}/2_assembly/
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out 

                                FILE_MLST=${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out
                                if [[ -f "$FILE_MLST" ]]; then
			            MLST=$(awk -F'\t' '{result=""; for(i=4; i<=NF; i++) result = result $i ";"; sub(/;$/, "", result); gsub(/,/, "-", result); print result}' "$FILE_MLST")
				    STVAR=$(awk '{if(NR==1) print $3}' "$FILE_MLST")
                                else
                                    MLST="-"
                                    STVAR="-"
                                fi

                                # ABRICATE
                                ##########
                                mkdir -p ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                2> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.log

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                --summary ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                > ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # Limpiamos un poco el fichero
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # awk '{print $6}' ${OUTPUT_PATH}"/typing/"${sample}"/MLST_VIRULENCE-PROFILE/"${sample}"_abricate.out" | uniq | grep -v "GENE" | awk 'BEGIN { ORS = " " } { print }'
                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/ariba
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHARIBA}:${RESOURCES}/ariba staphb/ariba:latest \
                                ariba run \
                                --threads ${THREADS} --force \
                                ${PATHARIBA}/Listeria_monocytogenes/ref_db \
                                ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                ${OUTPUT_PATH}/3_typing/PLUS/ariba/${sample} \
			                    > ${OUTPUT_PATH}/log/ariba/${sample}.log 2>&1

                                #---------------------------------------------------------
                                # 05.4  Allele calling
                                #---------------------------------------------------------
                                #mkdir -p ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}
                                #ln -sf ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}
                                #PATHMLST=/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/cgMLST_data

                                chewBBACA_BLAST_SCORE_RATIO=0.6
                                chewBBACA_MINIMUM_LENGTH=144
                                chewBBACA_TRANSLATION_TABLE=11
                                chewBBACA_SIZE_THRESHOLD=0.2
                                BACT=Listeria_monocytogenes

                                # Allele Calling
                                #################

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py AlleleCall \
                                --cpu ${THREADS} --force-continue \
                                -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ \
                                --bsr ${chewBBACA_BLAST_SCORE_RATIO} \
                                --l ${chewBBACA_MINIMUM_LENGTH} \
                                --t ${chewBBACA_TRANSLATION_TABLE} \
                                --st ${chewBBACA_SIZE_THRESHOLD} \
                                -g ${PATHcgMLST}/schema/${BACT} \
                                --ptf ${PATHcgMLST}/schema/${BACT}/${BACT}.trn \
                                -o ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_1_allelecall.log 2>&1

                                # Get PATH results
                                PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}/$(ls ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample})
                                #PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}
                                # Paralog detection
                                ###################
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py RemoveGenes -i ${PATHALLELE}/results_alleles.tsv -g ${PATHALLELE}/paralogous_counts.tsv -o ${PATHALLELE}/results_alleles_NoParalogs.tsv \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_2_removegenes.log 2>&1

                                # ExtractCgMLST - Determine the set of loci that constitute the core genome
                                ###################
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py ExtractCgMLST -i ${PATHALLELE}/results_alleles.tsv -o ${PATHALLELE} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_3_extractgcmlst.log 2>&1

                                # rm -r ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}

                                #echo -e "$sample,$fastq_1,$fastq_2,$VAR_C2_SPE,$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,$CONTROL_1_Q30,$C2_SPE,$CONTROL_2_BACT,$C3_warning,$CONTROL_3_CONT,$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,$CONTROL4_CHECKM,$vCONTROL4_CHECKM_quality,$STVAR,$MLST,$STANTI,$STPATHO,$ST_CGMLST,$ST_H1,$ST_H2,$ST_O_antigen,$SEROVAR,$RESFINDER_GENES,$RESFINDER_RESISTANT,$RESFINDER_MUTS" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
				echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",$C3_warning,\"$CONTROL_3_CONT\",$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,\"$CONTROL4_CHECKM\",\"$vCONTROL4_CHECKM_quality\",\"$STVAR\",\"$MLST\",\"$STANTI\",\"$STPATHO\",\"$ST_CGMLST\",\"$ST_H1\",\"$ST_H2\",\"$ST_O_antigen\",\"$SEROVAR\",\"$RESFINDER_GENES\",\"$RESFINDER_RESISTANT\",\"$RESFINDER_MUTS\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"

                            elif [[ "$VAR_C2_SPE" = "Escherichia coli" ]]; then

                                #---------------------------------------------------------
                                # 05.0 Declaramos variables que no recogemos en esta especie
                                #---------------------------------------------------------
                                ST_CGMLST="-"
                                ST_H1="-"
                                ST_H2="-"
                                ST_O_antigen="-"
                                SEROVAR="-"

                                #---------------------------------------------------------
                                # 05.1 Serotipado y patotipado
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/seq_typing
                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/patho_typing

                                #sample="17_STEC_00757"

                                # SEQ_typing
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/seq_typing:2.3-dev-2021-03 \
                                seq_typing.py reads \
                                --org Escherichia coli \
                                --fastq ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                -o ${OUTPUT_PATH}/3_typing/PLUS/seq_typing \
                                -j ${THREADS} >> ${OUTPUT_PATH}/log/seq_typing/${sample}_seq_typing.log 2>&1

                                # Renombramos SEQ_typing
                                # Lista de archivos a verificar
                                fp2rename=(
                                    "seq_typing.report.txt"
                                    "seq_typing.report_types.tab"
                                )

                                # Verificar cada archivo en la lista FILES
                                for fp2 in "${fp2rename[@]}"; do
                                    if [ -e "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${fp2}" ]; then
                                        mv "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${fp2}" "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${sample}_${fp2}"
                                    fi
                                done
                                # Verificar y renombrar archivos que coinciden con la expresión regular
                                log2rename=$(find ${OUTPUT_PATH}/3_typing/PLUS/seq_typing -type f -name "run.*.log" 2>/dev/null)

                                if [ -e "$log2rename" ]; then
                                    mv "${log2rename}" "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${sample}_$(basename ${log2rename})"
                                fi

                                if [ -f "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${sample}_seq_typing.report.txt" ]; then
                                    STVAR_1=$(head -n 1 "${OUTPUT_PATH}/3_typing/PLUS/seq_typing/${sample}_seq_typing.report.txt")
                                else
                                    STVAR_1="-"
                                fi

                                # PATHO_typing
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/patho_typing:0.3.0-1 \
                                patho_typing.py \
                                -s Escherichia coli \
                                -f ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                -o ${OUTPUT_PATH}/3_typing/PLUS/patho_typing \
                                -j ${THREADS} >> ${OUTPUT_PATH}/log/patho_typing/${sample}_patho_typing.log 2>&1

                                # Renombramos PATHO_typing

                                # Lista de archivos a verificar
                                fp2rename=(
                                    "patho_typing.extended_report.txt"
                                    "patho_typing.report.txt"
                                    "rematch"
                                )

                                # Verificar cada archivo en la lista FILES
                                for fp2 in "${fp2rename[@]}"; do
                                    if [ -e "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${fp2}" ]; then
                                        mv "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${fp2}" "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${sample}_${fp2}"
                                    fi
                                done
                                # Verificar y renombrar archivos que coinciden con la expresión regular
                                log2rename=$(find ${OUTPUT_PATH}/3_typing/PLUS/patho_typing -type f -name "run.*.log" 2>/dev/null)

                                if [ -e "$log2rename" ]; then
                                    mv "${log2rename}" "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${sample}_$(basename ${log2rename})"
                                fi

                                if [ -f "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${sample}_patho_typing.report.txt" ]; then
                                    STVAR_2=$(head -n 1 "${OUTPUT_PATH}/3_typing/PLUS/patho_typing/${sample}_patho_typing.report.txt")
                                else
                                    STVAR_2="-"
                                fi

                                STANTI=$STVAR_1
                                STPATHO=$STVAR_2

                                #---------------------------------------------------------
                                # 05.2 Detección de AMR
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/AMR
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} genomicepidemiology/resfinder:4.3.3 \
                                -ifa ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                -o ${OUTPUT_PATH}/3_typing/AMR/${sample} -s ecoli --acquired --point \
                                >> ${OUTPUT_PATH}/log/resfinder/${sample}_resfinder.log 2>&1


                                FILE_RESISTANT="${OUTPUT_PATH}/3_typing/AMR/${sample}/pheno_table.txt"
                                if [[ -f "$FILE_RESISTANT" ]]; then
                                    RESFINDER_RESISTANT=$(awk -F'\t' '{gsub(/ /, "_", $1); gsub(/ /, "_", $2)} $3 == "Resistant" && $4 == 3 {print $1"("$2")"}' "$FILE_RESISTANT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_RESISTANT="-"
                                fi

                                FILE_RESGENE="${OUTPUT_PATH}/3_typing/AMR/${sample}/ResFinder_results_tab.txt"
                                if [[ -f "$FILE_RESGENE" ]]; then
                                    RESFINDER_GENES=$(awk -F'\t' '$4 ~ /^[0-9]+(\.[0-9]+)?$/ && $4 >= 90 {gsub(/ /, "_", $1); print $1}' "$FILE_RESGENE" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_GENES="-"
                                fi

                                FILE_RESMUT="${OUTPUT_PATH}/3_typing/AMR/${sample}/PointFinder_results.txt"
                                if [[ -f "$FILE_RESMUT" ]]; then
                                    RESFINDER_MUTS=$(awk -F'\t' 'NR>1 {gsub(/ /, "_", $1); gsub(/ /, "_", $4); gsub(/,_/, ";", $4); print $1"("$4")"}' "$FILE_RESMUT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_MUTS="-"
                                fi

                                [ -z "$RESFINDER_GENES" ] && RESFINDER_GENES="-"
                                [ -z "$RESFINDER_RESISTANT" ] && RESFINDER_RESISTANT="-"
                                [ -z "$RESFINDER_MUTS" ] && RESFINDER_MUTS="-"
                                # docker run -it staphb/ncbi-amrfinderplus /bin/bash
                                # amrfinder -p ecoli.faa --plus -o AMRFinder_complete.tsv --threads 4 --ident_min $(echo "scale=2; 90/100" | bc -l ) \
                                # --coverage_min $(echo "scale=2; 80/100" | bc -l ) --name ecoli --protein_output ecoli_args.faa --database /amrfinder/data/2023-04-17.1 
                                # awk -F '	' '{ if ($3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;

                                #---------------------------------------------------------
                                # 05.3 MLST y perfil de virulencia
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/MLST

                                # stringMLST.py \
                                # --predict \
                                # -1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                                # -2 ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                # -p \
                                # --prefix ${RESOURCES}/string_mlst_db/Escherichia_coli_1/Escherichia_coli_1 \
                                # -k 35 \
                                # -o ${OUTPUT_PATH}/3_typing/MLST/${sample}_stringmlst_ecoli1.txt

                                # stringMLST.py \
                                # --predict \
                                # -1 ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz \
                                # -2 ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                # -p \
                                # --prefix ${RESOURCES}/string_mlst_db/Escherichia_coli_2/Escherichia_coli_2 \
                                # -k 35 \
                                # -o ${OUTPUT_PATH}/3_typing/MLST/${sample}_stringmlst_ecoli2.txt

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/mlst:2.23.0 \
                                mlst \
                                --threads ${THREADS} \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out \
                                2> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.log

                                # Limpiamos un poco el fichero
                                var2erase=${OUTPUT_PATH}/2_assembly/
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out 

                                FILE_MLST=${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out
                                if [[ -f "$FILE_MLST" ]]; then
				    MLST=$(awk -F'\t' '{result=""; for(i=4; i<=NF; i++) result = result $i ";"; sub(/;$/, "", result); gsub(/,/, "-", result); print result}' "$FILE_MLST")
				    STVAR=$(awk '{if(NR==1) print $3}' "$FILE_MLST")
                                else
                                    MLST="-"
                                    STVAR="-"
                                fi

                                # ABRICATE
                                ##########
                                mkdir -p ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)"
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                2> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.log

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                --summary ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                > ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # Limpiamos un poco el fichero
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # awk '{print $6}' ${OUTPUT_PATH}"/typing/"${sample}"/MLST_VIRULENCE-PROFILE/"${sample}"_abricate.out" | uniq | grep -v "GENE" | awk 'BEGIN { ORS = " " } { print }'
                                mkdir -p ${OUTPUT_PATH}/3_typing/PLUS/ariba/${sample}
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHARIBA}:${PATHARIBA} staphb/ariba:latest \
                                ariba run \
                                --threads ${THREADS} --force \
                                ${PATHARIBA}/Escherichia_coli1/ref_db \
                                ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                ${OUTPUT_PATH}/3_typing/PLUS/ariba/${sample}_ecoli1 \
			                    > ${OUTPUT_PATH}/log/ariba/${sample}_ecoli1.log 2>&1


                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHARIBA}:${PATHARIBA} staphb/ariba:latest \
                                ariba run \
                                --threads ${THREADS} --force \
                                ${PATHARIBA}/Escherichia_coli2/ref_db \
                                ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                ${OUTPUT_PATH}/3_typing/PLUS/ariba/${sample}_ecoli2 \
			                    > ${OUTPUT_PATH}/log/ariba/${sample}_ecoli2.log 2>&1

                                #---------------------------------------------------------
                                # 05.4  Allele calling
                                #---------------------------------------------------------
                                #PATHcgMLST=/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/cgMLST_data

                                chewBBACA_BLAST_SCORE_RATIO=0.6
                                chewBBACA_MINIMUM_LENGTH=0
                                chewBBACA_TRANSLATION_TABLE=11
                                chewBBACA_SIZE_THRESHOLD=0.2
                                BACT=Escherichia_coli

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py AlleleCall \
                                --cpu ${THREADS} --force-continue \
                                -i ${OUTPUT_PATH}/tmp/pipeline/${sample}/BACTpipe_results/shovill/ \
                                --bsr ${chewBBACA_BLAST_SCORE_RATIO} \
                                --l ${chewBBACA_MINIMUM_LENGTH} \
                                --t ${chewBBACA_TRANSLATION_TABLE} \
                                --st ${chewBBACA_SIZE_THRESHOLD} \
                                -g ${PATHcgMLST}/schema/${BACT} \
                                --ptf ${PATHcgMLST}/schema/${BACT}/${BACT}.trn \
                                -o ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_1_allelecall.log 2>&1

                                # Get PATH results
                                PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}/$(ls ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample})
                                # PATHALLELE=${OUTPUT_PATH}/3_typing/ALLELE_CALLING/${sample}
                                # Paralog detection
                                ###################
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHcgMLST}:${PATHcgMLST} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py RemoveGenes -i ${PATHALLELE}/results_alleles.tsv -g ${PATHALLELE}/paralogous_counts.tsv -o ${PATHALLELE}/results_alleles_NoParalogs.tsv \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_2_removegenes.log 2>&1

                                # ExtractCgMLST - Determine the set of loci that constitute the core genome
                                ###################
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} ummidock/chewbbaca:3.1.2 \
                                chewBBACA.py ExtractCgMLST -i ${PATHALLELE}/results_alleles.tsv -o ${PATHALLELE} \
                                > ${OUTPUT_PATH}/log/chewbbaca/${sample}_3_extractgcmlst.log 2>&1

                                # rm -r ${OUTPUT_PATH}/3_typing/ALLELE_CALLING/tmp_${sample}

                                #echo -e "$sample,$fastq_1,$fastq_2,$VAR_C2_SPE,$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,$CONTROL_1_Q30,$C2_SPE,$CONTROL_2_BACT,$C3_warning,$CONTROL_3_CONT,$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,$CONTROL4_CHECKM,$vCONTROL4_CHECKM_quality,$STVAR,$MLST,$STANTI,$STPATHO,$ST_CGMLST,$ST_H1,$ST_H2,$ST_O_antigen,$SEROVAR,$RESFINDER_GENES,$RESFINDER_RESISTANT,$RESFINDER_MUTS" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
                                echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",$C3_warning,\"$CONTROL_3_CONT\",$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,\"$CONTROL4_CHECKM\",\"$vCONTROL4_CHECKM_quality\",\"$STVAR\",\"$MLST\",\"$STANTI\",\"$STPATHO\",\"$ST_CGMLST\",\"$ST_H1\",\"$ST_H2\",\"$ST_O_antigen\",\"$SEROVAR\",\"$RESFINDER_GENES\",\"$RESFINDER_RESISTANT\",\"$RESFINDER_MUTS\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"

                            elif [[ "$VAR_C2_SPE" = "Campylobacter jenuni" ]]; then

                                #---------------------------------------------------------
                                # 05.0 Declaramos variables que no recogemos en esta especie
                                #---------------------------------------------------------
                                STANTI="-"
                                STPATHO="-"
                                ST_CGMLST="-"
                                ST_H1="-"
                                ST_H2="-"
                                ST_O_antigen="-"
                                SEROVAR="-"

                                #---------------------------------------------------------
                                # 05.2 Detección de AMR
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/AMR

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} genomicepidemiology/resfinder:4.3.3 \
                                -ifa ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                -o ${OUTPUT_PATH}/3_typing/AMR/${sample} -s cjejuni --acquired --point \
                                >> ${OUTPUT_PATH}/log/resfinder/${sample}_resfinder.log 2>&1

                                FILE_RESISTANT="${OUTPUT_PATH}/3_typing/AMR/${sample}/pheno_table.txt"
                                if [[ -f "$FILE_RESISTANT" ]]; then
                                    RESFINDER_RESISTANT=$(awk -F'\t' '{gsub(/ /, "_", $1); gsub(/ /, "_", $2)} $3 == "Resistant" && $4 == 3 {print $1"("$2")"}' "$FILE_RESISTANT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_RESISTANT="-"
                                fi

                                FILE_RESGENE="${OUTPUT_PATH}/3_typing/AMR/${sample}/ResFinder_results_tab.txt"
                                if [[ -f "$FILE_RESGENE" ]]; then
                                    RESFINDER_GENES=$(awk -F'\t' '$4 ~ /^[0-9]+(\.[0-9]+)?$/ && $4 >= 90 {gsub(/ /, "_", $1); print $1}' "$FILE_RESGENE" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_GENES="-"
                                fi

                                FILE_RESMUT="${OUTPUT_PATH}/3_typing/AMR/${sample}/PointFinder_results.txt"
                                if [[ -f "$FILE_RESMUT" ]]; then
                                    RESFINDER_MUTS=$(awk -F'\t' 'NR>1 {gsub(/ /, "_", $1); gsub(/ /, "_", $4); gsub(/,_/, ";", $4); print $1"("$4")"}' "$FILE_RESMUT" | sort | uniq | paste -sd';' -)
                                else
                                    RESFINDER_MUTS="-"
                                fi

                                [ -z "$RESFINDER_GENES" ] && RESFINDER_GENES="-"
                                [ -z "$RESFINDER_RESISTANT" ] && RESFINDER_RESISTANT="-"
                                [ -z "$RESFINDER_MUTS" ] && RESFINDER_MUTS="-"
                                # docker run -it staphb/ncbi-amrfinderplus /bin/bash
                                # amrfinder -p ecoli.faa --plus -o AMRFinder_complete.tsv --threads 4 --ident_min $(echo "scale=2; 90/100" | bc -l ) \
                                # --coverage_min $(echo "scale=2; 80/100" | bc -l ) --name ecoli --protein_output ecoli_args.faa --database /amrfinder/data/2023-04-17.1 
                                # awk -F '	' '{ if ($3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;

                                #---------------------------------------------------------
                                # 05.3 MLST y perfil de virulencia
                                #---------------------------------------------------------
                                mkdir -p ${OUTPUT_PATH}/3_typing/MLST
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/mlst:2.23.0 \
                                mlst \
                                --threads ${THREADS} \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out \
                                2> ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.log

                                # Limpiamos un poco el fichero
                                var2erase=${OUTPUT_PATH}/2_assembly/
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out

                                FILE_MLST=${OUTPUT_PATH}/3_typing/MLST/${sample}_mlst.out
                                if [[ -f "$FILE_MLST" ]]; then
				    MLST=$(awk -F'\t' '{result=""; for(i=4; i<=NF; i++) result = result $i ";"; sub(/;$/, "", result); gsub(/,/, "-", result); print result}' "$FILE_MLST")
				    STVAR=$(awk '{if(NR==1) print $3}' "$FILE_MLST")
                                else
                                    MLST="-"
                                    STVAR="-"
                                fi

                                # ABRICATE
                                ##########
                                # Para ver bases de datos disponibles
                                #docker run --rm staphb/abricate:1.0.1-insaflu-220727 abricate --list

                                # insaflu	34	nucl	2022-Nov-16
                                mkdir -p ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                ${OUTPUT_PATH}/2_assembly/${sample}.contigs.fa \
                                1> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                2> ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.log

                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} staphb/abricate:1.0.1-insaflu-220727 \
                                abricate \
                                --threads ${THREADS} \
                                --db vfdb \
                                --summary ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out \
                                > ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # Limpiamos un poco el fichero
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate.out
                                sed -i "s:^$var2erase::" ${OUTPUT_PATH}/3_typing/VIRULENCE-PROFILE/${sample}_abricate_summary.tsv

                                # awk '{print $6}' ${OUTPUT_PATH}"/typing/"${sample}"/MLST_VIRULENCE-PROFILE/"${sample}"_abricate.out" | uniq | grep -v "GENE" | awk 'BEGIN { ORS = " " } { print }'
                                mkdir -p ${OUTPUT_PATH}/3_typing/ariba/${sample}
                                docker run --cpus ${THREADS} --rm -u "$(id -u)":"$(id -g)" -v ${OUTPUT_PATH}:${OUTPUT_PATH} -v ${PATHARIBA}:${PATHARIBA} staphb/ariba:latest \
                                ariba run \
                                --threads ${THREADS} --force \
                                ${PATHARIBA}/Campylobacter_jejuni/ref_db \
                                ${OUTPUT_PATH}/0_fastq/${sample}_1.fastq.gz ${OUTPUT_PATH}/0_fastq/${sample}_2.fastq.gz \
                                ${OUTPUT_PATH}/3_typing/ariba/${sample} \
			                    > ${OUTPUT_PATH}/log/ariba/${sample}.log 2>&1

                                #echo -e "$sample,$fastq_1,$fastq_2,$VAR_C2_SPE,$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,$CONTROL_1_Q30,$C2_SPE,$CONTROL_2_BACT,$C3_warning,$CONTROL_3_CONT,$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,$CONTROL4_CHECKM,$vCONTROL4_CHECKM_quality,$STVAR,$MLST,$STANTI,$STPATHO,$ST_CGMLST,$ST_H1,$ST_H2,$ST_O_antigen,$SEROVAR,$RESFINDER_GENES,$RESFINDER_RESISTANT,$RESFINDER_MUTS" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
                                echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",$C3_warning,\"$CONTROL_3_CONT\",$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,\"$CONTROL4_CHECKM\",\"$vCONTROL4_CHECKM_quality\",\"$STVAR\",\"$MLST\",\"$STANTI\",\"$STPATHO\",\"$ST_CGMLST\",\"$ST_H1\",\"$ST_H2\",\"$ST_O_antigen\",\"$SEROVAR\",\"$RESFINDER_GENES\",\"$RESFINDER_RESISTANT\",\"$RESFINDER_MUTS\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
                            fi

                        fi
                    else
                        echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        printf '%s FAILED IN CONTAMINATION CHECK\n' "$sample" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                        # Si no se ejecuta checkm, tenemos que añadir los 12 campos vacios, para no tener problemas con el número de columnas finales
			CONTROL_3_CONT="FAIL"
			vcheckm="\"-\",-,-,-,\"-\",-,-,-,-,-,-,-,-"
			echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",$C3_warning,\"$CONTROL_3_CONT\",$VAR_C4_COMPLETENESS,$VAR_C4_CONTAMINATION,$vcheckm,\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
			
                        conda deactivate
                        
                    fi
                fi
            else
                echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                printf '%s FAILED IN SPECIES CHECK\n' "$sample" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
                # Si no se ejecuta checkm, tenemos que añadir los 12 campos vacios, para no tener problemas con el número de columnas finales
                CONTROL_2_BACT="FAIL"
		vcheckm="\"-\",-,-,-,\"-\",-,-,-,-,-,-,-,-"
                echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"$C2_SPE\",\"$CONTROL_2_BACT\",-,\"-\",-,-,$vcheckm,\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
                conda deactivate
            fi

        else
            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            printf '%s FAILED IN Q30 CHECK\n' "$sample" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
            # Si no se ejecuta checkm, tenemos que añadir los 12 campos vacios, para no tener problemas con el número de columnas finales
            CONTROL_1_Q30="FAIL"
	    vcheckm="\"-\",-,-,-,\"-\",-,-,-,-,-,-,-,-"
            echo -e "\"$sample\",\"$fastq_1\",\"$fastq_2\",\"$VAR_C2_SPE\",$VAR_C3_GENOME,$VAR_C4_GENOME_MIN,$VAR_C4_GENOME_MAX,$VAR_C3_SNV,$VAR_C4_CONTIGS,$VAR_C1_LENGTH,$BASESQ30,$COVQ30,$VALUEQ30,\"$CONTROL_1_Q30\",\"-\",\"-\",-,\"-\",-,-,$vcheckm,\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\",\"-\"" >> "${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv"
            
            conda deactivate
        fi

        # RESETAMOS VARIABLES A "-"
        var_reset=(sample fastq_1 fastq_2 VAR_C2_SPE VAR_C3_GENOME VAR_C4_GENOME_MIN VAR_C4_GENOME_MAX VAR_C3_SNV VAR_C4_CONTIGS
            BASESQ30 COVQ30 VALUEQ30 CONTROL_1_Q30 C2_SPE CONTROL_2_BACT C3_warning CONTROL_3_CONT
            vcheckm CONTROL4_CHECKM vCONTROL4_CHECKM_quality STVAR STPATHO ST_CGMLST ST_H1 ST_H2 ST_O_antigen SEROVAR RESFINDER_GENES 
            RESFINDER_RESISTANT RESFINDER_MUTS MLST STANTI)
        for varr in "${var_reset[@]}"
        do
            if [[ $varr == "vcheckm" ]]; then
                vcheckm="\"-\",-,-,-,\"-\",-,-,-,-,-,-,-,-"
            else
                eval "$varr='-'"
            fi
            #unset $varr
        done

    done
} <"${SAMPLESHEET}"

conda activate ${CONDAPATH}/dgsp_efsa_sp
multiqc -o ${OUTPUT_PATH}/4_results/${DATE}_multiqc ${OUTPUT_PATH}
conda deactivate

# Reporte final
FILE_SUMMARY=${OUTPUT_PATH}/4_results/${DATE}_summary_${RUN}.csv

if [[ -f "$FILE_SUMMARY" ]]; then
    conda activate dgsp_efsa_report

    Rscript $EFSA_PATH/scripts/summary_to_excel.R $FILE_SUMMARY

    echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
    printf '%s DONE !\n' "$sample" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
    echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
else
    echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
    printf '%s FAILED FINAL REPORT\n' "$sample" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
    echo "****************************************" | tee -a ${OUTPUT_PATH}/log/efsa/${sample}.log
fi
conda deactivate

