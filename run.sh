# bash /home/RNA-seq_data/codes/run.sh GSE.../E-MTAB-... MesAur/GRCh38/GRCm39

pjname=$1
#MesAur, GRCh38, GRCm39
strain=$2

mkdir -p /home/C/working/${pjname}
cd /home/C/working/${pjname}

mkdir -p /home/D/${pjname}

if [[ ${pjname} =~ ^E-MTAB.* ]]
then
    wget https://www.ebi.ac.uk/arrayexpress/files/${pjname}/${pjname}.sdrf.txt
    awk -F "\t" -vcol=$(head -n1 ${pjname}.sdrf.txt | xargs -n1 -d\\t | sed -n "/Comment\[FASTQ_URI\]/=") 'NR>1{print $col}' ${pjname}.sdrf.txt | sort > URL.list

    for URL in `cat URL.list`
    do

        ERRfile=${URL##*/}
        ERR=${ERRfile%%.*}

        if [[ ${ERR} =~ .*_1 ]]
        then
            echo 1
            ERR1=${ERR}
            ERRfile1=${ERRfile}
            URL1=${URL}
        elif [[ ${ERR} == ${ERR1%%_1}_2 ]]
        then
            echo 2
            ERRfile2=${ERRfile}
            ERR2=${ERR}
            ERR=${ERR1%%_1}
            echo ${ERR} >> RUN.list

            URL2=${URL}

            mkdir ${ERR}
            cd ${ERR}
            wget ${URL1} ${URL2}
            
            fastp -i  ${ERRfile1}  -I ${ERRfile2} -o ${ERR1}_fastp.fastq.gz -O ${ERR2}_fastp.fastq.gz -w 16 -q 15 -n 10
            unpigz ${ERR1}_fastp.fastq.gz ${ERR2}_fastp.fastq.gz
            STAR --genomeDir /home/C/references/${strain}/STAR_index_${strain}  \
                --readFilesIn ${ERR1}_fastp.fastq ${ERR2}_fastp.fastq \
                --genomeLoad NoSharedMemory \
                --outFilterMultimapNmax 1 \
                --runThreadN 16 \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix  ${ERR}_
                # --limitBAMsortRAM 540000000000
            
            samtools sort ${ERR}_Aligned.out.bam > ${ERR}_Aligned.sortedByCoord.out.bam
            rm  ${ERR}_Aligned.out.bam

            rsem-calculate-expression --alignments -p 16 --paired-end \
                ${ERR}_Aligned.toTranscriptome.out.bam \
                /home/C/references/${strain}/RSEM_index_${strain}/RSEM_index_${strain} \
                ${ERR}_RSEM
            multiqc .
            cd ..

            mv ${ERR} /home/D/${pjname}
        else
            echo 3
            echo ${ERR} >> RUN.list

            mkdir ${ERR}
            cd ${ERR}
            wget ${URL}


            fastp -i  ${ERRfile} -o ${ERR}_fastp.fastq.gz -w 16 -q 15 -n 10
            unpigz ${ERR}_fastp.fastq.gz
            STAR --genomeDir /home/C/references/${strain}/STAR_index_${strain}  \
                --readFilesIn ${ERR}_fastp.fastq \
                --genomeLoad NoSharedMemory \
                --outFilterMultimapNmax 1 \
                --runThreadN 16 \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix  ${ERR}_
                # --limitBAMsortRAM 540000000000
            
            samtools sort ${ERR}_Aligned.out.bam > ${ERR}_Aligned.sortedByCoord.out.bam
            rm  ${ERR}_Aligned.out.bam

            rsem-calculate-expression --alignments -p 16 \
                ${ERR}_Aligned.toTranscriptome.out.bam \
                /home/C/references/${strain}/RSEM_index_${strain}/RSEM_index_${strain} \
                ${ERR}_RSEM
            multiqc .
            cd ..

            mv ${ERR} /home/D/${pjname}
        fi
    done

elif [[ ${pjname} =~ ^GSE.* ]]
then
    esearch -db gds -query " ${pjname}[ACCN] AND GSM[ETYP]" |elink -target SRA |efetch -format runinfo > ${pjname}_metadata.csv

    IFS=$'\n'
    for i in `cat ${pjname}_metadata.csv | awk -F ',' 'NR>1'`
    do
        SRR=`echo $i | cut -d ',' -f 1`
        LAYOUT=`echo $i | cut -d ',' -f 16`
        
        echo ${SRR} >> RUN.list

        prefetch ${SRR}
        cd ${SRR}
        fasterq-dump -e 16 -p --split-files ${SRR}.sra
        if [[ "${LAYOUT}" == "SINGLE" ]]
        then
            fastp -i  ${SRR}.fastq -o ${SRR}_fastp.fastq.gz -w 16 -q 15 -n 10
            unpigz ${SRR}_fastp.fastq.gz
            STAR --genomeDir /home/C/references/${strain}/STAR_index_${strain}  \
                --readFilesIn ${SRR}_fastp.fastq \
                --genomeLoad NoSharedMemory \
                --outFilterMultimapNmax 1 \
                --runThreadN 16 \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix  ${SRR}_
                # --limitBAMsortRAM 540000000000

            samtools sort ${SRR}_Aligned.out.bam > ${SRR}_Aligned.sortedByCoord.out.bam
            rm  ${SRR}_Aligned.out.bam
            rsem-calculate-expression --alignments -p 16 \
                ${SRR}_Aligned.toTranscriptome.out.bam \
                /home/C/references/${strain}/RSEM_index_${strain}/RSEM_index_${strain} \
                ${SRR}_RSEM
        elif [[ "${LAYOUT}" == "PAIRED" ]]
        then
            fastp -i  ${SRR}_1.fastq  -I ${SRR}_2.fastq -o ${SRR}_1_fastp.fastq.gz -O ${SRR}_2_fastp.fastq.gz -w 16 -q 15 -n 10
            unpigz ${SRR}_1_fastp.fastq.gz ${SRR}_2_fastp.fastq.gz
            STAR --genomeDir /home/C/references/${strain}/STAR_index_${strain}  \
                --readFilesIn ${SRR}_1_fastp.fastq ${SRR}_2_fastp.fastq \
                --genomeLoad NoSharedMemory \
                --outFilterMultimapNmax 1 \
                --runThreadN 16 \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix  ${SRR}_
                # --limitBAMsortRAM 540000000000

            samtools sort ${SRR}_Aligned.out.bam > ${SRR}_Aligned.sortedByCoord.out.bam
            rm  ${SRR}_Aligned.out.bam
            rsem-calculate-expression --alignments -p 16 --paired-end \
                ${SRR}_Aligned.toTranscriptome.out.bam \
                /home/C/references/${strain}/RSEM_index_${strain}/RSEM_index_${strain} \
                ${SRR}_RSEM
        else
            echo "Layout error"
        fi
        multiqc .
        cd ..

        mv ${SRR} /home/D/${pjname}
    done
else
    echo "Invalid name"
fi

# ls -l /home/C/working/${pjname} | grep ^- | awk '{print $9}' | xargs -I {} mv {} /home/D/RNA-seq_archives/${pjname}
# mv -r /home/C/working/${pjname}/* /home/D/RNA-seq_archives/${pjname}
rm -r /home/C/working/${pjname}

echo "result.csv exported to /mnt/d/RNA-seq_data/${pjname}"