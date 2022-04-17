FROM debian:latest

WORKDIR /home
RUN apt update -y && apt upgrade -y
ENV TZ=Asia/Tokyo

RUN apt install -y r-base python3 wget pigz
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm -r Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda3/bin:$PATH

RUN conda update -y conda
RUN conda config --add channels conda-forge
RUN conda config --add channels defaults
RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar xvzf v1.3.3.tar.gz
RUN cd RSEM-1.3.3 && make install

RUN conda install -y STAR samtools sra-tools multiqc fastp entrez-direct
RUN conda update -y --all

# RUN mkdir -p references/HSapiens && cd references/HSapiens && \
#     wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
#         http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz && \ 
#     unpigz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.106.gtf.gz 
# RUN cd references/HSapiens && mkdir STAR_index_HSapiens && \
#     STAR --runMode genomeGenerate \
#         --genomeDir STAR_index_HSapiens \
#         --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#         --runThreadN 16 \
#         --sjdbGTFfile Homo_sapiens.GRCh38.106.gtf && \
#     mkdir RSEM_index_HSapiens && \
#     rsem-prepare-reference -p 16 \
#         --gtf Homo_sapiens.GRCh38.106.gtf \
#         Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#         RSEM_index_HSapiens/RSEM_index_HSapiens
# 
# RUN mkdir -p references/MMusculus && cd references/MMusculus && \
#     wget http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
#         http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz && \
#     unpigz Mus_musculus.GRCm39.dna.primary_assembly.fa.gz Mus_musculus.GRCm39.106.gtf.gz
# RUN cd references/MMusculus && mkdir STAR_index_MMusculus && \
#     STAR --runMode genomeGenerate \
#         --genomeDir STAR_index_MMusculus \
#         --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \
#         --runThreadN 16 \
#         --sjdbGTFfile Mus_musculus.GRCm39.106.gtf && \
#     mkdir RSEM_index_MMusculus && \
#     rsem-prepare-reference -p 16 \
#         --gtf Mus_musculus.GRCm39.106.gtf \
#         Mus_musculus.GRCm39.dna.primary_assembly.fa \
#         RSEM_index_MMusculus/RSEM_index_MMusculus
# 
# RUN mkdir -p references/MAuratus && cd references/MAuratus && \
#     wget http://ftp.ensembl.org/pub/release-106/fasta/mesocricetus_auratus/dna/Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa.gz \
#         http://ftp.ensembl.org/pub/release-106/gtf/mesocricetus_auratus/Mesocricetus_auratus.MesAur1.0.106.gtf.gz && \
#         unpigz Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa.gz Mesocricetus_auratus.MesAur1.0.106.gtf.gz
# RUN cd references/MAuratus && mkdir STAR_index_MAuratus && \
#     STAR --runMode genomeGenerate \
#         --genomeDir STAR_index_MAuratus \
#         --genomeFastaFiles Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa \
#         --runThreadN 16 \
#         --sjdbGTFfile Mesocricetus_auratus.MesAur1.0.106.gtf && \
#     mkdir RSEM_index_MAuratus && \
#     rsem-prepare-reference -p 16 \
#         --gtf Mesocricetus_auratus.MesAur1.0.106.gtf \
#         Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa \
#         RSEM_index_MAuratus/RSEM_index_MAuratus

RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('tidyverse')"
RUN R -e "BiocManager::install('biomaRt')"
