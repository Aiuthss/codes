library(biomaRt)
library(tidyverse)

args <- commandArgs(trailingOnly = T)
setwd(paste0("D:/RNA-seq_data/", args[1]))
RUN.list <- read.table("RUN.list", header = 1)[[1]]
df = read_table(paste0(RUN.list[1], "/", RUN.list[1], "_RSEM.genes.results"))[c('gene_id', 'TPM')]
colnames(df) = c('gene_id', RUN.list[1])

for (i in RUN.list[2:length(RUN.list)]){
  df_tmp = read_table(paste0(i, "/", i, "_RSEM.genes.results"))[c('gene_id', 'TPM')]
  colnames(df_tmp) = c('gene_id', i)
  df = full_join(df, df_tmp, by='gene_id')
}
head(df)

if (args[2]=="MAuratus"){
  dataset="mauratus_gene_ensembl"
} else if (args[2]=="HSapiens"){
  dataset="hsapiens_gene_ensembl"
} else if (args[2]=="MMusculus"){
  dataset="mmusculus_gene_ensembl"
} else {
  print("invalid strain specified!")
}
db <- useMart("ensembl")
hd <- useDataset(dataset, mart = db)
res <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
             filters = "ensembl_gene_id", values = df$gene_id, 
             mart = hd, useCache = FALSE)

result <- right_join(res, df, by=c("ensembl_gene_id"="gene_id"))
write_csv(result, "results.csv")