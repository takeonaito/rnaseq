library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(data.table)
# make data table for gene length (hg19)

hg19.ens <- makeTxDbFromUCSC(genome="hg19", tablename="ensGene")

## extract transcript and gene ID 
trans <- transcriptsBy(hg19.ens,by = "gene" )

## convert trans into dataframe
df <- as.data.frame(trans)


## extract gene ID 
genes <- df %>% 
  distinct(group_name) %>% 
  pull(group_name)


## read Ensmble DB objects()
## from '~/rnaseq/WashU/R/processing_dataset/for_disese_duration_and_onsetage.R'
DB <- readRDS("~/rnaseq/WashU/data/DB.obj")
EDB <- EnsDb(DB)

symbol_df <- ensembldb::select(EDB, key=genes,
                            columns=c("SYMBOL","SEQNAME"), 
                            keytype="GENEID")


## merge symbol_df and df
df1 <- df %>% 
  left_join(symbol_df,by = c("group_name" = "GENEID")) %>% 
  dplyr::select(tx_name,SYMBOL)

## export the result


df1 %>% 
  write_tsv("/mnt/share6/FOR_Takeo/temporary/tx_gene_anot.txt")





symbol <- mapIds(org.Hs.eg.db,
       keys=genes,
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first")




df_ano <- data.table(gene = symbol, id = names(symbol))


df1 <- df %>% 
  left_join(df_ano,by = c("group_name" = "id")) %>% 
  dplyr::select(tx_name,gene)


df1 %>% 
  write_tsv("/mnt/share6/FOR_Takeo/temporary/tx_gene_anot.txt")



