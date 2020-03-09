library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(readr)

# make data table for gene length (hg19)
hg19.ens <- makeTxDbFromUCSC(genome="hg19", tablename="ensGene")

exonic <- exonsBy(hg19.ens, by="gene")
red.exonic <- GenomicRanges::reduce(exonic) # this reduce remove duplicated region and merge it.
exon.lengths <- sum(width(red.exonic))

df <- data.frame(id = names(exon.lengths), len = exon.lengths) %>% 
  mutate(id = as.character(id))


# read 
count <- read_tsv("/home/takeo/rnaseq/WashU/data/all.gene_counts.xls")

count <- count %>% 
  left_join(df,by = c("ensembl_gene_id" = "id"))

count %>% 
  filter(is.na(len)) %>% 
  dplyr::select(ensembl_gene_id)


# calculate tpm, fpkm and cpm
tpm <- count %>% 
  drop_na() %>% 
  mutate_at(vars(-ensembl_gene_id,-len), list(~ ./len)) %>% 
  mutate_at(vars(-ensembl_gene_id,-len), list(~ (./ sum(.)*1e6)))

fpkm <- count %>% 
  drop_na() %>% 
  mutate_at(vars(-ensembl_gene_id,-len), list(~ (./ sum(.))*1e9 )) %>% 
  mutate_at(vars(-ensembl_gene_id,-len), list(~ ./len)) 

cpm <- count %>% 
  mutate_at(vars(-ensembl_gene_id,-len), list( ~ (./ sum(.))*1e6))


# compare original file and confirm concordance.
dir="/mnt/share6/SHARED_DATASETS/Gene_Expression_Transcriptomics_DATA/htcf.wustl.edu/files/mMgOrrMY/Stappenbeck_2819_12/"

otehon <- read_tsv(paste0(dir,"all.gene_CPM.xls"))



otehon1 <- otehon %>% 
  left_join(cpm,by = "ensembl_gene_id") %>% 
  drop_na()

otehon1 %>% 
  ggplot(aes(x = sample.cc_29.x, y = sample.cc_29.y)) +
  geom_point()


# export each file

fpkm %>%  write_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myfpkm.txt")
cpm %>% write_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myfcpm.txt")
tpm %>% write_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myftmp.txt")


