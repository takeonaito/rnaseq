library(tidyverse)
library(data.table)
library(readxl)
library(readr)
library(rlang)
library(ggbeeswarm)

# read id link file of Washi U rnaseq
sample <- read_xlsx("/mnt/share6/FOR_Takeo/RNAseq_WashU/WashU_BMI_RNAseq_IDlink.xlsx")
sample$Genetic_ID <- str_replace(sample$Genetic_ID,"10-0441/10-1045","10-0441")

# extract necessary column 
target <- sample %>%
  dplyr::select(Genetic_ID,RNAseq_ID)


# read necessary files (disease type)

disease <- read_xls("/mnt/share6/FOR_Takeo/phenotypdata/Copy of Genetics 01_02_2019.xls",
                    col_types = "text")
colnames(disease) <- make.names(colnames(disease))



# merge key file and disease file 
target1 <- target %>% 
  left_join(disease,by = c("Genetic_ID" = "Genetic.ID")) %>% 
  data.frame()





wasU <- read_tsv('/mnt/share6/FOR_Takeo/RNAseq_WashU/all.transcript_counts.xls')



ACE2 <- c('ENST00000427411', 'ENST00000252519', 'ENST00000471548', 
'ENST00000473851','ENST00000484756')


wasU <- read_tsv('/mnt/share6/FOR_Takeo/RNAseq_WashU/all.transcript_counts.xls')

answer <- data.frame(seqid = colnames(wasU[,-1])) %>% 
  mutate(seqid = as.character(seqid))

count_tmp <- wasU %>% 
  filter(ensembl_transcript_id %in% ACE2) %>% 
  rownames_to_column() %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) 

colnames(count_tmp) <- as.vector(count_tmp[1,])


count_tmp <- count_tmp %>% 
  dplyr::slice(-1)

# make RNAseq_ID column of target1  match that of answer 
target2 <- target1 %>% 
  mutate(RNAseq_ID = paste0("sample.",RNAseq_ID)) %>% 
  mutate(RNAseq_ID = str_replace_all(RNAseq_ID,"CC","cc")) %>% 
  mutate(RNAseq_ID = str_replace_all(RNAseq_ID,"-","_")) 

target3 <- target2 %>% 
  left_join(count_tmp,by = c('RNAseq_ID' = 'ensembl_transcript_id'))

target4 <- target3 %>% 
  dplyr::select(Genetic_ID,Age,Gender,contains('ENST')) %>% 
  gather(t_name,count_num,-Genetic_ID,-Age,-Gender) %>% 
  mutate(count_num = as.numeric(count_num)) %>% 
  mutate(Age = as.numeric(Age))


p <- target4 %>% 
  ggplot(aes(x = Gender ,y = count_num ,color= Gender)) + geom_beeswarm(cex = 1,
                                                                        size = 0.5) +
  geom_smooth(method = 'glm', formula = y ~x,color = "black",fullrange=TRUE) 

p + facet_wrap(~ t_name,scales="free")


p <- target4 %>% 
  ggplot(aes(x = Age ,y = count_num )) + geom_point() 

p + facet_wrap(~ t_name,scales="free")




