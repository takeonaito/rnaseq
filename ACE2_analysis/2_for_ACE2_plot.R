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


# read count related files. and extract only ACE2 gene (ENSG00000130234)
count <- read_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/all.gene_counts.xls") 
tpm <- read_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myftmp.txt")
fpkm <- read_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myfpkm.txt")
cpm <- read_tsv("/mnt/share6/FOR_Takeo/RNAseq_WashU/myfcpm.txt")

types <- c("count","tpm","fpkm","cpm")

answer <- data.frame(seqid = colnames(count[,-1])) %>% 
  mutate(seqid = as.character(seqid))

for ( i in 1:length(types)) {
  count_name <- paste0("ACE2_",types[i])
  
  count_tmp <- eval(parse(text = types[i])) %>% 
    filter(ensembl_gene_id == "ENSG00000130234") %>% 
    data.frame()
  df_count <- data.frame(seqid = colnames(count_tmp[,-1])) %>% 
    mutate(UQ(rlang::sym(count_name)) := unlist(count_tmp[1,-1], use.names=FALSE)) %>% 
    mutate(seqid = as.character(seqid))
  
  answer <- answer %>% 
    left_join(df_count,by = "seqid")
  
}


# make RNAseq_ID column of target1  match that of answer 
target2 <- target1 %>% 
  mutate(RNAseq_ID = paste0("sample.",RNAseq_ID)) %>% 
  mutate(RNAseq_ID = str_replace_all(RNAseq_ID,"CC","cc")) %>% 
  mutate(RNAseq_ID = str_replace_all(RNAseq_ID,"-","_")) 

length(intersect(target2$RNAseq_ID,answer$seqid))

# join target df and answer df 
# mutate dammy Gender column for regression
target3 <- target2 %>% 
  left_join(answer,by = c("RNAseq_ID" = "seqid")) %>% 
  dplyr::select(Genetic_ID,RNAseq_ID,contains("ACE2"),Age.at.Collection,Gender) %>% 
  mutate(Age.at.Collection = as.numeric(Age.at.Collection)) %>% 
  mutate(Gender1 = ifelse(Gender == "M",1,0))

# tidy data frame for making facet plot
target4 <- target3 %>% 
  gather(contains("ACE"),key = ACE2_sort, value= count_num) 


# make facet plots for gender and age at collection with regression line.
p <- target4 %>% 
  ggplot(aes(x = Gender1 ,y = count_num,color= Gender)) + geom_beeswarm(cex = 5,
                                                                      size = 0.5) +
  geom_smooth(method = 'glm', formula = y ~x,color = "black",fullrange=TRUE) +
  scale_x_continuous(breaks=c(0,1),
                   labels=c("F", "M"), expand = c(0, 0.5))

p + facet_wrap(~ ACE2_sort,scales="free")


z <- target4 %>% 
  ggplot(aes(x = Age.at.Collection,y = count_num)) + geom_point(size = 0.5) +
  scale_x_continuous( breaks =seq(0, 74, 5))+
  geom_smooth(method = 'glm', formula = y ~x,color = "black") 

z + facet_wrap(~ ACE2_sort,scales="free_y")



