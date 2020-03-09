
library(ggrepel)
library(DESeq2)
library(BiocParallel)
library(tidyverse)
library(readr)
library(data.table)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggbeeswarm)

# read necessary files (count, sample and serology)
count <- read_tsv("/home/takeo/rnaseq/WashU/data/all.gene_counts.xls")
sample <- read_xlsx("/home/takeo/rnaseq/WashU/data/WashU_BMI_RNAseq_IDlink.xlsx")
sample$Genetic_ID <- str_replace(sample$Genetic_ID,"10-0441/10-1045","10-0441")
serology<- read_tsv("/home/takeo/rnaseq/WashU/data/serology_updated_dalin.txt")

# read necessary files (disease type and phenotype)
location <- read_tsv("/home/takeo/rnaseq/WashU/data/cd_clean.txt")
disease <- read_xls("/home/takeo/rnaseq/WashU/data/Copy of Genetics 01_02_2019.xls")
colnames(disease) <- make.names(colnames(disease))

# merge disease type and phenotype data to target file 
target <- sample %>%
  dplyr::select(Genetic_ID,RNAseq_ID)

# exclude non CD and caucasian
target <- target %>% left_join(disease,by = c("Genetic_ID" = "Genetic.ID")) %>% 
  dplyr::filter(Race == "Caucasian")

target1 <- target %>% 
  left_join(location,by = c("Genetic_ID" = "genetic_id"))

target1$RNAseq_ID <- str_replace(target1$RNAseq_ID,"-","_")

# make RNAseq_ID in sample file match to count files.
target1$RNAseq_ID <- paste0("sample.",target1$RNAseq_ID) %>% 
  base::tolower()


# confirm complete match
table( target1$RNAseq_ID %in% colnames(count))


# make sampleTable which contain serology and Genetic ID
serology1 <- serology %>% 
  drop_na(Genetic.ID)
sampleTable <- target1 %>%  
  inner_join(serology1,by = c("Genetic_ID" = "Genetic.ID")) 
## there is a duplication in cc -89( genetcid = 97-0329) --> ask dalin.
sampleTable <- sampleTable[-73,] # I excluded old data of 97-0329


# omit subjects whose serology are NA
sampleTable$hensuu <- as.numeric(sampleTable$Age.at.Collection)
sampleTable$hensuu <- as.factor(sampleTable$Gender)

sampleTable1 <- sampleTable %>% 
  drop_na(hensuu) %>% 
  dplyr::select(-Genetic_ID)


# extract count data whose serology data are available
kouho <- sampleTable1$RNAseq_ID
count1 <- count %>% 
  dplyr::select(ensembl_gene_id,kouho) %>% 
  as.data.frame()
row.names(count1) <- count1$ensembl_gene_id
count1 <- count1[,-1]

# confirm sample order match
identical(colnames(count1),sampleTable1$RNAseq_ID)

# make dds object for DESeq2

dds <- DESeqDataSetFromMatrix(countData = count1,
                              colData = sampleTable1,
                              design = ~hensuu )
keep <- rowSums(counts(dds)>10) > 10  
dds <- dds[ keep, ]

nrow(dds)


# do variance stabilizing transformation (VST) for visualization of data.
vsd <- vst(dds, blind = FALSE)

# make PCA for visualizing data set.
plotPCA(vsd, intgroup = c("hensuu"))
pcaData <- plotPCA(vsd, intgroup = c("hensuu"), returnData = TRUE)
pcaData$order <- c(1:dim(vsd)[2])
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = hensuu)) +
  geom_point(size =3)
p  + geom_text_repel(data = pcaData,aes(label = order))

# exclude outliner CC_89
exclude <- which(row.names(colData(dds)) == "sample.cc_89")
dds <- dds[,-exclude]

# do variance stabilizing transformation (VST) for visualization of data after exclusion
vsd <- vst(dds, blind = FALSE)

# make PCA for visualizing data set after exclude outliner.
plotPCA(vsd, intgroup = c("hensuu"))

# do analysis of DEG --> this will take a few minutes
register(MulticoreParam(workers = 10))

dds1 <- DESeq(dds,parallel = TRUE,minReplicatesForReplace = Inf)
res <- results(dds1,alpha = 0.05,contrast=c("hensuu","F","M"))
res <- results(dds1,alpha = 0.05)
res$ensemble = row.names(res)


res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res %>% 
  data.frame(res) %>% 
  filter(ensemble == 'ENSG00000130234')

# subset only significant genes
resSig <- subset(res, padj < 0.05) 
resSig <- subset(resSig,abs(log2FoldChange) > 1 )
resSig <- resSig[order(resSig$pvalue),]
resSig

# check the cooks distance of resSig
kouho <- which(row.names(res) %in% row.names(resSig))
round(apply(assays(dds1)[["cooks"]][kouho,],1,max),2)




# make plot of top hit gene

topGene <- rownames(res)[which.min(res$padj)]

geneCounts <- plotCounts(dds1, gene =topGene, intgroup = c("hensuu"),
                         returnData = T,normalized = T)

geneCounts <- plotCounts(dds1, gene ="ENSG00000130234", intgroup = c("hensuu"),
                         returnData = T,normalized = T)
ggplot(geneCounts, aes(x = hensuu, y = count,color = hensuu)) +  scale_y_log10() +
  geom_beeswarm(cex = 3)

ggplot(geneCounts,aes(x = hensuu, y = count)) +scale_y_log10() + geom_point()

mode(geneCounts)
