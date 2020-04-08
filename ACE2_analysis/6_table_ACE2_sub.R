
## read variant information
vari <- read_xlsx("/mnt/share6/FOR_Takeo/WES/candidate_genes/ACE2_TMPRSS2_sub_rev.xlsx")

## extract ACE2 variants and tidy data and export
vari %>% 
  filter(grepl("ACE2",genes)) %>% 
  distinct(allele_locus,.keep_all = T) %>% 
  dplyr::select(allele_locus,CADD,MAF_gnomad,vari_function,AA_change) %>% 
  separate(allele_locus,into = c("CHR","POS"),sep = ":") %>% 
  separate(POS, into = c("pos","REF","ALT"), sep = "_") %>% 
  arrange(pos) %>% 
  write_xlsx('/mnt/share6/FOR_Takeo/WES/candidate_genes/ACE2_vari.xlsx')



## read phenotype file of high CADD ACE2 subjects
pheno <- read_xlsx('/mnt/share6/FOR_Takeo/phenotypdata/ACE2_Variants.xlsx',skip = 2)
colnames(pheno) <- make.names(colnames(pheno))


## merge pheno and vari and export
pheno %>% 
  left_join(vari,by = c("Genetic.ID" = "GeneticID")) %>% 
  write_xlsx('/mnt/share6/FOR_Takeo/WES/candidate_genes/ACE2_pheno_geno.xlsx')
