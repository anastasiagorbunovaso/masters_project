molecular.subtypes <- PanCancerAtlas_subtypes()
molecular.subtypes <- molecular.subtypes[molecular.subtypes$cancer.type=="BRCA",]
molecular.subtypes <- molecular.subtypes %>%
  mutate(pan.samplesID = str_trunc(pan.samplesID, 16, ellipsis = "", side = 'right'))
molecular.subtypes <- molecular.subtypes[ , (names(molecular.subtypes) %in% c('pan.samplesID','Subtype_mRNA'))]


load('E:/meDNA_BRCA.RData')

meDNA <- as.data.frame(t(meDNA$X))  

meDNA <- meDNA %>%
  mutate(ROWNAMES = rownames(meDNA))  %>%
  mutate(pan.samplesID1 = str_trunc(ROWNAMES, 28, ellipsis = "", side = 'left')) %>%
  mutate(pan.samplesID = str_trunc(pan.samplesID1, 16, ellipsis = "", side = 'right'))

meDNA <- meDNA[ , !(names(meDNA) %in% c("ROWNAMES",'pan.samplesID1'))] %>%
  inner_join(molecular.subtypes)

labels <- meDNA$Subtype_mRNA

meDNA <- meDNA[ , !(names(meDNA) %in% c("pan.samplesID", 'cancer.type', 'Subtype_mRNA',
                                     'Subtype_DNAmeth','Subtype_protein','Subtype_miRNA',
                                     'Subtype_CNA','Subtype_Integrative','Subtype_other',
                                     'Subtype_Selected'))]

meDNA <- meDNA[ ,colSums(is.na(meDNA)) == 0]

sd_meDNA <- as.data.frame(as.data.frame(apply(meDNA, MARGIN = 2, FUN = sd)) > 0)

meDNA <- meDNA %>% rbind(t(sd_meDNA))

meDNA <- meDNA[, tail(meDNA, 1) != 0]

meDNA <- meDNA[-nrow(meDNA), ]
rm(sd_meDNA)

fwrite(meDNA,paste0(data_dir,'data_meDNA_BRCA.csv'))
fwrite(as.data.frame(labels),paste0(data_dir,'labels_meDNA_BRCA.csv'))

load('data/mRNA_TCGA.RData')

mRNA <- as.data.frame(t(mRNA$data))  

mRNA <- mRNA %>%
  mutate(ROWNAMES = rownames(mRNA))  %>%
  filter(grepl("BRCA", ROWNAMES, fixed = TRUE)) %>%
  mutate(pan.samplesID = str_trunc(ROWNAMES, 16, ellipsis = "", side = 'left'))

mRNA <- mRNA[ , !(names(mRNA) %in% c("ROWNAMES",'pan.samplesID1'))] %>%
  inner_join(molecular.subtypes)

labels <- mRNA$Subtype_mRNA

mRNA <- mRNA[ , !(names(mRNA) %in% c("pan.samplesID", 'cancer.type', 'Subtype_mRNA',
                                     'Subtype_DNAmeth','Subtype_protein','Subtype_miRNA',
                                     'Subtype_CNA','Subtype_Integrative','Subtype_other',
                                     'Subtype_Selected'))]

mRNA <- mRNA[ ,colSums(is.na(mRNA)) == 0]

sd_mRNA <- as.data.frame(as.data.frame(apply(mRNA, MARGIN = 2, FUN = sd)) > 0)

mRNA <- mRNA %>% rbind(t(sd_mRNA))

mRNA <- mRNA[, tail(mRNA, 1) != 0]

mRNA <- mRNA[-nrow(mRNA), ]
rm(sd_mRNA)



fwrite(mRNA,paste0(data_dir,'data_mRNA_BRCA.csv'))
fwrite(as.data.frame(labels),paste0(data_dir,'labels_mRNA_BRCA.csv'))


load('data/miRNA_TCGA.RData')

miRNA <- as.data.frame(t(miRNA$data))  

miRNA <- miRNA %>%
  mutate(ROWNAMES = rownames(miRNA))  %>%
  filter(grepl("BRCA", ROWNAMES, fixed = TRUE)) %>%
  mutate(pan.samplesID = str_trunc(ROWNAMES, 16, ellipsis = "", side = 'left'))

miRNA <- miRNA[ , !(names(miRNA) %in% c("ROWNAMES",'pan.samplesID1'))] %>%
  inner_join(molecular.subtypes)

labels <- miRNA$Subtype_mRNA

miRNA <- miRNA[ , !(names(miRNA) %in% c("pan.samplesID", 'cancer.type', 'Subtype_mRNA',
                                        'Subtype_DNAmeth','Subtype_protein','Subtype_miRNA',
                                        'Subtype_CNA','Subtype_Integrative','Subtype_other',
                                        'Subtype_Selected'))]

miRNA <- miRNA[ ,colSums(is.na(miRNA)) == 0]

sd_miRNA <- as.data.frame(as.data.frame(apply(miRNA, MARGIN = 2, FUN = sd)) > 0)

miRNA <- miRNA %>% rbind(t(sd_miRNA))

miRNA <- miRNA[, tail(miRNA, 1) != 0]

miRNA <- miRNA[-nrow(miRNA), ]
rm(sd_miRNA)

fwrite(miRNA,paste0(data_dir,'data_miRNA_BRCA.csv'))
fwrite(as.data.frame(labels),paste0(data_dir,'labels_miRNA_BRCA.csv'))