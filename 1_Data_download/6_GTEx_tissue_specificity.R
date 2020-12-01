# Load libraries
library(yarn)
library(tidyverse)

grch  <- readr::read_delim("../GRCh38/EXON_ANNOTATION", delim="\t")
gtex8 <- readRDS("../GTExV8/gtex8.rds")

# We extract the relevant annotation information from grch38
chromosome_name <- data.frame(sapply(grch$chromsome_name, as.character), stringsAsFactors = FALSE)
gene_start      <- data.frame(grch$gene_start)
gene_end        <- data.frame(grch$gene_end)
strand          <- data.frame(grch$strand)
gene_id         <- data.frame(sapply(grch$gene_id, as.character), stringsAsFactors = FALSE)
gene_name       <- data.frame(sapply(grch$gene_name,as.character), stringsAsFactors = FALSE)
gene_biotype    <- data.frame(sapply(grch$gene_biotype, as.character), stringsAsFactors = FALSE)

# Not all the gene ensembl names are overlapping between the hg38 data and the gtex data,
# so we first check which ones do, and then subset the data accordingly

intersecting_ensembl_names <- intersect(rownames(gtex8@featureData@data), grch$gene_id)

# Combine all the individual columns ( In the future, this very weird legacy code can be changed to much less lines...I just took it from Yiming)
grch38_annotation <- data.frame(chromosome_name=chromosome_name,
gene_start=gene_start, gene_end=gene_end, strand=strand, gene_name=gene_name, gene_id=gene_id, gene_biotype=gene_biotype)

colnames(grch38_annotation)[1] <- "chromosome_name"
colnames(grch38_annotation)[5] <- "gene_name"
colnames(grch38_annotation)[6] <- "gene_id"
colnames(grch38_annotation)[7] <- "gene_biotype"
rownames(grch38_annotation)    <- NULL

# Subset the annotation file by the intersecting ensembl IDs
grch38_annotation <- grch38_annotation %>%
filter(gene_id %in% intersecting_ensembl_names)

gtex8@featureData@data <- dplyr::bind_cols(gtex8@featureData@data,grch38_annotation)


# remove the duplicates and rename the columns to have the correct names
gtex8@featureData@data<-gtex8@featureData@data[,-c(1,2,3,4,5,6,7)]

colnames(gtex8@featureData@data)[1] <- "chromosome_name"
colnames(gtex8@featureData@data)[2] <- "start_position"
colnames(gtex8@featureData@data)[3] <- "end_position"
colnames(gtex8@featureData@data)[4] <- "strand"
colnames(gtex8@featureData@data)[5] <- "gene_name"
colnames(gtex8@featureData@data)[6] <- "gene_id"
colnames(gtex8@featureData@data)[7] <- "gene_biotype"

# Save the data object so we do not have to repeat this step
saveRDS(gtex8, file="../data/gtexV8_annotated_with_GRCh38.rds")


gtex8 <- readRDS("../data/gtexV8_annotated_withGRCh38.rds")

pData(gtex8)$NormGroup <- as.character(pData(gtex8)$SMTS)

pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Adipose - Subcutaneous"] <- "Adipose - Subcutaneous"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Adipose - Visceral (Omentum)"] <- "Adipose - Visceral (Omentum)"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Cells - EBV-transformed lymphocytes"] <- "Cells - EBV-transformed lymphocytes"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Whole Blood"] <- "Whole Blood"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Artery - Aorta"] <- "Artery - Aorta"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Artery - Coronary"] <- "Artery - Coronary"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Artery - Tibial"] <- "Artery - Tibial"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD %in% c("Brain - Amygdala","Brain - Anterior cingulate cortex (BA24)","Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Hippocampus","Brain - Hypothalamus","Brain - Spinal cord (cervical c-1)","Brain - Substantia nigra")] <- "Brain - Other"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD %in% c("Brain - Cerebellar Hemisphere","Brain - Cerebellum")] <- "Brain - Cerebellum"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD %in% c("Brain - Caudate (basal ganglia)","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)")] <- "Brain - Basal ganglia"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Cells - Cultured fibroblasts"] <- "Cells - Cultured fibroblasts"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Colon - Sigmoid"] <- "Colon - Sigmoid"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Colon - Transverse"] <- "Colon - Transverse"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Esophagus - Gastroesophageal Junction"] <- "Esophagus - Gastroesophageal Junction"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Esophagus - Mucosa"] <- "Esophagus - Mucosa"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Esophagus - Muscularis"] <- "Esophagus - Muscularis"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Heart - Atrial Appendage"] <- "Heart - Atrial Appendage"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD == "Heart - Left Ventricle"] <- "Heart - Left Ventricle"
#pData(gtex8)$NormGroup[pData(gtex8)$SMTSD %in% c("Kidney - Cortex","Kidney - Medulla")] <- "Kidney"
pData(gtex8)$NormGroup[pData(gtex8)$SMTSD %in% c("Skin - Not Sun Exposed (Suprapubic)","Skin - Sun Exposed (Lower leg)")] <- "Skin"

# Filter genes with not enough samples
gtex8.filtered <- filterLowGenes(gtex8, groups="NormGroup", minSamples=9)

# Tissue aware normalization based on yarn package
gtex8.filtered <- normalizeTissueAware(gtex8.filtered, groups="NormGroup")

saveRDS(gtex8.filtered, "../data/gtex8_normalises_and_filtered.rds")


# Load data
gtex8.filtered <- readRDS("../data/gtex8_normalises_and_filtered.rds")

# Calculate the mean
gtex8.mean <- by(data=t(assayData(gtex8.filtered)[["normalizedMatrix"]]),
                INDICES=pData(gtex8.filtered)$NormGroup,
                FUN=function(x) apply(x, 2, mean))

gtex8.mean <- do.call(cbind, gtex8.mean)
gtex8.mean <- cbind(gtex8.mean, gtex8.filtered@featureData@data)

# Rename the tissues so they work with R etc as variable names

colnames(gtex8.mean)[colnames(gtex8.mean) == "Adipose - Subcutaneous"] <-"Adipose_Subcutaneous"
colnames(gtex8.mean)[colnames(gtex8.mean) =="Adipose - Visceral (Omentum)" ] <-"Adipose_Visceral_Omentum"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Cells - EBV-transformed lymphocytes"] <-"Cells_EBV_transformed_lymphocytes"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Whole Blood"] <-"Whole_Blood"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Artery - Aorta"] <-"Artery_Aorta"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Artery - Coronary"] <-"Artery_Coronary"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Artery - Tibial"] <- "Artery_Tibial"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Brain - Other"] <- "Brain_Other"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Brain - Cerebellum"] <- "Brain_Cerebellum"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Brain - Basal ganglia"] <-  "Brain_Basal_ganglia"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Colon - Sigmoid"] <- "Colon_Sigmoid"
colnames(gtex8.mean)[colnames(gtex8.mean) ==  "Colon - Transverse"] <-  "Colon_Transverse"
colnames(gtex8.mean)[colnames(gtex8.mean) =="Esophagus - Gastroesophageal Junction" ] <-  "Esophagus_Gastroesophageal_Junction"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Esophagus - Mucosa" ] <- "Esophagus_Mucosa"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Esophagus - Muscularis"] <-  "Esophagus_Muscularis"
colnames(gtex8.mean)[colnames(gtex8.mean) =="Heart - Atrial Appendage" ] <- "Heart_Atrial_Appendage"
colnames(gtex8.mean)[colnames(gtex8.mean) =="Heart - Left Ventricle" ] <- "Heart_Left_Ventricle"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Cells - Cultured fibroblasts"] <- "Cells_Cultured_fibroblasts"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Adrenal Gland"] <- "Adrenal_Gland"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Cervix Uteri"] <- "Cervix_Uteri"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Fallopian Tube"] <- "Fallopian_Tube"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Salivary Gland"] <- "Salivary_Gland"
colnames(gtex8.mean)[colnames(gtex8.mean) == "Small Intestine"] <- "Small_Intestine"


#Calculate nonparametric-expression specificity score
gtex8_spec <- data.frame(t(apply(gtex8.mean[gtex8.mean$gene_biotype == "protein_coding", !colnames(gtex8.mean) %in% c("gene_id","gene_name","chromosome_name","start_position","end_position", "strand", "gene_biotype")], 1,
function(row) row/sqrt(sum(row^2)))), gtex8.mean[gtex8.mean$gene_biotype == "protein_coding", colnames(gtex8.mean) %in% c("gene_id","gene_name","chromosome_name","start_position","end_position", "strand", "gene_biotype")], check.names=FALSE)
gtex8_specificity_percentile <- data.frame(apply(gtex8_spec[ , !colnames(gtex8_spec) %in% c("gene_id","gene_name","chromosome_name","start_position","end_position","strand","gene_biotype")], 2,
function(col) rank(-col)/length(col)), gtex8_spec[ , colnames(gtex8_spec) %in% c("gene_id","gene_name","chromosome_name","start_position","end_position","strand","gene_biotype")], check.names=FALSE)

# Write the expression and the specificity tables
write.csv(gtex8.mean, file = "../data/GTExV8_expression.csv", row.names = FALSE, quote=FALSE)
write.csv(gtex8_specificity_percentile, file = "../data/GTExV8_specificity.csv", row.names = FALSE, quote=FALSE)
