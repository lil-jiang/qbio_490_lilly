# Change working directory
setwd("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data")

# Load neccessary packages
library(BiocManager)
library(TCGAbiolinks)
library(maftools)

# Read in clinical data csv
clinical <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_clinical_data.csv")

# Initialize maf_object with clinical annotations
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)

# merge maf_object@data and maf_object@clinical.data data frames
clinical_data_merge <- merge(maf_object@data, maf_object@clinical.data, by = "Tumor_Sample_Barcode")

# eliminate NA and indeterminate values in the breast carcinoma progesterone receptor status column
na_indeterminate_mask <- ifelse(is.na(clinical_data_merge$breast_carcinoma_progesterone_receptor_status) | clinical_data_merge$breast_carcinoma_progesterone_receptor_status == "Indeterminate", F, T)
clinical_data_merge <- clinical_data_merge[na_indeterminate_mask,]

# rewrite the breast carcinoma progesterone receptor status column as a factor
clinical_data_merge$breast_carcinoma_progesterone_receptor_status <- factor(clinical_data_merge$breast_carcinoma_progesterone_receptor_status)

# subset MAF object - positive breast carcinoma progesterone receptor status
positive_mask <- ifelse(clinical_data_merge$breast_carcinoma_progesterone_receptor_status == "Positive", T, F)
positive_patient_barcodes <- clinical_data_merge$Tumor_Sample_Barcode[positive_mask]
positive_maf <- subsetMaf(maf = maf_object,
                        tsb = positive_patient_barcodes)

# subset MAF object - negative breast carcinoma progesterone receptor status
negative_mask <- ifelse(clinical_data_merge$breast_carcinoma_progesterone_receptor_status == "Negative", T, F)
negative_patient_barcodes <- clinical_data_merge$Tumor_Sample_Barcode[negative_mask]
negative_maf <- subsetMaf(maf = maf_object,
                      tsb = negative_patient_barcodes)

# merge two subset MAFs to get top 10 most mutated genes for both 
m1.genes = getGeneSummary(positive_maf)[1:5]
m2.genes = getGeneSummary(negative_maf)[1:5]
mdt = merge(m1.genes[,.(Hugo_Symbol, MutatedSamples)], m2.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
mdt = mdt[order(max, decreasing = TRUE)]

# create a co-oncoplot with the top 10 most mutated genes for patients with positive receptor status vs.patients with negative receptor status 
coOncoplot(m1 = positive_maf, 
           m2 = negative_maf,
           m1Name = "Patients With + Breast Carcinoma Progesterone Receptor Status",
           m2Name = "Patients With - Breast Carcinoma Progesterone Receptor Status",
           borderCol = NA,
           genes = mdt[,Hugo_Symbol])
ggsave("Users/lilly/Desktop/QBIO490/qbio_490_lilly/week7_maf/PartIICooncoplot.png")

# From the cooncoplot, we can observe that there is a large discrepancy in the percentage of TP53 mutations between patients who have a positive receptor status and those who have a negative receptor status in regards to the breast carcinoma progesterone receptor.
# There are a significantly greater percentage of mutations in patients with the negative receptor status. 
# Because TP53 codes for a tumor suppressor protein, a higher percentage of mutations in the gene would implicate a faster growth of cancer. Similarly, patients who lack the breast carcinoma progesterone receptor are more prone to faster cancer cell growth, since treatment with hormone therapy is ineffective. 
# Therefore, it makes sense that the patients who have a negative receptor status have higher mutation numbers in the TP53 gene. 

# Create a contingency table with receptor status and TP53
TP53_mutations <- ifelse(clinical_data_merge$Hugo_Symbol == "TP53", "TP53 Mutation", "No TP53 Mutation")
TP53_mutations <- factor(TP53_mutations)
clinical$breast_carcinoma_estrogen_receptor_status <- factor(clinical$breast_carcinoma_progesterone_receptor_status, levels = c("Positive",  "Negative"))
contig <- table(clinical_data_merge$breast_carcinoma_progesterone_receptor_status, TP53_mutations)

# Run a Fisher's exact test
fisher_test <- fisher.test(contig)
fisher_test

# From the results of the Fisher's exact test, there is a 0.508x chance for a patient to have a negative breast carcinoma progesterone receptor status if the patient does not have any TP53 mutations compared to if they did. 
# The p-value is less than a significance level of 0.05, so we can conclude that the above statement is true. 

# Create a mosaic plot
mosaicplot(contig)
ggsave("Users/lilly/Desktop/QBIO490/qbio_490_lilly/week7_maf/PartIIMosaicPlot")

# Create a co-lollipop plot of TP53 divided between positive receptor status vs. negative receptor status
lollipopPlot2(m1 = positive_maf, 
              m2 = negative_maf, 
              m1_name = "Patients With + Breast Carcinoma Progesterone Receptor Status",
              m2_name = "Patients With - Breast Carcinoma Progesterone Receptor Status",
              gene = "TP53")
ggsave("Users/lilly/Desktop/QBIO490/qbio_490_lilly/week7_maf/PartIIColollipopPlot")

# Patients with a negative breast carcinoma progesterone receptor status seem to have mutations throughout the whole TP53 gene.
# Patients with a positive breast carcinoma progesterone receptor status seem to only have mutations that center primarily around the P53 domain. 

# Create a mafSurvival KM plot based on mutations in the TP53 gene
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == "Alive", T, F)

mafSurvival(maf = maf_object,
            genes = "TP53", 
            time = "days_to_last_followup", 
            Status = "Overall_Survival_Status", 
            isTCGA = TRUE)
ggsave("Users/lilly/Desktop/QBIO490/qbio_490_lilly/week7_maf/PartIIMAFSurvivalKMPlot")

# There is a slightly higher survival probability in patients who have TP53 mutations at around 1000-3000 days, however this doesn't make sense when reconsidering the function of the TP53 gene mentioned above (tumor-suppressing).
# Because the difference is by a very small amount and is only seen at around 1000-3000 days (start and end at around the same survival probability), we can conclude that the difference in survival probability between patients who have a normal TP53 gene compared to those who have mutated TP53 genes over time isn't significant 

