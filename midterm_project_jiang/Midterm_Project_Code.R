setwd("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/midterm_project_jiang")
dir.create("outputs")
setwd("outputs")

library(SummarizedExperiment)
library(BiocManager)
library(TCGAbiolinks)
library(ggplot2)
library(maftools)
library(DESeq2)
library(EnhancedVolcano)

if (!require(survival)){
  install.packages("survival")
}
library(survival)

if(!require(survminer)){
  install.packages("survminer")
}

# read in necessary data 
clinical <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_clinical_data.csv")
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)
rna_counts <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_rna_count_data.csv")
rna_genes <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_rna_gene_data.csv")
rna_clinical <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_rna_clincial_data.csv")


# adjusting columns and changing row names because of a weird bug in the Intro to Transcriptomics HW 
rna_clinical <- rna_clinical[, -1]
rna_counts <- rna_counts[, -1]
rownames(rna_clinical) <- rna_clinical$barcode
rownames(rna_genes) <- rna_genes$gene_id
rownames(rna_counts) <- rownames(rna_genes)
colnames(rna_counts) <- rownames(rna_clinical)

# remove NA values and "Other" from surgical procedure column 
clinical <- clinical[ifelse(is.na(clinical$breast_carcinoma_surgical_procedure_name), F, T), ]
clinical <- clinical[ifelse(clinical$breast_carcinoma_surgical_procedure_name == "Other", F, T), ]

# make a KM survival plot for the type of surgical procedure
# make a survival time column 
clinical$survival_time <- ifelse(is.na(clinical$days_to_death),
                                 clinical$survival_time <- clinical$days_to_last_followup,
                                 clinical$survival_time <- clinical$days_to_death)

# remove NA and -Inf values from survival time column
survival_na_mask <- ifelse(is.na(clinical$survival_time), F, T)
clinical <- clinical[survival_na_mask, ]
inf_mask <- ifelse(clinical$survival_time == "-Inf", F, T)
clinical <- clinical[inf_mask, ]

# make a death event (T/F) column for survival plots
clinical$death_event <- ifelse(clinical$vital_status == "Alive", clinical$death_event <- FALSE, clinical$death_event <- TRUE)

# initialize survival object and create fit object
survival_object <- Surv(time = clinical$survival_time, event = clinical$death_event)
surgery_fit <- survfit(survival_object ~ clinical$breast_carcinoma_surgical_procedure_name, data = clinical)

# format and create KM plot
survplot <- ggsurvplot(surgery_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20),
                                              axis.text = element_text(size = 16),
                                              legend.title = element_text(size = 14),
                                              legend.text = element_text(size = 12))
KM_plot
jpeg("KM_plot_surgery.jpg")
dev.off()

# make a co-oncoplot comparing the mutations of patients who underwent lumpectomy vs. simple mastectomy
# merge maf_object@data and maf_object@clinical.data data frames
clinical_data_merge <- merge(maf_object@data, maf_object@clinical.data, by = "Tumor_Sample_Barcode")

# eliminate NA and "Other" values in the surgical procedure column
na_other_mask <- ifelse(is.na(clinical_data_merge$breast_carcinoma_surgical_procedure_name) | clinical_data_merge$breast_carcinoma_surgical_procedure_name == "Other", F, T)
clinical_data_merge <- clinical_data_merge[na_other_mask,]

# rewrite the surgical procedure column as a factor
clinical_data_merge$breast_carcinoma_surgical_procedure_name <- factor(clinical_data_merge$breast_carcinoma_surgical_procedure_name)

# subset MAF object - lumpectomy
lumpectomy_mask <- ifelse(clinical_data_merge$breast_carcinoma_surgical_procedure_name == "Lumpectomy", T, F)
lumpectomy_barcodes <- clinical_data_merge$Tumor_Sample_Barcode[lumpectomy_mask]
lumpectomy_maf <- subsetMaf(maf = maf_object,
                            
                            tsb = lumpectomy_barcodes)

# subset MAF object - simple mastectomy
simple_mastectomy_mask <- ifelse(clinical_data_merge$breast_carcinoma_surgical_procedure_name == "Simple Mastectomy", T, F)
simple_mastectomy_barcodes <- clinical_data_merge$Tumor_Sample_Barcode[simple_mastectomy_mask]
simple_mastectomy_maf <- subsetMaf(maf = maf_object,
                                   tsb = simple_mastectomy_barcodes)

# create oncoplot
coOncoplot(m1 = lumpectomy_maf, 
           m2 = simple_mastectomy_maf,
           m1Name = "Lumpectomy",
           m2Name = "Simple Mastectomy",
           borderCol = NA)
ggsave("cooncoplot.png")
dev.off()

# create a lollipop plot of the TP53 gene comparing lumpectomy patients with simple mastectomy patients
lollipopPlot2(m1 = lumpectomy_maf, 
              m2 = simple_mastectomy_maf, 
              m1_name = "Lumpectomy",
              m2_name = "Simple Mastectomy",
              gene = "TP53")
ggsave("lollipop_plot.png")
dev.off()

# make a contingency table and run a Fisher test to compare the age of the patient and their surgical procedure
clinical_surg_compare <- clinical[ifelse(clinical$breast_carcinoma_surgical_procedure_name == "Lumpectomy" | clinical$breast_carcinoma_surgical_procedure_name == "Simple Mastectomy", T, F), ]
clinical_surg_compare$breast_carcinoma_surgical_procedure_name <- factor(clinical_surg_compare$breast_carcinoma_surgical_procedure_name)
clinical_surg_compare$age_category <- ifelse(clinical_surg_compare$age_at_initial_pathologic_diagnosis <= 58, "young", "old")
clinical_surg_compare$race_list <- factor(clinical_surg_compare$age_category)
contig <- table(clinical_surg_compare$age_category, clinical_surg_compare$breast_carcinoma_surgical_procedure_name)
mosaicplot(contig)
ggsave("mosaic_plot.png")

# Run a Fisher's exact test
fisher_test <- fisher.test(contig)
fisher_test

# make a volcano plot to compare the counts of each gene to the age of the patient 
# convert the categorical variable and covariables into factor data type
rna_clinical$age_category <- factor(rna_clinical$age_category)
rna_clinical$race <- factor(rna_clinical$race)
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)

# check if there are existing NA values in any of the columns and eliminate them
sum(is.na(rna_clinical$age_category))
sum(is.na(rna_clinical$race))
sum(is.na(rna_clinical$ajcc_pathologic_stage))
na_mask <- ifelse(is.na(rna_clinical$ajcc_pathologic_stage), F, T)
rna_clinical <- rna_clinical[na_mask, ] 
rna_counts <- rna_counts[ ,na_mask]

# remove all genes where the total number of counts (across all patients) is less than 10
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums >= 10, T, F)
rna_counts <- rna_counts[low_counts_mask, ]
rna_genes <- rna_genes[low_counts_mask, ]

# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~race + ajcc_pathologic_stage + age_category)
dds_obj <- DESeq(dds) 

# see what comparisons got run
resultsNames(dds_obj)  
results <- results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

# rewrite results with the rows and columns of interest and rename the columns 
gene_name_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
results <- data.frame(rna_genes$gene_name[gene_name_mask], results@rownames, results$log2FoldChange, results$pvalue, results$padj, -log10(results$padj))
colnames(results) <- c("gene_name", "gene_id", "log2FoldChange", "pvalue", "padj", "-log10(padj)")

# select rows (genes) that have a padj value < 0.05
padj_mask <- ifelse(results$padj < 0.05, T, F)
sig_results <- results[padj_mask, ]

# sort the data frame by descending log2foldchange and include only the genes where the log2foldchange is > 1
up_reg_results <- sig_results[order(-sig_results$log2FoldChange), ]
up_reg_results <- up_reg_results[ifelse(up_reg_results$log2FoldChange > 1, T, F), ]

# sort the data frame by ascending log2foldchange and include only the genes where the log2foldchange is < -1
down_reg_results <- sig_results[order(sig_results$log2FoldChange), ]
down_reg_results <- down_reg_results[ifelse(down_reg_results$log2FoldChange < -1, T, F), ]

# create new row names that are more descriptive
.rowNamesDF(up_reg_results, make.names=TRUE) <- up_reg_results$gene_id
.rowNamesDF(down_reg_results, make.names=TRUE) <- down_reg_results$gene_id

# check to see if the TP53 gene is an upregulated gene, downregulated gene, or neither
"TP53" %in% up_reg_results$gene_name
"TP53" %in% down_reg_results$gene_name

# create volcano plot 
par(mar=c(1,1,1,1))
EnhancedVolcano(results,
                lab = results$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                legendLabels = c("Not sig.", "Log2FC", "pvalue", "pvalue & Log2FC"),
                legendPosition = "right")
ggsave("enhanced_volcano")
dev.off()

