# change working directory
setwd("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data")

# download necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)

if (!require(survival)){
  install.packages("survival")
}
library(survival)

if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

# read in clinical data csv
clinical <- read.csv("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/brca_clinical_data.csv") 

# prepare drug and radiation data
clin_query <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical", 
                       file.type = "xml")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")

# 1. age_at_initial_pathologic_diagnosis
# 2. continuous
# 3. clinical_rad anatomic_treatment_site; which site of the patient's tumor is treated 
# 4. categorical
# 5. i) There is some correlation between the age of the patient and their anatomic treatment site.
#   ii) The patients who are older are more likely to have lower survival rates in breast cancer. 
#  iii) Certain anatomic treatment sites increase patients' survival rates more significantly than others.

# merge the two data frames together into one
clinical_rad_merge <- merge(clinical, clinical_rad, by = "bcr_patient_barcode")

# create a boxplot with anatomic treatment site on the x-axis (categorical) and age on the y-axis (numerical)
# boxplot instead of scatter plot because one categorical and one numerical variable instead of two numerical variables
boxplot(formula = clinical_rad_merge$age_at_initial_pathologic_diagnosis ~ clinical_rad_merge$anatomic_treatment_site,
        data = clinical_rad_merge,
        xlab = "Anatomic Treatment Site",
        ylab = "Age at Initial Diagnosis",
        main = "Age at Initial Diagnosis vs. Anatomic Treatment Site")
# The boxplots all somewhat overlap with each other, indicating that there isn't any significant differences between certain anatomic treatment sites and the patient's age at initial diagnosis.  


# 1st variable: age!
age_clinical <- clinical_rad_merge

# no need to remove NA values in age column because there are none

# create column in age_cleaned_clinical where patients are labeled "Young", "Middle", or "Old"
young_mask <- ifelse(clinical_rad_merge$age_at_initial_pathologic_diagnosis <= 35, T, F)
middle_mask <- ifelse(clinical_rad_merge$age_at_initial_pathologic_diagnosis >= 35 & clinical_rad_merge$age_at_initial_pathologic_diagnosis <= 50, T, F)
old_mask <- ifelse(clinical_rad_merge$age_at_initial_pathologic_diagnosis >= 35, T, F)
age_clinical$age_at_initial_pathologic_diagnosis <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))

# make a survival time column 
age_clinical$survival_time <- ifelse(is.na(age_clinical$days_to_death),
                                     age_clinical$survival_time <- age_clinical$days_to_last_followup,
                                     age_clinical$survival_time <- age_clinical$days_to_death)

# remove NA values in survival time column
survival_na_mask_age <- ifelse(is.na(age_clinical$survival_time), F, T)
age_cleaned_clinical <- age_clinical[survival_na_mask_age, ]

# remove -Inf in survival time column
inf_mask_age <- ifelse(age_cleaned_clinical$survival_time == "-Inf", F, T)
age_cleaned_clinical <- age_cleaned_clinical[inf_mask_age, ]

# make a death event (T/F) column for survival plots
age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$death_event <- FALSE, age_cleaned_clinical$death_event <- TRUE)

# initialize survival object
survival_object_age <- Surv(time = age_cleaned_clinical$survival_time, event = age_cleaned_clinical$death_event)

# create fit object
age_fit <- survfit(survival_object_age ~ age_cleaned_clinical$age_at_initial_pathologic_diagnosis, data = age_cleaned_clinical)

# format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size = 20), 
                                                      axis.text = element_text(size = 16), 
                                                      legend.title = element_text(size = 14), 
                                                      legend.text = element_text(size = 12))

# save plot as png/jpeg on local computer 
jpeg("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/KM_plot_age.jpg")
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size = 20), 
                                                      axis.text = element_text(size = 16), 
                                                      legend.title = element_text(size = 14), 
                                                      legend.text = element_text(size = 12))
KM_plot_age
dev.off()

# 2nd variable: anatomic treatment site!
treatment_site_clinical <- clinical_rad_merge

# remove NA values from anatomic treatment site column
treatment_site_na_mask <- ifelse(is.na(treatment_site_clinical$anatomic_treatment_site), F, T)
treatment_cleaned_clinical <- treatment_site_clinical[treatment_site_na_mask, ]

# variable already categorical so don't need to make any further alterations to the column

# make a survival time column 
treatment_cleaned_clinical$survival_time <- ifelse(is.na(treatment_cleaned_clinical$days_to_death),
                                     treatment_cleaned_clinical$survival_time <- treatment_cleaned_clinical$days_to_last_followup,
                                     treatment_cleaned_clinical$survival_time <- treatment_cleaned_clinical$days_to_death)

# remove NA values in survival time column
survival_na_mask_treatment <- ifelse(is.na(treatment_cleaned_clinical$survival_time), F, T)
treatment_cleaned_clinical <- treatment_cleaned_clinical[survival_na_mask_treatment, ]

# remove -Inf in survival time column
inf_mask_treatment <- ifelse(treatment_cleaned_clinical$survival_time == "-Inf", F, T)
treatment_cleaned_clinical <- treatment_cleaned_clinical[inf_mask_treatment, ]

# make a death event (T/F) column for survival plots
treatment_cleaned_clinical$death_event <- ifelse(treatment_cleaned_clinical$vital_status == "Alive", treatment_cleaned_clinical$death_event <- FALSE, treatment_cleaned_clinical$death_event <- TRUE)

# initialize survival object
survival_object_treatment <- Surv(time = treatment_cleaned_clinical$survival_time, event = treatment_cleaned_clinical$death_event)

# create fit object
treatment_fit <- survfit(survival_object_treatment ~ treatment_cleaned_clinical$anatomic_treatment_site, data = treatment_cleaned_clinical)

# format and create KM plot
survplot_treatment <- ggsurvplot(treatment_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot_treatment <- survplot_treatment$plot + theme_bw() + theme(axis.title = element_text(size = 10), 
                                                      axis.text = element_text(size = 8), 
                                                      legend.title = element_text(size = 7), 
                                                      legend.text = element_text(size = 6))

# save plot as png/jpeg on local computer 
jpeg("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/analysis_data/KM_plot_treatment.jpg")
KM_plot_treatment <- survplot_treatment$plot + theme_bw() + theme(axis.title = element_text(size = 10), 
                                                      axis.text = element_text(size = 8), 
                                                      legend.title = element_text(size = 7), 
                                                      legend.text = element_text(size = 6))
KM_plot_treatment
dev.off()

# From the KM plot of anatomic treatment site, it is evident that certain anatomic treatment sites significantly decrease the survival probability of patients over time, such as distant recurrence and local recurrence. The distant treatment site seem to yield higher survival probability in patients over time. 
# From the KM plot of age, there isn't really that much of a distinction between young, middle-aged, and old patients in terms of survival probability. In fact, over time, it seems like the old patients have a slightly more stable and higher survival probability than both the younger and middle-aged patients, disproving my original hypothesis. 
# The p-value for the KM plot of age is 0.19 and the p-value for the KM plot of anatomic treatment site is less than 0.0001. Assuming a significance level of 0.05, we can conclude from these values that anatomic treatment sites have a significant impact on breast cancer survival (0.0001 < 0.05). 
# The differences in survival shown in the KM plot of anatomic treatment site do appear to be significant (~0.25-0.3 in survival probability for local recurrence, 1 for distant recurrence).
# The differences in survival shown in the KM plot of age appear to be less significant, as the patterns over time of the three age groups are very similar. 

