#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: BMI analysis

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("BMI_Results")


##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)
# args <- "~/git/Diet_EWL_BariatricSurgery_2022"
args <- c(args, "weight_update_BLonly_excluded.txt")

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  message("gitRoot = ", gitRoot)
  
  today <- as.character(format(Sys.Date(), "%Y%b%d"))
  root <- paste0("~/BioLockJ_pipelines/")
  dir.create(root, showWarnings = FALSE)
  root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
  dir.create(root, showWarnings = FALSE)
  rootInput <- paste0(root, "input/")
  dir.create(rootInput, showWarnings = FALSE)
  
  gitInput <- file.path(gitRoot, "analysis", "data", "metadataTables")
  message("gitInput = ", gitInput)
  
  # file.copy(gitInput,
  #           rootInput,
  #           recursive = TRUE)
  
  dir.create(paste0(root, module, "/"), showWarnings = FALSE)
  message(paste0(root, module, "/"))
  
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")
  message("gitScripts = ", gitScripts)
  
  dir.create(paste0(root, module, "/script/"), showWarnings = FALSE)
  script = paste0(gitScripts,"/",str_subset(dir(gitScripts), module))
  file.copy(script,
            paste0(root, module, "/script/"),
            recursive = TRUE)
  
  dir.create(paste0(root, module, "/output/"), showWarnings = FALSE)
  
  dir.create(paste0(root, module, "/resources/"), showWarnings = FALSE)
  script = paste0(gitScripts,"/functions.R"); script
  file.copy(script,
            paste0(root, module, "/resources/"),
            recursive = TRUE)
  
  setwd(paste0(root, module, "/script/"))
  
}
rm(ANALYSIS, gitInput, gitRoot, gitScripts, root, rootInput, script, today)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

##### Set up functions file #####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

##### Set up input #####
inputDir = paste0(pipeRoot, "/input/metadataTables/")
inputFile <- args[2]

##### Set up output #####
outputDir = file.path(moduleDir,"output/")
RYGB.bmi.lm.file <- "RYGB_BMI_LM_changes_over_time.tsv"
SG.bmi.lm.file <- "SG_BMI_LM_changes_over_time.tsv"
surgerytype.bmi.wilcox.file <- "SurgeryType_BMI_wilcox_at_each_timepoint.tsv"
surgerytype.bmi.loss.wilcox.file <- "SurgeryType_BMI_Loss_from_BL_wilcox_at_each_timepoint.tsv"
Site.bmi.wilcox.file <- "Site_BMI_wilcox_at_each_timepoint.tsv"
updated.weightFile <- "updated_weight.tsv"

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
myTable$time = as.numeric(gsub("BIO-.-...-","",myTable$SampleID))
myTable$Surgery <- myTable$TypeofSurgery
myTable$Surgery <- ifelse(myTable$Surgery == "Gastric Bypass", "RYGB",
                          ifelse(myTable$Surgery == "Sleeve Gastrectomy", "SG",
                                 NA))
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))
row.names(myTable) = myTable$SampleID
# BL.df <- myTable[myTable$Timepoint == "BL",]
myTable <- myTable[!is.na(myTable$Surgery),]


##### RYGB patient BMI over time - Mixed Linear Model #####
Filt <- myTable[myTable$Surgery == "RYGB",]
Filt <- na.omit(Filt)

mlm <- lme(BMI_kgm2 ~ Timepoint, method = 'REML', random = ~1 | PatientID, data = Filt)
smry <- summary(mlm)
lm.df <- data.frame(smry$tTable)
lm.df$Comparison <- row.names(lm.df)
lm.df$Comparison <- gsub(pattern = "Timepoint", replacement = "BL v ", lm.df$Comparison)
lm.df$Metric <- "BMI"
lm.df$Adj_pvalue <- p.adjust(lm.df$p.value, method = "BH")
lm.df$p_value <- roundP(lm.df$Adj_pvalue)

file.path <- paste0(outputDir,RYGB.bmi.lm.file)
write.table(lm.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### SG patient BMI over time - Mixed Linear Model #####
Filt <- myTable[myTable$Surgery == "SG",]
Filt <- na.omit(Filt)

mlm <- lme(BMI_kgm2 ~ Timepoint, method = 'REML', random = ~1 | PatientID, data = Filt)
smry <- summary(mlm)
lm.df <- data.frame(smry$tTable)
lm.df$Comparison <- row.names(lm.df)
lm.df$Comparison <- gsub(pattern = "Timepoint", replacement = "BL v ", lm.df$Comparison)
lm.df$Metric <- "BMI"
lm.df$Adj_pvalue <- p.adjust(lm.df$p.value, method = "BH")
lm.df$p_value <- roundP(lm.df$Adj_pvalue)


file.path <- paste0(outputDir,SG.bmi.lm.file)
write.table(lm.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

myTable$Timepoint <- factor(myTable$Timepoint)

model <- aov(BMI_kgm2 ~ Surgery * Timepoint + Error(PatientID / Timepoint), data = myTable)
##### BMI differences in surgery type at each timepoint - wilcox #####
months <- c(0, 1, 6, 12, 18, 24)

Timepoint <- vector()
Avg_BMI_RYGB <- vector()
SD_BMI_RYGB <- vector()
Avg_BMI_SG <- vector()
SD_BMI_SG <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_BMI_RYGB[index] <- mean(Filt$BMI_kgm2[Filt$Surgery=="RYGB"], na.rm = TRUE)
  SD_BMI_RYGB[index] <- sd(Filt$BMI_kgm2[Filt$Surgery=="RYGB"], na.rm = TRUE)
  
  Avg_BMI_SG[index] <- mean(Filt$BMI_kgm2[Filt$Surgery=="SG"], na.rm = TRUE)
  SD_BMI_SG[index] <- sd(Filt$BMI_kgm2[Filt$Surgery=="SG"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$BMI_kgm2 ~ Filt$Surgery)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_BMI_RYGB, SD_BMI_RYGB, Avg_BMI_SG, SD_BMI_SG, p_value)

result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

m <- c(0, 1, 6, 12, 18, 24)
bio.conc <- vector()
index <- 1
for (i in 1:nrow(result.df)) {
  
  RYGB.mean <- paste0("(", round(result.df$Avg_BMI_RYGB[i], 2), " +/- ", round(result.df$SD_BMI_RYGB[i], 2), " kg/m2)")
  SG.mean <- paste0("(", round(result.df$Avg_BMI_SG[i], 2), " +/- ", round(result.df$SD_BMI_SG[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI in RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI in RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}

file.path <- paste0(outputDir,surgerytype.bmi.wilcox.file)
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### BMI loss over time between surgery types - wilcox #####
# Find baseline rows
bl <- which(myTable$Timepoint == "BL")

# Extract baseline BMI
blb <- myTable[bl, "BMI_kgm2"]
names(blb) = myTable[bl, "PatientID"]

# Add info to each row
myTable$Baseline_BMI = blb[myTable$PatientID]

# Calculate average loss from baseline
myTable$Loss_from_BL_BMI <- myTable$Baseline_BMI - myTable$BMI_kgm2

months <- c(1, 6, 12, 18, 24)
Timepoint <- vector()
Avg_Loss_RYGB <- vector()
SD_Loss_RYGB <- vector()
Avg_Loss_SG <- vector()
SD_Loss_SG <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_Loss_RYGB[index] <- mean(Filt$Loss_from_BL_BMI[Filt$Surgery=="RYGB"], na.rm = TRUE)
  SD_Loss_RYGB[index] <- sd(Filt$Loss_from_BL_BMI[Filt$Surgery=="RYGB"], na.rm = TRUE)
  
  Avg_Loss_SG[index] <- mean(Filt$Loss_from_BL_BMI[Filt$Surgery=="SG"], na.rm = TRUE)
  SD_Loss_SG[index] <- sd(Filt$Loss_from_BL_BMI[Filt$Surgery=="SG"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$Loss_from_BL_BMI ~ Filt$Surgery)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_Loss_RYGB, SD_Loss_RYGB, Avg_Loss_SG, SD_Loss_SG, p_value)
result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

bio.conc <- vector()
index <- 1

for (i in 1:nrow(result.df)) {
  
  RYGB.mean <- paste0("(", round(result.df$Avg_Loss_RYGB[i], 2), " +/- ", round(result.df$SD_Loss_RYGB[i], 2), " kg/m2)")
  SG.mean <- paste0("(", round(result.df$Avg_Loss_SG[i], 2), " +/- ", round(result.df$SD_Loss_SG[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI loss between RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI loss between RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}

file.path <- paste0(outputDir,surgerytype.bmi.loss.wilcox.file)
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### BMI differences in location at each timepoint - wilccox #####
months <- c(0, 1, 6, 12, 18, 24)

Timepoint <- vector()
Avg_Fargo <- vector()
SD_Fargo <- vector()
Avg_Cleveland <- vector()
SD_Cleveland <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_Fargo[index] <- mean(Filt$BMI_kgm2[Filt$Site=="Fargo"], na.rm = TRUE)
  SD_Fargo[index] <- sd(Filt$BMI_kgm2[Filt$Site=="Fargo"], na.rm = TRUE)
  
  Avg_Cleveland[index] <- mean(Filt$BMI_kgm2[Filt$Site=="Cleveland"], na.rm = TRUE)
  SD_Cleveland[index] <- sd(Filt$BMI_kgm2[Filt$Site=="Cleveland"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$BMI_kgm2 ~ Filt$Site)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_Fargo, SD_Fargo, Avg_Cleveland, SD_Cleveland, p_value)
result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

m <- c(0, 1, 6, 12, 18, 24)
bio.conc <- vector()
index <- 1
for (i in 1:nrow(result.df)) {
  
  Fargo.mean <- paste0("(", round(result.df$Avg_Fargo[i], 2), " +/- ", round(result.df$SD_Fargo[i], 2), " kg/m2)")
  Cleveland.mean <- paste0("(", round(result.df$Avg_Cleveland[i], 2), " +/- ", round(result.df$SD_Cleveland[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI between Fargo ", Fargo.mean, " and Cleveland patients ", Cleveland.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI between Fargo ", Fargo.mean, " and Cleveland patients ", Cleveland.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}

file.path <- paste0(outputDir,Site.bmi.wilcox.file)
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, updated.weightFile)
write.table(myTable, file.path,sep="\t",quote = FALSE, row.names = FALSE)
