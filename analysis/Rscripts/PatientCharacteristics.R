#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Summarize patient demographics

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("PatientCharacteristics")
date <- "2022Sep06"


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
today <- Sys.Date()
today <- format(today, "%Y%b%d")
today <- as.character(today)


if (date == today) {
  root <- paste0("~/BioLockJ_pipelines/",ANALYSIS,"_analysis_", today)
} else {
  root <- paste0("~/BioLockJ_pipelines/",ANALYSIS,"_analysis_", date)
}
root <- dir(root, pattern=module, full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  args <- commandArgs(trailingOnly = TRUE)
  inputFile <- args[1]
} else {
  setwd(paste0(root, "/script"))
  inputFile <- "weight_update.txt"
}
rm(date, today, root, ANALYSIS, module)

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

##### Set up output #####
outputDir = file.path(moduleDir,"output/")
sumTableName <- "SummaryTable.tsv"

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
myTable$time = as.numeric(gsub("BIO-.-...-","",myTable$SampleID))
myTable$Surgery <- myTable$TypeofSurgery
myTable$Surgery <- ifelse(myTable$Surgery == "Gastric Bypass", "RYGB",
                          ifelse(myTable$Surgery == "Sleeve Gastrectomy", "SG",
                                 NA))

row.names(myTable) = myTable$SampleID
BL.df <- myTable[myTable$Timepoint == "BL",]

##### Start summary #####
BL.df$Sex[is.na(BL.df$Sex)] <- "Not recorded"
var <- table(BL.df$Sex)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Sex"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- varTable

BL.df$Ethnicity[is.na(BL.df$Ethnicity)] <- "Not recorded"
var <- table(BL.df$Ethnicity)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Ethnicity"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

BL.df$TypeofSurgery[is.na(BL.df$TypeofSurgery)] <- "Not recorded"
var <- table(BL.df$TypeofSurgery)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "TypeofSurgery"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

BL.df$Site[is.na(BL.df$Site)] <- "Not recorded"
var <- table(BL.df$Site)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Site"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

ages <- BL.df$Age
Variable <- "Age"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

ages <- BL.df$Weight_kg
Variable <- "Weight (kg)"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

ages <- BL.df$BMI_kgm2
Variable <- "BMI (kg/m2)"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

BL.df$Ideal_BL_kg <- BL.df$Baseline_m^2 * 25
BL.df$Excess_kg <- BL.df$Weight_kg - BL.df$Ideal_BL_kg

wgt <- BL.df$Excess_kg
Variable <- "Excess weight (kg)"
Type <- "Avg & St Dev"
Count <- mean(wgt, na.rm = TRUE)
Percent <- sd(wgt, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

file.path <- paste0(outputDir, sumTableName)
write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
