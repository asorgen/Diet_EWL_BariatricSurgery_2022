#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Summarize patient demographics

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("PatientCharacteristics")

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
  
  file.copy(gitInput,
            rootInput,
            recursive = TRUE)
  
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
# message("inputFile = ", inputFile)

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
