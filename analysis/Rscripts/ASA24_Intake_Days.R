#Author: Alicia Sorgen
#Date: 2022 September 23
#Description: ASA24 intake data

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("ASA24_Intake_Days")
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
  inputFile <- "TNS_Master_file_enrolled_5-31-22.txt"
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
inputDir = paste0(pipeRoot,"/input/metadataTables/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.delim(file.path, sep="\t",header = TRUE)

# Assign weekdays and weekends
myTable$DayClass <- ifelse(myTable$DayoftheWeek == "Saturday" | myTable$DayoftheWeek == "Sunday", "Weekend",
                           "Weekday")
# Add filler 0s to timepoint
Timepoint <- myTable$Timepoint
Timepoint <- str_pad(Timepoint, 2, pad = "0")
myTable$SampleID <- paste0(myTable$ParticipantID., "-", Timepoint)

SampleIDs <- unique(myTable$SampleID)

timepoint <- vector()
numWeekday <- vector()
numWeekend <- vector()
CriteriaMet <- vector()
index <- 1

for (i in SampleIDs) {
  
  df <- myTable[ myTable$SampleID == i, ]
  
  timepoint[index] <- df$Timepoint[1]
  numWeekday[index] <- length(which(df$DayClass == "Weekday"))
  numWeekend[index] <- length(which(df$DayClass == "Weekend"))
  CriteriaMet[index] <- ifelse(length(which(df$DayClass == "Weekday")) > 0 & length(which(df$DayClass == "Weekend")) > 0, "Met",
                               "Not met")
  index <- index + 1
} # for (i in SampleIDs)

dFrame <- data.frame(timepoint, numWeekday, numWeekend, CriteriaMet)
file.path <- paste0(outputDir, "recall_days.tsv")
write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)

sumTable <- as.data.frame(table(dFrame$CriteriaMet, dFrame$timepoint))
names(sumTable)[names(sumTable) == "Var2"] <- "Timepoint"
df <- spread(sumTable, Var1, Freq)
df$Total <- df$Met + df$`Not met`
df$Met_Percent <- ( df$Met / df$ Total ) * 100

file.path <- paste0(outputDir, "recall_days_summary.tsv")
write.table(df, file.path, sep="\t",quote = FALSE, row.names = FALSE)
