#Author: Alicia Sorgen
#Date: 10-19-21
#Description: #Description: Average multiple data points from ASA24 data into single SampleID


## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(data.table); message("data.table:", packageVersion("data.table"))

rm(list=ls())

############# To edit #############
ANALYSIS <- "ASA24"
module <- paste0("ASA24Average")
date <- "2022Sep06"
###################################

##### Environment set up #####
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
  metaVariables <- args[2:length(args)]
} else {
  setwd(paste0(root, "/script"))
  
}
rm(date, today, root, ANALYSIS, module)

##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = file.path(pipeRoot, "input/metadataTables/")
moduleDir <- dirname(getwd())
output = file.path(pipeRoot, "input/metadataTables/")



## Read in ASA24 table
args <- commandArgs(trailingOnly = TRUE)
asa24Name <- args[1]
message(asa24Name)
asa24 <- read.delim(paste0(inputDir, asa24Name), sep="\t",header = TRUE)



## Add filler 0s to timepoint
Timepoint <- asa24$Timepoint
Timepoint <- str_pad(Timepoint, 2, pad = "0")
asa24$SampleID <- paste0(asa24$ParticipantID., "-", Timepoint)




## ID unique data points
ParticipantIDs <- unique(asa24$ParticipantID.)
SampleIDs <- unique(asa24$SampleID)

# Unique ParticipantIDs
message("Unique ParticipantIDs: ", length(ParticipantIDs))

# Unique SampleIDs
message("Unique SampleIDs: ", length(SampleIDs))

dFrame <- data.frame()

for (SampleID in SampleIDs) {
  
  asa24_SampleID <- asa24[ asa24$SampleID == SampleID, ]
  ParticipantID <- unique(asa24_SampleID$ParticipantID.)
  
  ## Assign column numbers for first to last ASA24 variable
  metaStart <- which(colnames(asa24_SampleID) == "KCAL")
  metaEnd <- which(colnames(asa24_SampleID) == "SampleID")-1
  
  
  ## Maintain only applicable columns
  avgCols <- asa24_SampleID[,metaStart:metaEnd]
  
  # Calculate column averages
  AVG <- as.data.frame(colMeans(avgCols))
  AVG <- t(AVG)
  rownames(AVG) <- SampleID
  # AVG$ParticipantID <- ParticipantID
  
  dFrame <- rbind(dFrame, AVG)
}

SampleID <- rownames(dFrame)

dFrame <- cbind(SampleID, dFrame)

outputFile <- paste0(output, "ASA24_metadata.tsv")
write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
