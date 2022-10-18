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
###################################

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)
# args <- "~/git/Diet_EWL_BariatricSurgery_2022"
args <- c(args, "TNS_Master_file_enrolled_5-31-22.txt")

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

##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = file.path(pipeRoot, "input/metadataTables/")
moduleDir <- dirname(getwd())
output = file.path(pipeRoot, "input/metadataTables/")



## Read in ASA24 table
asa24Name <- args[2]
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
  metaEnd <- which(colnames(asa24_SampleID) == "DayoftheWeek")-1
  
  
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
