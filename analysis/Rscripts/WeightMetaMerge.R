#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())



############# To edit #############

ANALYSIS <- "microbiome"
module <- "WeightMetaMerge"
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
  metaFile <- args[1]
  HEI.fileName <- args[2]
} else {
  setwd(paste0(root, "/script"))
  metaFile <- "ASA24_metadata.tsv"
  HEI.fileName <- "HEI_data.txt"
}
rm(date, today, root, ANALYSIS, module)


##### Prep #####

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
output = file.path(moduleDir,"output/")

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)


##### Set input parameters #####
prevModule <- "ExcessWeightLoss"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
weightFile <- "updated_weight2.tsv"
weightTable <- read.delim(paste0(inputDir, weightFile), sep="\t",header = TRUE)


if (ANALYSIS == "microbiome") {
  
  mergeTable <- weightTable
  
} else {
  
  inputDir = paste0(pipeRoot, "/input/metadataTables/")
  metaTable <- read.delim(paste0(inputDir, metaFile), sep="\t",header = TRUE)
  
  metaIDs <- unique(metaTable$SampleID)
  weightIDs <- unique(weightTable$SampleID)
  combinedIDs <- c(metaIDs, weightIDs)
  
  duplicateIDs <- combinedIDs[duplicated(combinedIDs)]
  
  not_in_meta <- weightIDs[! weightIDs %in% duplicateIDs]
  not_in_weight <- metaIDs[! metaIDs %in% duplicateIDs]
  
  
  mergeTable <- merge(metaTable, weightTable, by = "SampleID", all = TRUE)
  row.names(mergeTable) = mergeTable$SampleID
  
  
  if (ANALYSIS == "ASA24") {
    ##### Edit HEI data for merging #####
    file.path <- paste0(inputDir, HEI.fileName)
    HEI.df <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
    HEI.df <- na.omit(HEI.df)
    
    HEI.df$Timepoint <- gsub(pattern = "0", replacement = "00", HEI.df$Timepoint)
    HEI.df$Timepoint <- gsub(pattern = "1$", replacement = "01", HEI.df$Timepoint)
    HEI.df$Timepoint <- gsub(pattern = "6", replacement = "06", HEI.df$Timepoint)
    
    SampleID <- paste0(HEI.df$ParticipantID_, "-", HEI.df$Timepoint)
    HEI <- HEI.df$HEI2015_TOTAL_SCORE
    
    df <- data.frame(SampleID, HEI)
    
    
    mergeTable <- merge(df, mergeTable, by = "SampleID", all = TRUE)
    
  }
  
  mergeTable <- mergeTable[!(mergeTable$SampleID %in% not_in_weight),]
  
}

write.table(mergeTable, paste0(output,"metadata.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
