#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())



############# To edit #############

ANALYSIS <- "ASA24"
module <- "WeightMetaMerge"
###################################

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
metaFile <- args[2]
HEI.fileName <- args[3]

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
