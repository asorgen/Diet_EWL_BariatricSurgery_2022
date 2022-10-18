#Author: Alicia Sorgen
#Date: 2022 June 15
#Description: Analysis of nutreint intake over time.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("Nutrient_Intake_Over_Time")
included <- c(  0
                , 1
                , 6
                , 12
                , 18
                , 24
)

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
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

resultsFileName <- "Nutrient_intake_differences_between_responders_24mo_patients.tsv"
summaryFileName <- "Average_nutrient_intake_by_responder_24mo_patients.tsv"
energyRatioFileName <- "Energy_ratio_differences_between_responders_24mo_patients.tsv"

##### Read in table and data prep #####
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
macros <- c("KCAL", "CARB", "PROT", "TFAT")

##### Set up energy ratios #####
# 4.184 kJ in 1 kcal
myTable$CARB_in_kcal <- (myTable$CARB * 16.7) / 4.184 # 16.7 kJ in 1 gram of carbs
myTable$PROT_in_kcal <- (myTable$PROT * 16.7) / 4.184 # 16.7 kJ in 1 gram of protein
myTable$TFAT_in_kcal <- (myTable$TFAT * 37.7) / 4.184 # 37.7 kJ in 1 gram of fat

myTable$CARB_energy_ratio <- (myTable$CARB_in_kcal / myTable$KCAL) * 100
myTable$PROT_energy_ratio <- (myTable$PROT_in_kcal / myTable$KCAL) * 100
myTable$TFAT_energy_ratio <- (myTable$TFAT_in_kcal / myTable$KCAL) * 100

##### Responder/Non-responder intake over time #####
results <- data.frame()
summary <- data.frame()

months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1

for (month in months) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  KCAL <- myTable$KCAL
  CARB <- myTable$CARB
  PROT <- myTable$PROT
  TFAT <- myTable$TFAT
  PEWL_kg <- myTable$PEWL_kg
  
  df <- data.frame(PatientID, Timepoint, PEWL_kg)
  df <- na.omit(df)
  df$Timepoint <- paste0("M", df$Timepoint)
  
  df <- spread(df, Timepoint, PEWL_kg)
  m <- which(colnames(df) == paste0("M", month))
  
  df$ResponderStatus <- ifelse(df[,m] >= 50, "Responder",
                               ifelse(df[,m] < 50, "Non-responder",
                                      NA))
  df <- df[!is.na(df$ResponderStatus),]
  bl <- which(colnames(df) == "M0")
  end <- which(colnames(df) == "ResponderStatus")-1
  df <- df %>%
    gather("Timepoint", "PEWL_kg", bl:end)
  df$Timepoint <- gsub("M", "", df$Timepoint)
  
  included <- c(included, month)
  
  df <- df[df$Timepoint %in% included,]
  df$Timepoint <- as.factor(df$Timepoint)
  df$Timepoint <- factor(df$Timepoint, levels = included)
  df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))
  
  names(KCAL) = myTable[, "SampleID"]
  names(CARB) = myTable[, "SampleID"]
  names(PROT) = myTable[, "SampleID"]
  names(TFAT) = myTable[, "SampleID"]
  
  # Add info to each row
  df$KCAL = KCAL[df$SampleID]
  df$CARB = CARB[df$SampleID]
  df$PROT = PROT[df$SampleID]
  df$TFAT = PROT[df$SampleID]
  
  for (macro in macros) {
    
    macroCol <- df[,which(colnames(df) == macro)]
    pID <- df$PatientID
    time <- df$Timepoint
    status <- df$ResponderStatus
    
    df2 <- data.frame(pID, time, status, macroCol)
    df2 <- na.omit(df2)
    
    for (aStatus in unique(status)) {
      
      df3 <- df2[df2$status == aStatus,]
      
      # Linear model
      fit <- anova(lme(macroCol ~ time, method = 'REML', random = ~1 | pID, data = df3))
      fit <- lme(macroCol ~ time, method = 'REML', random = ~1 | pID, data = df2)
      smry <- summary(fit)
      tTable <- as.data.frame(smry$tTable)
      tTable$time <- rownames(tTable)
      tTable$ResponderMonth <- month
      tTable$Variable <- macro
      tTable$Group <- aStatus
      rownames(tTable) <- NULL
      
      results <- rbind(results, tTable)
      
    } # for (aStatus in unique(status))
    
    
  } # for (macro in macros)
  
  
} # for (month in months)


file.path <- paste0(outputDir, "Outcome_intake_over_time.tsv")
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)




