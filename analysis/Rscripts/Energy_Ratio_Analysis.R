#Author: Alicia Sorgen
#Date: 2022 Oct 04
#Description: Analysis macronutrient energy ratios

rm(list=ls())

##### To edit #####
ANALYSIS <- "ASA24"

moduleRoot <- paste0("Energy_Ratio_Analysis")
params <- "~/git/Diet_EWL_BariatricSurgery_2022"

##### Global setup #####
R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(ggplot2); message("ggplot2:", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra:", packageVersion("gridExtra"))
library(scales); message("scales:", packageVersion("scales"))
library(rstatix); message("rstatix:", packageVersion("rstatix"))
library(ggpubr); message("ggpubr:", packageVersion("ggpubr"))
library(nlme); message("nlme:", packageVersion("nlme"))
library(dplyr); message("dplyr:", packageVersion("dplyr"))
library(data.table); message("data.table:", packageVersion("data.table"))

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
  rm(params)
}

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts"); rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }; rm(pipeline, ANALYSIS)
  
  if (any(dir(root) == "input") == FALSE) {
    # rootInput <- paste0(root, "input/")
    # dir.create(rootInput, showWarnings = FALSE)
    
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }; rm(gitInput)
  
  module <- moduleRoot
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }; rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R")
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  }; rm(scriptDir, moduleRoot)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files); rm(files, outputDir)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }; rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  
}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

##### Prep #####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]

roundP <- function(pvalue) {
  
  p <- ifelse(pvalue == 0, "p < 2e-16", 
              ifelse(pvalue < 0.0001, "p < 0.0001",
                     ifelse(pvalue >= 0.0001 & pvalue < 0.001, "p < 0.001",
                            ifelse(pvalue >= 0.001 & pvalue < 0.05, paste0("p = ", round(pvalue, digits = 3)),
                                   "NS"))))
  
  
  
}




##### Input file names #####
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

##### Output file names #####
outputDir = file.path(moduleDir,"output/")

##### Script variables #####
macros <- c("CARB", "PROT", "TFAT")
macro.labs <- c("Carbohydrate energy ratio (%)", "Protein energy ratio (%)", "Total fat energy ratio (%)")
tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.",
          "m.", "n.", "o.", "p.", "q.", "r.", "s", "t.", "u.", "v.", "w.", "x.", "y.", "z.")

##### Read in tables #####
# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Order Timepoint as factor
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))

# Assign SampleID to as dataframe row names
row.names(myTable) = myTable$SampleID

# Convert the timepoint data to numeric values.
myTable$time = as.numeric(myTable$time)

xLabels <- c("Energy (kcal)", "Carbohydrates (g)", "Protein (g)", "Total fat (g)")

##### Set up energy ratios #####
# 4.184 kJ in 1 kcal
myTable$CARB_in_kcal <- (myTable$CARB * 16.7) / 4.184 # 16.7 kJ in 1 gram of carbs
myTable$PROT_in_kcal <- (myTable$PROT * 16.7) / 4.184 # 16.7 kJ in 1 gram of protein
myTable$TFAT_in_kcal <- (myTable$TFAT * 37.7) / 4.184 # 37.7 kJ in 1 gram of fat

myTable$CARB_energy_ratio <- (myTable$CARB_in_kcal / myTable$KCAL) * 100
myTable$PROT_energy_ratio <- (myTable$PROT_in_kcal / myTable$KCAL) * 100
myTable$TFAT_energy_ratio <- (myTable$TFAT_in_kcal / myTable$KCAL) * 100

##### Macronutrient energy ratios between responders and non-responders #####
results <- data.frame()
months <- c("12", "18", "24")
included <- c("0", "1", "6")

for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    ratioCol <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]
    PEWL_kg <- myTable$PEWL_kg
    
    df <- data.frame(PatientID, Timepoint, PEWL_kg)
    df <- na.omit(df)
    df$Timepoint <- paste0("M", df$Timepoint)
    
    df <- spread(df, Timepoint, PEWL_kg)
    m <- which(colnames(df) == paste0("M", month))

    df$ResponderStatus <- ifelse(df[,m] >= 50, "Responder",
                                 ifelse(df[,m] < 50, "Suboptimal responder",
                                        NA))
    df <- df[!is.na(df$ResponderStatus),]
    bl <- which(colnames(df) == "M0")
    end <- which(colnames(df) == "ResponderStatus")-1
    df <- df %>%
      gather("Timepoint", "PEWL_kg", bl:end)
    df$Timepoint <- gsub("M", "", df$Timepoint)


    df <- df[df$Timepoint %in% included,]
    df$Timepoint <- as.factor(df$Timepoint)
    df$Timepoint <- factor(df$Timepoint, levels = included)
    df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))

    names(ratioCol) = myTable[, "SampleID"]

    # Add info to each row
    df$ratioCol = ratioCol[df$SampleID]
    df <- na.omit(df)


    # Averages
    averages <- df %>%
      group_by(Timepoint, ResponderStatus) %>%
      get_summary_stats(ratioCol, type = "mean_sd")
    averages$variable <- paste0(macro, " energy ratio")

    rRows <- which(averages$ResponderStatus == "Responder")
    nrRows <- which(averages$ResponderStatus == "Suboptimal responder")
    
    averages$final <- paste0(round(averages$mean, 2), " ± ", round(averages$sd, 2))
    
    # t-test
    stats.summ <- df %>%
      group_by(Timepoint) %>%
      t_test(ratioCol ~ ResponderStatus) %>%
      # adjust_pvalue(method = "BH") %>%
      add_significance()
    stats.summ$p_value <- roundP(stats.summ$p)
    stats.summ$ResponderMonth <- month
    stats.summ$Model <- "t-test"
    stats.summ <- stats.summ[,-(which(colnames(stats.summ) == "df"))]
    stats.summ$.y. <- paste0(macro, " energy ratio")
    stats.summ$Responder <- averages$final[rRows]
    stats.summ$Non.Responder <- averages$final[nrRows]
    results <- rbind(results, stats.summ)
    

    # Wilcoxon
    stats.summ <- df %>%
      group_by(Timepoint) %>%
      wilcox_test(ratioCol ~ ResponderStatus) %>%
      # adjust_pvalue(method = "BH") %>%
      add_significance()
    stats.summ$p_value <- roundP(stats.summ$p)
    stats.summ$ResponderMonth <- month
    stats.summ$Model <- "Wilcoxon"
    stats.summ$.y. <- paste0(macro, " energy ratio")
    stats.summ$Responder <- averages$final[rRows]
    stats.summ$Non.Responder <- averages$final[nrRows]
    results <- rbind(results, stats.summ)

    # index <- index + 1
    
  } # for (macro in macros)
} # for (month in months)


resultsFileName <- "Macro_Energy_Ratio_by_6M_12M_18M_24M_Outcome_wilcox_ttest.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Macronutrient energy ratios over time #####
results <- data.frame()
plotList <- list()
index <- 1

for (macro in macros) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  ratioCol <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]

  df <- data.frame(PatientID, Timepoint, ratioCol)
  df <- na.omit(df)
  
  df$Timepoint <- as.factor(df$Timepoint)
  df$Timepoint <- factor(df$Timepoint, levels = included)
  
  # Averages
  averages <- df %>%
    group_by(Timepoint) %>%
    get_summary_stats(ratioCol, type = "mean_sd")
  averages$variable <- paste0(macro, " energy ratio")

  # nrRows <- which(averages$ResponderStatus == "Suboptimal responder")

  averages$final <- paste0(round(averages$mean, 2), " ± ", round(averages$sd, 2))
  T0 <- averages$final[which(averages$Timepoint == 0)]
  T1 <- averages$final[which(averages$Timepoint == 1)]
  T6 <- averages$final[which(averages$Timepoint == 6)]
  T12 <- averages$final[which(averages$Timepoint == 12)]
  T18 <- averages$final[which(averages$Timepoint == 18)]
  T24 <- averages$final[which(averages$Timepoint == 24)]
  
  # summary <- rbind(summary, averages)
  
    # Tukey HSD
  stats.summ <- df %>%
    # group_by(Timepoint) %>%
    tukey_hsd(ratioCol ~ Timepoint) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  stats.summ$p_value <- roundP(stats.summ$p.adj)
  stats.summ$Macronutrient <- paste0(macro, " energy ratio")
  stats.summ$group1_Mean <- ifelse(stats.summ$group1 == 0, T0,
                                   ifelse(stats.summ$group1 == 1, T1,
                                          ifelse(stats.summ$group1 == 6, T6,
                                                 ifelse(stats.summ$group1 == 12, T12,
                                                        ifelse(stats.summ$group1 == 18, T18,
                                                               T24)))))
  stats.summ$group2_Mean <- ifelse(stats.summ$group2 == 0, T0,
                                   ifelse(stats.summ$group2 == 1, T1,
                                          ifelse(stats.summ$group2 == 6, T6,
                                                 ifelse(stats.summ$group2 == 12, T12,
                                                        ifelse(stats.summ$group2 == 18, T18,
                                                               T24)))))
  results <- rbind(results, stats.summ)
  
  stats.summ <- stats.summ %>%
    add_xy_position(x = "Timepoint", fun = "mean_se")
  
  plot <- ggbarplot(df, "Timepoint", "ratioCol",
                    # fill = "Status",
                    # color = "black", 
                    # palette = c("tomato", "steelblue"),
                    add = "mean_se",
                    label = FALSE, position = position_dodge()); plot
  
  x <- which(macros == macro)
  # tag <- tags[tagIndex]
  plot <- plot +
    labs(x = "Time (months)", y = macro.labs[x]); plot
  
  plot <- plot +
    stat_pvalue_manual(
      stats.summ,
      bracket.nudge.y = 0.5,
      # bracket.shorten = 1,
      bracket.size = 0.5,
      tip.length = 0.01,
      # remove.bracket = TRUE,
      size = 6,
      hide.ns = TRUE,
      step.increase = 0.01,
      label = "p.adj.signif"
    ); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
} # for (macro in macros)

resultsFileName <- "Macro_Energy_Ratio_over_time_tukey.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

plotFileName <- "energy_ratio_barplots_over_time_tukey.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], 
               # plotList[[i+1]],
               # plotList[[i+2]], plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 1, nrow = 1)
}
dev.off()

##### Macronutrient energy ratios over time (grouped by ResponderStatus) #####
results <- data.frame()
months <- c("12", "18", "24")
included <- c("0", "1", "6")
plotList <- list()
index <- 1
for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    ratioCol <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]
    PEWL_kg <- myTable$PEWL_kg
    
    df <- data.frame(PatientID, Timepoint, PEWL_kg)
    df <- na.omit(df)
    df$Timepoint <- paste0("M", df$Timepoint)
    
    df <- spread(df, Timepoint, PEWL_kg)
    m <- which(colnames(df) == paste0("M", month))
    
    df$ResponderStatus <- ifelse(df[,m] >= 50, "Responder",
                                 ifelse(df[,m] < 50, "Suboptimal responder",
                                        NA))
    df <- df[!is.na(df$ResponderStatus),]
    bl <- which(colnames(df) == "M0")
    end <- which(colnames(df) == "ResponderStatus")-1
    df <- df %>%
      gather("Timepoint", "PEWL_kg", bl:end)
    df$Timepoint <- gsub("M", "", df$Timepoint)
    
    
    df <- df[df$Timepoint %in% included,]
    df$Timepoint <- as.factor(df$Timepoint)
    df$Timepoint <- factor(df$Timepoint, levels = included)
    df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))
    
    names(ratioCol) = myTable[, "SampleID"]
    
    # Add info to each row
    df$ratioCol = ratioCol[df$SampleID]
    df <- na.omit(df)
    
    # Averages
    averages <- df %>%
      group_by(Timepoint, ResponderStatus) %>%
      get_summary_stats(ratioCol, type = "mean_sd")
    averages$variable <- paste0(macro, " energy ratio")
    
    averages$final <- paste0(round(averages$mean, 2), " ± ", round(averages$sd, 2))
    
    # Tukey HSD
    stats.summ <- df %>%
      group_by(ResponderStatus) %>%
      tukey_hsd(ratioCol ~ Timepoint) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    stats.summ$p_value <- roundP(stats.summ$p.adj)
    stats.summ$Macronutrient <- paste0(macro, " energy ratio")
    stats.summ$ResponderMonth <- month
    
    stats.summ$group1_Mean <- NA
    stats.summ$group2_Mean <- NA
    for (x in 1:length(included)) {
      
      Tx_R <- averages$final[which(averages$Timepoint == included[x] & averages$ResponderStatus == "Responder")]
      Tx_NR <- averages$final[which(averages$Timepoint == included[x] & averages$ResponderStatus == "Suboptimal responder")]
      
      stats.summ$group1_Mean[which(stats.summ$group1 == included[x] & stats.summ$ResponderStatus == "Responder")] <- Tx_R
      stats.summ$group1_Mean[which(stats.summ$group1 == included[x] & stats.summ$ResponderStatus == "Suboptimal responder")] <- Tx_NR
      
      stats.summ$group2_Mean[which(stats.summ$group2 == included[x] & stats.summ$ResponderStatus == "Responder")] <- Tx_R
      stats.summ$group2_Mean[which(stats.summ$group2 == included[x] & stats.summ$ResponderStatus == "Suboptimal responder")] <- Tx_NR
      
    } # for (x in 1:length(included))

    results <- rbind(results, stats.summ)
    
    stats.summ <- stats.summ %>%
      add_xy_position(x = "Timepoint", fun = "mean_se")
    
    plot <- ggbarplot(df, "Timepoint", "ratioCol",
                      fill = "ResponderStatus",
                      color = "black",
                      palette = c("tomato", "steelblue"),
                      add = "mean_se",
                      label = FALSE, position = position_dodge(),
                      facet.by = c("ResponderStatus")); plot
    
    x <- which(macros == macro)
    # tag <- tags[tagIndex]
    plot <- plot +
      labs(x = "Time (months)", y = macro.labs[x]); plot
    
    plot <- plot + labs(tag = tags[index]); plot
    
    plot <- plot + theme(legend.position="none"); plot
    
    plot <- plot +
      stat_pvalue_manual(
        stats.summ,
        bracket.nudge.y = 0.5,
        # bracket.shorten = 1,
        bracket.size = 0.5,
        tip.length = 0.01,
        # remove.bracket = TRUE,
        size = 6,
        hide.ns = TRUE,
        step.increase = 0.01,
        label = "p.adj.signif"
      ); plot
    
    plotList[[index]] <- plot
    index <- index + 1
    
    
  } # for (macro in macros)
} # for (month in months)

resultsFileName <- "Macro_Energy_Ratio_over_time_by_Outcome_tukey.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

plotFileName <- "energy_ratio_barplots_over_time_by_Outcome_tukey.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in seq(1, length(plotList), 3)) {
  grid.arrange(plotList[[i]], 
    plotList[[i+1]],
    plotList[[i+2]],
    # plotList[[i+3]],
    # plotList[[i+4]],
    ncol = 3, nrow = 1)
}
dev.off()
