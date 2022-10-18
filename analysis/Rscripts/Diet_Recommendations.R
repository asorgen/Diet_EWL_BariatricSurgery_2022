#Author: Alicia Sorgen
#Date: 2022 Oct 04
#Description: Analysis macronutrient energy ratios

rm(list=ls())

##### To edit #####
ANALYSIS <- "ASA24"

module <- paste0("Diet_Recommendations")

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
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]

##### Input file names #####
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

##### Output file names #####
outputDir = file.path(moduleDir,"output/")

##### Script variables #####
macros <- c("CARB", "PROT", "TFAT")
macro.labs <- c("Carbohydrate energy ratio (%)", "Protein energy ratio (%)", "Total fat energy ratio (%)")
xLabels <- c("Carbohydrates", "Protein", "Total fat")
Colors <- c("red", "green3")
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
macroColors <- c("cornflowerblue", "orange", "red2")
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

##### Set up energy ratios #####
# 4.184 kJ in 1 kcal
myTable$CARB_in_kcal <- (myTable$CARB * 16.7) / 4.184 # 16.7 kJ in 1 gram of carbs
myTable$PROT_in_kcal <- (myTable$PROT * 16.7) / 4.184 # 16.7 kJ in 1 gram of protein
myTable$TFAT_in_kcal <- (myTable$TFAT * 37.7) / 4.184 # 37.7 kJ in 1 gram of fat
# myTable$Macros_in_kcal <- myTable$CARB_in_kcal + myTable$PROT_in_kcal + myTable$TFAT_in_kcal

myTable$CARB_energy_ratio <- (myTable$CARB_in_kcal / myTable$KCAL) * 100
myTable$PROT_energy_ratio <- (myTable$PROT_in_kcal / myTable$KCAL) * 100
myTable$TFAT_energy_ratio <- (myTable$TFAT_in_kcal / myTable$KCAL) * 100

##### Set up recommended protein intake #####
# ~1.1-1.5 g / kg ideal body weight / day

myTable$PROT_1.1g <- 1.1 * myTable$Ideal_kg
myTable$PROT_1.5g <- 1.5 * myTable$Ideal_kg

##### Analyze actual v recommended protein intake at each timepoint #####
results <- data.frame()
summary <- data.frame()
plotList <- list()
index <- 1

months <- c("6", "12", "18", "24")
included <- c("0", "1")

for (month in months) {
  
  included <- c(included, month)
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  PROT <- myTable$PROT
  PROT_1.1g <- myTable$PROT_1.1g
  PROT_1.5g <- myTable$PROT_1.5g
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
  
  
  df <- df[df$Timepoint %in% included,]
  df$Timepoint <- as.factor(df$Timepoint)
  df$Timepoint <- factor(df$Timepoint, levels = included)
  df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))
  
  names(PROT) = myTable[, "SampleID"]
  names(PROT_1.1g) = myTable[, "SampleID"]
  names(PROT_1.5g) = myTable[, "SampleID"]

  # Add info to each row
  df$PROT = PROT[df$SampleID]
  df$PROT_1.1g = PROT_1.1g[df$SampleID]
  df$PROT_1.5g = PROT_1.5g[df$SampleID]
  df <- na.omit(df)
  
  df$PROT_needs_met <- ifelse(df$PROT > df$PROT_1.5g, "Met",
                              ifelse(df$PROT <= df$PROT_1.5g & df$PROT >= df$PROT_1.1g, "Met",
                                     "Deficient"))
  
  
  for (i in unique(df$Timepoint)) {
    
    df2 <- df[ df$Timepoint == i, ]
    
    freq <- as.data.frame(table(df2$ResponderStatus, df2$PROT_needs_met))
    freq2 <- as.data.frame(table(df2$ResponderStatus))
    
    n_NR <- freq2$Freq[which(freq2$Var1 == "Non-responder")]
    n_R <- freq2$Freq[which(freq2$Var1 == "Responder")]
    
    freq$Percent <- ifelse(freq$Var1 == "Non-responder", (freq$Freq / n_NR) * 100,
                           (freq$Freq / n_R) * 100)
    freq$ResponderMonth <- month
    freq$Timepoint <- i
    names(freq)[names(freq) == "Var1"] <- "ResponderStatus"
    names(freq)[names(freq) == "Var2"] <- "Protein"
    
    summary <- rbind(summary, freq)
    
    # Chi-squared test
    stats.summ <- chisq_test(df2$ResponderStatus, df2$PROT_needs_met)
    stats.summ$ResponderMonth <- month
    stats.summ$Timepoint <- i
    p_value <- roundP(stats.summ$p)
    
    results <- rbind(results, stats.summ)
    
    
    title.lab <- paste0("Diet at ", month.labs[which(included == i)])
    subtitle.lab <- paste0("Chi-squared ", p_value)
    x.lab <- paste0("Patients classified at ", month, "-months")
    
    plot <- ggplot(freq, aes(fill=Protein, y=Freq, x=ResponderStatus)) + 
      geom_bar(position="stack", stat="identity", color = "black")+
      labs(x = x.lab, y = "Number of patients", title = title.lab, subtitle = subtitle.lab)+
      scale_fill_manual(values=Colors); plot
    
    plot <- plot + geom_text(aes(label = paste0(round(Percent, 1), "%")), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "white"); plot
    
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (i in unique(df$Timepoint))
  
} # for (month in months)


resultsFileName <- "PROT_needs_met_by_ResponderStatus_chi_squared.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)


summaryFileName <- "PROT_needs_met_by_ResponderStatus_frequency.tsv"
file.path <- paste0(outputDir, summaryFileName)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)


plotFileName <- "PROT_needs_met_by_ResponderStatus.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in seq(1, length(plotList), 1)) {
  grid.arrange(plotList[[i]], 
               # plotList[[i+1]],
               # plotList[[i+2]],
               # plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 1, nrow = 1)
}
dev.off()

##### Diet Stacked Barplots #####
summary <- data.frame()
plotList <- list()
index <- 1

months <- c("6", "12", "18", "24")
included <- c("0", "1")

for (month in months) {
  
  included <- c(included, month)
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  CARB <- myTable$CARB_energy_ratio
  PROT <- myTable$PROT_energy_ratio
  TFAT <- myTable$TFAT_energy_ratio
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
  
  
  df <- df[df$Timepoint %in% included,]
  df$Timepoint <- as.factor(df$Timepoint)
  df$Timepoint <- factor(df$Timepoint, levels = included)
  df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))
  
  names(CARB) = myTable[, "SampleID"]
  names(PROT) = myTable[, "SampleID"]
  names(TFAT) = myTable[, "SampleID"]
  
  # Add info to each row
  df$CARB = CARB[df$SampleID]
  df$PROT = PROT[df$SampleID]
  df$TFAT = TFAT[df$SampleID]
  # df$Sum <- df$CARB + df$PROT + df$TFAT
  df <- na.omit(df)
  
  
  bl <- which(colnames(df) == "CARB")
  df <- df %>%
    gather("Diet", "Intake", bl:ncol(df))
  
  for (i in 1:length(included)) {
    
    df2 <- df[ df$Timepoint == included[i], ]
    
    averages <- df2 %>%
      group_by(ResponderStatus, Diet) %>%
      get_summary_stats(Intake, type = "mean")
    
    n_NR <- averages$n[which(averages$ResponderStatus == "Non-responder")][1]
    n_R <- averages$n[which(averages$ResponderStatus == "Responder")][1]
    
    averages$ResponderMonth <- month
    averages$Timepoint <- included[i]

    summary <- rbind(summary, averages)
    
    title.lab <- paste0("Diet at ", month.labs[which(included == included[i])])
    # subtitle.lab <- paste0("Chi-squared ", p_value)
    x.lab <- paste0("Patients classified at ", month, "-months")
    
    plot <- ggplot(averages, aes(fill=Diet, y=mean, x=ResponderStatus)) +
      # geom_bar(position="stack", stat="identity", color = "black")+
      geom_bar(position="fill", stat="identity", color = "black")+
      labs(x = x.lab, y = "Energy Ratio")+
      scale_fill_manual(values=macroColors); plot
    
    plot <- plot + labs(title = title.lab)
    # plot <- plot + labs(subtitle = subtitle.lab)
    
    # plot <- plot + geom_text(aes(label = paste0(round(mean, 1), "%")), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "white"); plot
    
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (i in 1:length(included))
  
  
} # for (month in months)


summaryFileName <- "Energy_Ratio_by_ResponderStatus_summary.tsv"
file.path <- paste0(outputDir, summaryFileName)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)


plotFileName <- "Energy_Ratio_by_ResponderStatus_stacked_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in seq(1, length(plotList), 1)) {
  grid.arrange(plotList[[i]],
               # plotList[[i+1]],
               # plotList[[i+2]],
               # plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 1, nrow = 1)
}
dev.off()

##### Analyze actual v recommended protein intake - post-op timepoints combined #####
results <- data.frame()
summary <- data.frame()
plotList <- list()
index <- 1

months <- c("12", "18", "24")
included <- c("0", "1", "6")

for (month in months) {
  
  included <- c(included, month)
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  PROT <- myTable$PROT
  PROT_1.1g <- myTable$PROT_1.1g
  PROT_1.5g <- myTable$PROT_1.5g
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
  
  
  df <- df[df$Timepoint %in% included,]
  df$Timepoint <- as.factor(df$Timepoint)
  df$Timepoint <- factor(df$Timepoint, levels = included)
  df$SampleID <- paste0(df$PatientID, "-", str_pad(df$Timepoint, width=2, pad="0"))
  
  names(PROT) = myTable[, "SampleID"]
  names(PROT_1.1g) = myTable[, "SampleID"]
  names(PROT_1.5g) = myTable[, "SampleID"]
  
  # Add info to each row
  df$PROT = PROT[df$SampleID]
  df$PROT_1.1g = PROT_1.1g[df$SampleID]
  df$PROT_1.5g = PROT_1.5g[df$SampleID]
  df <- na.omit(df)
  
  df$PROT_needs_met <- ifelse(df$PROT > df$PROT_1.5g, "Met",
                              ifelse(df$PROT <= df$PROT_1.5g & df$PROT >= df$PROT_1.1g, "Met",
                                     "Deficient"))
  
  
  # for (i in unique(df$Timepoint)) {
  
  df <- df[ !df$Timepoint == 0, ]
  
  freq <- as.data.frame(table(df$ResponderStatus, df$PROT_needs_met))
  freq2 <- as.data.frame(table(df$ResponderStatus))
  
  n_NR <- freq2$Freq[which(freq2$Var1 == "Non-responder")]
  n_R <- freq2$Freq[which(freq2$Var1 == "Responder")]
  
  freq$Percent <- ifelse(freq$Var1 == "Non-responder", (freq$Freq / n_NR) * 100,
                         (freq$Freq / n_R) * 100)
  freq$ResponderMonth <- month
  # freq$Timepoint <- i
  names(freq)[names(freq) == "Var1"] <- "ResponderStatus"
  names(freq)[names(freq) == "Var2"] <- "Protein"
  
  summary <- rbind(summary, freq)
  
  # Chi-squared test
  stats.summ <- chisq_test(df$ResponderStatus, df$PROT_needs_met)
  stats.summ$ResponderMonth <- month
  # stats.summ$Timepoint <- i
  p_value <- roundP(stats.summ$p)
  
  results <- rbind(results, stats.summ)
  
  
  title.lab <- paste0("All post-op timepoints")
  subtitle.lab <- paste0("Chi-squared ", p_value)
  x.lab <- paste0("Patients classified at ", month, "-months")
  
  plot <- ggplot(freq, aes(fill=Protein, y=Freq, x=ResponderStatus)) + 
    geom_bar(position="stack", stat="identity", color = "black")+
    labs(x = x.lab, y = "Number of patients")+
    scale_fill_manual(values=Colors); plot
  
  plot <- plot + labs(tag = tags[index])
  plot <- plot + labs(title = title.lab)
  plot <- plot + labs(subtitle = subtitle.lab)
  
  plot <- plot + geom_text(aes(label = paste0(round(Percent, 1), "%")), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "white"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  # } # for (i in unique(df$Timepoint))
  
} # for (month in months)


resultsFileName <- "PROT_needs_met_by_ResponderStatus_chi_squared_all_timepoints.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)


summaryFileName <- "PROT_needs_met_by_ResponderStatus_frequency_all_timepoints.tsv"
file.path <- paste0(outputDir, summaryFileName)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)


plotFileName <- "PROT_needs_met_by_ResponderStatus_all_timepoints.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in seq(1, length(plotList), 1)) {
  grid.arrange(plotList[[i]], 
               # plotList[[i+1]],
               # plotList[[i+2]],
               # plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 1, nrow = 1)
}
dev.off()

