#Author: Alicia Sorgen
#Date: 2022 September 07
#Description: Generate summary barplots.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
moduleRoot <- paste0("Barplot_summaries")
params <- "."


##### Libraries #####
library(stringr)
library(ggplot2)
library(gridExtra)
library(scales)
library(rstatix)
library(ggpubr)
library(nlme)

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

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Prep data table for analysis #####
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

# 4.184 kJ in 1 kcal
myTable$CARB_in_kcal <- (myTable$CARB * 16.7) / 4.184 # 16.7 kJ in 1 gram of carbs
myTable$PROT_in_kcal <- (myTable$PROT * 16.7) / 4.184 # 16.7 kJ in 1 gram of protein
myTable$TFAT_in_kcal <- (myTable$TFAT * 37.7) / 4.184 # 37.7 kJ in 1 gram of fat

myTable$CARB_energy_ratio <- (myTable$CARB_in_kcal / myTable$KCAL) * 100
myTable$PROT_energy_ratio <- (myTable$PROT_in_kcal / myTable$KCAL) * 100
myTable$TFAT_energy_ratio <- (myTable$TFAT_in_kcal / myTable$KCAL) * 100

myTable$MacroTotal <- myTable$PROT + myTable$CARB + myTable$TFAT


tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
macros <- c("KCAL", "CARB", "PROT", "TFAT", "MacroTotal")
macro.labs <- c("Calorie intake (kcal)", "Carbohydrate intake (g)", "Protein intake (g)", "Total fat intake (g)", "Total macronutrient intake (g)")

##### Bar plots intake between responder status groups at each timepoint #####
months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1
plotList <- list()
tagIndex <- 1
for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    MACRO <- myTable[,which(colnames(myTable) == macro)]

    # Create smaller data table with relevant columns
    df <- data.frame(PatientID, Timepoint, MACRO)
    
    # Pick out the patients designated at a particular responder status
    rIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month responder"))]
    nrIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month non-responder"))]
    
    # Filter table to only contain designated IDs
    df <- df[df$PatientID %in% c(rIDs, nrIDs),]
    
    # Filter table to only contain needed timepoints
    df <- df[df$Timepoint %in% included,]
    
    df$Status <- ifelse(df$PatientID %in% rIDs, "Responder",
                                "Suboptimal responder")
    
    stat.test <- df %>%
      group_by(Timepoint) %>%
      wilcox_test(MACRO ~ Status) %>%
      # adjust_pvalue(method = "BH") %>%
      add_significance()
    
    stat.test <- stat.test %>%
      add_xy_position(x = "Timepoint", fun = "mean_se")
    
    plot <- ggbarplot(df, "Timepoint", "MACRO",
              fill = "Status",
              color = "black", 
              palette = c("tomato", "steelblue"),
              add = "mean_se",
              label = FALSE, position = position_dodge())
    
    x <- which(macros == macro)
    tag <- tags[tagIndex]
    plot <- plot +
      labs(x = "Time (months)", y = macro.labs[x], tag = tag)
    
    plot <- plot +
      stat_pvalue_manual(
        stat.test,
        bracket.nudge.y = 0.5,
        # bracket.shorten = 1,
        bracket.size = 0.5,
        tip.length = 0.01,
        # remove.bracket = TRUE,
        size = 6,
        hide.ns = TRUE,
        # label = "p.adj.signif"
        label = "p.signif"
      )
    plot <- plot +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
      theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
    plot
    plotList[[index]] <- plot
    
    
    avg <- df %>%
      group_by(Timepoint, Status) %>%
      get_summary_stats(MACRO, type = "mean_sd")
    avg$variable <- macro
    
    index <- index + 1
    
  } # for (macro in macros)
    
  tagIndex <- tagIndex + 1
} # for (month in months)

plotFileName <- "response_intake_barplots_individual.pdf"
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

##### Bar plots intake between surgery types at each timepoint #####
index <- 1
plotList <- list()
final <- data.frame()
myTable$KCAL_energy_ratio <- 1
macros2 <- c("KCAL", "CARB", "PROT", "TFAT")

for (macro in macros2) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  MACRO <- myTable[,which(colnames(myTable) == macro)]
  Surgery <- myTable$Surgery
  Ratio <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]
  
  # Create smaller data table with relevant columns
  df <- data.frame(PatientID, Timepoint, MACRO, Surgery, Ratio)
  df <- na.omit(df)
  
  # Averages
  averages <- df %>%
    group_by(Timepoint, Surgery) %>%
    get_summary_stats(MACRO, type = "mean_sd")
  
  
  stat.test <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(MACRO ~ Surgery) %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stat.test$RYGB_Avg <- averages$mean[which(averages$Surgery == "RYGB")]
  stat.test$RYGB_SD <- averages$sd[which(averages$Surgery == "RYGB")]
  
  stat.test$SG_Avg <- averages$mean[which(averages$Surgery == "SG")]
  stat.test$SG_SD <- averages$sd[which(averages$Surgery == "SG")]
  
  averages_ratio <- df %>%
    group_by(Timepoint, Surgery) %>%
    get_summary_stats(Ratio, type = "mean_sd")
  
  stat.test$RYGB_Ratio_Avg <- averages_ratio$mean[which(averages$Surgery == "RYGB")]
  stat.test$RYGB_Ratio_SD <- averages_ratio$sd[which(averages$Surgery == "RYGB")]
  stat.test$SG_Ratio_Avg <- averages_ratio$mean[which(averages$Surgery == "SG")]
  stat.test$SG_Ratio_SD <- averages_ratio$sd[which(averages$Surgery == "SG")]
  
  stat.test$.y. <- macro
  
  if (macro == "KCAL") {
    stat.test$RYGB <- paste0(round(stat.test$RYGB_Avg, 2), " ± ", round(stat.test$RYGB_SD, 2))
    stat.test$SG <- paste0(round(stat.test$SG_Avg, 2), " ± ", round(stat.test$SG_SD, 2))
  } else {
    stat.test$RYGB <- paste0(round(stat.test$RYGB_Avg, 2), " ± ", round(stat.test$RYGB_SD, 2), " (", round(stat.test$RYGB_Ratio_Avg, 1), "%)")
    stat.test$SG <- paste0(round(stat.test$SG_Avg, 2), " ± ", round(stat.test$SG_SD, 2), " (", round(stat.test$SG_Ratio_Avg, 1), "%)")
  }
  
  stat.test$BH_adj <- p.adjust(stat.test$p, method = "BH")
  final <- rbind(final, stat.test)
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "Timepoint", fun = "mean_se")
  
  plot <- ggbarplot(df, "Timepoint", "MACRO",
                    fill = "Surgery",
                    color = "black",
                    palette = c("#3171BC", "#E8C241"),
                    add = "mean_se",
                    label = FALSE, position = position_dodge())

  x <- which(macros2 == macro)
  tag <- tags[index]

  plot <- plot +
    labs(x = "Time (months)", y = macro.labs[x],
         tag = tag)

  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.5,
      # bracket.shorten = 1,
      bracket.size = 0.5,
      tip.length = 0.01,
      # remove.bracket = TRUE,
      size = 6,
      hide.ns = TRUE,
      # label = "p.adj.signif"
      label = "p.signif"
    )
  plot <- plot +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
    theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
  plot
  plotList[[index]] <- plot
  
  index <- index + 1
  
} # for (macro in macros2)

plotFileName <- "surgery_intake_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 7.5, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 2, nrow = 2)
}
dev.off()

##### Make tables for surgery type and nutrients #####

for (macro in macros2) {
  
  df2 <- final[final$.y. == macro,]
  
  TP <- vector()
  RYGB <- vector()
  SG <- vector()
  pValue <- vector()
  pAdj <- vector()
  index <- 1
  
  TP[index] <- NA
  RYGB[index] <- paste0("(n = ", df2$n1[1], ")")
  SG[index] <- paste0("(n = ", df2$n2[1], ")")
  pValue[index] <- NA
  pAdj[index] <- NA
  
  for (i in 1:nrow(df2)) {
    index <- index + 1
    
    TP[index] <- df2$Timepoint[i]
    RYGB[index] <- df2$RYGB[i]
    SG[index] <- df2$SG[i]
    pValue[index] <- round(df2$p[i], 3)
    pAdj[index] <- round(df2$BH_adj[i], 3)
    
  } # for (i in 1:nrow(df2))
  
  temp <- data.frame(TP, RYGB, SG, pValue, pAdj)
  names(temp)[names(temp) == "RYGB"] <- paste0("RYGB_", macro)
  names(temp)[names(temp) == "SG"] <- paste0("SG_", macro)
  names(temp)[names(temp) == "pValue"] <- paste0("pValue_", macro)
  names(temp)[names(temp) == "pAdj"] <- paste0("pAdj_", macro)
  
  if (macro == macros2[1]) {
    resultsTable <- temp
  } else {
    temp <- temp[,-1]
    resultsTable <- cbind(resultsTable, temp)
  } # if (macro == macros2[1])
  
} # for (macro in macros2)

file.path <- paste0(outputDir, "Wilcoxon_Nutrient_Table_Surgery_Type.tsv")
write.table(resultsTable, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### Bar plots PEWL between responder status groups at each timepoint #####
months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1
plotList <- list()

for (month in months) {
  
  included <- c(included, month)
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  PEWL <- myTable$PEWL_kg
  
  # Create smaller data table with relevant columns
  df <- data.frame(PatientID, Timepoint, PEWL)
  df <- na.omit(df)
  
  # Pick out the patients designated at a particular responder status
  rIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month responder"))]
  nrIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month non-responder"))]
  
  # Filter table to only contain designated IDs
  df <- df[df$PatientID %in% c(rIDs, nrIDs),]
  
  # Filter table to only contain needed timepoints
  df <- df[df$Timepoint %in% included[2:length(included)],]
  
  df$Status <- ifelse(df$PatientID %in% rIDs, "Responder",
                      "Suboptimal responder")
  
  stat.test <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(PEWL ~ Status) %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stat.test <- stat.test %>%
    add_xy_position(x = "Timepoint", fun = "mean_se")
  
  plot <- ggbarplot(df, "Timepoint", "PEWL",
                    fill = "Status",
                    color = "black", 
                    palette = c("tomato", "steelblue"),
                    add = "mean_se",
                    label = FALSE, position = position_dodge())
  
  x <- which(macros == macro)
  
  plot <- plot +
    labs(x = "Time (months)", y = "Excess weight loss (%)")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.5,
      # bracket.shorten = 1,
      bracket.size = 0.5,
      tip.length = 0.01,
      # remove.bracket = TRUE,
      size = 6,
      hide.ns = TRUE,
      # label = "p.adj.signif"
      label = "p.signif"
    )
  plot <- plot +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
    theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
  plot
  plotList[[index]] <- plot
  
  index <- index + 1
  
} # for (month in months)

plotFileName <- "response_PEWL_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], 
               ncol = 3, nrow = 1)
}
dev.off()

##### Bar plots PEWL between surgery types at each timepoint #####
index <- 1
plotList <- list()

PatientID <- myTable$PatientID
Timepoint <- myTable$time
PEWL <- myTable$PEWL_kg
Surgery <- myTable$Surgery

# Create smaller data table with relevant columns
df <- data.frame(PatientID, Timepoint, PEWL, Surgery)
df <- na.omit(df)
df <- df[!(df$Timepoint == 0),]

averages$variable <- macro

stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(PEWL ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "PEWL",
                  fill = "Surgery",
                  color = "black", 
                  # palette = c("#3171BC", "#E8C241"),
                  add = "mean_se",
                  label = FALSE, position = position_dodge())

plot <- plot + scale_fill_grey()

x <- which(macros == macro)

plot <- plot +
  labs(x = "Time (months)", y = "Excess weight loss (%)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot <- plot +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
plot
plotList[[index]] <- plot

plotFileName <- "surgery_PEWL_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()


##### Bar plots intake between responder status groups at each timepoint (formatted pdf) #####
months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1
plotList <- list()
tagIndex <- 1
for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    # Create smaller data table with relevant columns
    df <- data.frame(PatientID, Timepoint, MACRO)
    
    # Pick out the patients designated at a particular responder status
    rIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month responder"))]
    nrIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month non-responder"))]
    
    # Filter table to only contain designated IDs
    df <- df[df$PatientID %in% c(rIDs, nrIDs),]
    
    # Filter table to only contain needed timepoints
    df <- df[df$Timepoint %in% included,]
    
    df$Status <- ifelse(df$PatientID %in% rIDs, "Responder",
                        "Suboptimal responder")
    
    stat.test <- df %>%
      group_by(Timepoint) %>%
      wilcox_test(MACRO ~ Status) %>%
      # adjust_pvalue(method = "BH") %>%
      add_significance()
    
    stat.test <- stat.test %>%
      add_xy_position(x = "Timepoint", fun = "mean_se")
    
    plot <- ggbarplot(df, "Timepoint", "MACRO",
                      fill = "Status",
                      color = "black", 
                      palette = c("tomato", "steelblue"),
                      add = "mean_se",
                      label = FALSE, position = position_dodge())
    
    x <- which(macros == macro)
    plot <- plot +
      labs(x = "Time (months)", y = macro.labs[x])
    
    # plot <- plot + labs(tag = tags[tagIndex])
    
    plot <- plot +
      stat_pvalue_manual(
        stat.test,
        bracket.nudge.y = 0.5,
        # bracket.shorten = 1,
        bracket.size = 0.5,
        tip.length = 0.01,
        # remove.bracket = TRUE,
        size = 6,
        hide.ns = TRUE,
        # label = "p.adj.signif"
        label = "p.signif"
      )
    plot <- plot +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
      theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
    plot
    plotList[[index]] <- plot
    
    
    avg <- df %>%
      group_by(Timepoint, Status) %>%
      get_summary_stats(MACRO, type = "mean_sd")
    avg$variable <- macro
    
    index <- index + 1
    
  } # for (macro in macros)
  
  tagIndex <- tagIndex + 1
} # for (month in months)

plotFileName <- "response_intake_barplots_macros.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 5)) {
  grid.arrange(plotList[[i]], 
               plotList[[i+1]],
               plotList[[i+2]],
               plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 2, nrow = 2)
}
dev.off()

plotFileName <- "response_intake_barplots_macros_only.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in seq(1, length(plotList), 5)) {
  grid.arrange(#plotList[[i]], 
    plotList[[i+1]],
    plotList[[i+2]],
    plotList[[i+3]],
    # plotList[[i+4]],
    ncol = 3, nrow = 1)
}
dev.off()

##### Macronutrients over time (grouped by ResponderStatus) #####
results <- data.frame()
months <- c("6", "12", "18", "24")
included <- c("0", "1")
plotList <- list()
index <- 1
for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    macroCol <- myTable[,which(colnames(myTable) == macro)]
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
    
    names(macroCol) = myTable[, "SampleID"]
    
    # Add info to each row
    df$macroCol = macroCol[df$SampleID]
    df <- na.omit(df)
    
    # Averages
    averages <- df %>%
      group_by(Timepoint, ResponderStatus) %>%
      get_summary_stats(macroCol, type = "mean_sd")
    averages$variable <- paste0(macro, " energy ratio")
    
    averages$final <- paste0(round(averages$mean, 2), " ± ", round(averages$sd, 2))
    
    # Tukey HSD
    stats.summ <- df %>%
      group_by(ResponderStatus) %>%
      tukey_hsd(macroCol ~ Timepoint) %>%
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
    
    plot <- ggbarplot(df, "Timepoint", "macroCol",
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

resultsFileName <- "macronutrient_over_time_by_Outcome_tukey.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

plotFileName <- "macronutrient_barplots_over_time_by_Outcome_tukey.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 5)) {
  grid.arrange(plotList[[i]], 
               plotList[[i+1]],
               plotList[[i+2]],
               plotList[[i+3]],
               plotList[[i+4]],
               ncol = 2, nrow = 3)
}
dev.off()


plotFileName <- "macronutrient_barplots_over_time_by_Outcome_tukey_no_total.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 5)) {
  grid.arrange(plotList[[i]], 
               plotList[[i+1]],
               plotList[[i+2]],
               plotList[[i+3]],
               # plotList[[i+4]],
               ncol = 2, nrow = 2)
}
dev.off()

plotFileName <- "macronutrient_barplots_over_time_by_Outcome_tukey_macros_only.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in seq(1, length(plotList), 5)) {
  grid.arrange(#plotList[[i]], 
    plotList[[i+1]],
    plotList[[i+2]],
    plotList[[i+3]],
    # plotList[[i+4]],
    ncol = 3, nrow = 1)
}
dev.off()
##### Bar plots energy ratios between responder status groups at each timepoint #####
months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1
plotList <- list()
macros <- c("CARB", "PROT", "TFAT")
macro.labs <- c("Energy derived from carbohydrates (%)", "Energy derived from protein (%)", "Energy derived from fat (%)")

for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    MACRO <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]
    
    # Create smaller data table with relevant columns
    df <- data.frame(PatientID, Timepoint, MACRO)
    
    # Pick out the patients designated at a particular responder status
    rIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month responder"))]
    nrIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month non-responder"))]
    
    # Filter table to only contain designated IDs
    df <- df[df$PatientID %in% c(rIDs, nrIDs),]
    
    # Filter table to only contain needed timepoints
    df <- df[df$Timepoint %in% included,]
    
    df$Status <- ifelse(df$PatientID %in% rIDs, "Responder",
                        "Suboptimal responder")
    
    stat.test <- df %>%
      group_by(Timepoint) %>%
      wilcox_test(MACRO ~ Status) %>%
      # adjust_pvalue(method = "BH") %>%
      add_significance()
    
    stat.test <- stat.test %>%
      add_xy_position(x = "Timepoint", fun = "mean_se")
    
    plot <- ggbarplot(df, "Timepoint", "MACRO",
                      fill = "Status",
                      color = "black", 
                      palette = c("tomato", "steelblue"),
                      add = "mean_se",
                      label = FALSE, position = position_dodge())
    
    x <- which(macros == macro)
    
    plot <- plot +
      labs(x = "Time (months)", y = macro.labs[x],
           tag = tags[index])
    
    plot <- plot +
      stat_pvalue_manual(
        stat.test,
        bracket.nudge.y = 0.5,
        # bracket.shorten = 1,
        bracket.size = 0.5,
        tip.length = 0.01,
        # remove.bracket = TRUE,
        size = 6,
        hide.ns = TRUE,
        # label = "p.adj.signif"
        label = "p.signif"
      )
    plot <- plot +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
      theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1))
    plot
    plotList[[index]] <- plot
    
    index <- index + 1
    
  } # for (macro in macros)
  
} # for (month in months)

plotFileName <- "response_energy_ratio_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in c(1, 4, 7)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]],
               ncol = 3, nrow = 1)
}
dev.off()


