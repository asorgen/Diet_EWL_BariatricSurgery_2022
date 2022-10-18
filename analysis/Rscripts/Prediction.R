##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("Prediction")


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
# args <- "~/git/Diet_EWL_BariatricSurgery_2022"
# args <- c(args, "weight_update.txt")

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
prevModule <- "ExcessWeightLoss"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "updated_weight2.tsv"

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Prep data table for analysis #####
# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Assign SampleID to as dataframe row names
row.names(myTable) = myTable$SampleID

# Convert the timepoint data to numeric values.
myTable$time = as.numeric(myTable$time)


PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
Percent_Loss_kg <- myTable$Percent_Loss_kg
PEWL_kg <- myTable$PEWL_kg
Surgery <- myTable$Surgery

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")

##### Plots for %WL correlations (Kendall) #####
df <- data.frame(PatientID, Timepoint, Percent_Loss_kg)
df <- na.omit(df)

df <- spread(df, Timepoint, Percent_Loss_kg)

MONTHS <- c("ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(1,6,12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)


plotList <- list()
index <- 1
for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "kendall"),
                    add.params = list(fill = "lightgray")
  )
  
  x.lab <- paste0(mo[col1], "-month weight loss (%)")
  y.lab <- paste0(mo[col2], "-month weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot

}


file.path <- paste0(outputDir, "weight_loss_kendall_plots.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

##### Plots for %EWL correlations (Kendall) #####
df <- data.frame(PatientID, Timepoint, PEWL_kg)
df <- na.omit(df)

df <- spread(df, Timepoint, PEWL_kg)

MONTHS <- c("ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(1,6,12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)


plotList <- list()
index <- 1
for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "kendall"),
                    add.params = list(fill = "lightgray")
  )
  
  x.lab <- paste0(mo[col1], "-month excess weight loss (%)")
  y.lab <- paste0(mo[col2], "-month excess weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  
}


file.path <- paste0(outputDir, "excess_weight_loss_kendall_plots.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

##### Plots for %EWL correlations (Kendall) - 12, 18, and 24 months only #####
df <- data.frame(PatientID, Timepoint, PEWL_kg)
df <- na.omit(df)

df <- spread(df, Timepoint, PEWL_kg)

MONTHS <- c("TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)
all_combinations <- all_combinations[,-1]

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "kendall"),
                    add.params = list(fill = "lightgray")
  )
  tag <- tags[tagIndex]
  x.lab <- paste0(mo[col1], "-month excess weight loss (%)")
  y.lab <- paste0(mo[col2], "-month excess weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab, tag = tag)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  tagIndex <- tagIndex + 1
  
}


file.path <- paste0(outputDir, "excess_weight_loss_kendall_plots_12_18_24_only.pdf")
pdf(file.path, width = 7, height = 3.5)
for (i in seq(1, length(plotList), 2)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], 
               ncol = 2, nrow = 1)
}
dev.off()

##### Plots for %EWL correlations (Kendall) - 12, 18, and 24 months only (color) #####
df <- data.frame(PatientID, Timepoint, PEWL_kg, Surgery)
df <- na.omit(df)

df <- spread(df, Timepoint, PEWL_kg)

df$Outcome_12M <- ifelse( df$TWELVE > 50, "Responder",
                          ifelse( df$TWELVE <= 50, "Non-responder", 
                                  NA))

df$Outcome_18M <- ifelse( df$EIGHTEEN > 50, "Responder",
                          ifelse( df$EIGHTEEN <= 50, "Non-responder", 
                                  NA))

df$Outcome_24M <- ifelse( df$TWENTY_FOUR > 50, "Responder",
                          ifelse( df$TWENTY_FOUR <= 50, "Non-responder", 
                                  NA))

df$Outcome_12M_18M <- ifelse((df$Outcome_12M == "Responder" & df$Outcome_18M == "Responder"), "Consistent Responder",
                             ifelse((df$Outcome_12M == "Responder" & df$Outcome_18M == "Non-responder"), "12-month responder only",
                                    ifelse((df$Outcome_12M == "Non-responder" & df$Outcome_18M == "Responder"), "18-month responder only",
                                           ifelse((df$Outcome_12M == "Non-responder" & df$Outcome_18M == "Non-responder"), "Consistent Non-responder",
                                                  NA))))

df$Outcome_12M_24M <- ifelse((df$Outcome_12M == "Responder" & df$Outcome_24M == "Responder"), "Consistent Responder",
                             ifelse((df$Outcome_12M == "Responder" & df$Outcome_24M == "Non-responder"), "12-month responder only",
                                    ifelse((df$Outcome_12M == "Non-responder" & df$Outcome_24M == "Responder"), "24-month responder only",
                                           ifelse((df$Outcome_12M == "Non-responder" & df$Outcome_24M == "Non-responder"), "Consistent Non-responder",
                                                  NA))))

df$Outcome_18M_24M <- ifelse((df$Outcome_18M == "Responder" & df$Outcome_24M == "Responder"), "Consistent Responder",
                             ifelse((df$Outcome_18M == "Responder" & df$Outcome_24M == "Non-responder"), "18-month responder only",
                                    ifelse((df$Outcome_18M == "Non-responder" & df$Outcome_24M == "Responder"), "24-month responder only",
                                           ifelse((df$Outcome_18M == "Non-responder" & df$Outcome_24M == "Non-responder"), "Consistent Non-responder",
                                                  NA))))

df$Consistency <- ifelse((df$Outcome_12M == "Responder" & df$Outcome_18M == "Responder" & df$Outcome_24M == "Responder"), "Consistent Responder",
                             ifelse((df$Outcome_12M == "Non-responder" & df$Outcome_18M == "Non-responder" & df$Outcome_24M == "Non-responder"), "Consistent Non-responder",
                                    "Inconsistent"))

MONTHS <- c("TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)
# all_combinations <- all_combinations[,-1]

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  outcome <- paste0("Outcome_", mo[col1], "M_", mo[col2], "M")
  df2 <- df[!is.na(df[,which(colnames(df) == outcome)]),]
  
  colors <- c("royalblue", "forestgreen", "darkgoldenrod2", "firebrick3")
  plot <- ggscatter(df2, x = m1, y = m2, color = outcome,
                    shape = "Surgery", size = 3.5, # Points shape and size
                    palette = colors
  ); plot
  
  tag <- tags[tagIndex]
  x.lab <- paste0(mo[col1], "-month excess weight loss (%)")
  y.lab <- paste0(mo[col2], "-month excess weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab, tag = tag); plot
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "black", linetype="solid"); plot
  
  plot <- plot + geom_hline(yintercept=50, linetype="dashed", 
                            color = "black", size=0.25); plot
  
  plot <- plot + geom_vline(xintercept=50, linetype="dashed", 
                            color = "black", size=0.25); plot
  
  plot <- plot + theme(legend.title = element_blank(),
                       legend.position = "right"); plot
  
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  tagIndex <- tagIndex + 1
  
}

file.path <- paste0(outputDir, "excess_weight_loss_plots_12_18_24_only_colored.pdf")
pdf(file.path, width = 7, height = 5)
i=1
for (i in 1:length(plotList)) {
  # message(i)
  grid.arrange(plotList[[i]], 
               ncol = 1, nrow = 1)
}
dev.off()



file.path <- paste0(outputDir, "Outcome_changes_12M_18M_24M.tsv")
write.table(df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

summary <- data.frame(table(df$Outcome_12M_18M))
file.path <- paste0(outputDir, "Outcome_summary_12M_18M.tsv")
write.table(summary, file.path,sep="\t",quote = FALSE, row.names = FALSE)

summary <- data.frame(table(df$Outcome_12M_24M))
file.path <- paste0(outputDir, "Outcome_summary_12M_24M.tsv")
write.table(summary, file.path,sep="\t",quote = FALSE, row.names = FALSE)

summary <- data.frame(table(df$Outcome_18M_24M))
file.path <- paste0(outputDir, "Outcome_summary_18M_24M.tsv")
write.table(summary, file.path,sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "outcome_consistency.tsv")
write.table(df, file.path,sep="\t",quote = FALSE, row.names = FALSE)


##### Plots for %WL correlations (Spearman) #####
df <- data.frame(PatientID, Timepoint, Percent_Loss_kg)
df <- na.omit(df)

df <- spread(df, Timepoint, Percent_Loss_kg)

MONTHS <- c("ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(1,6,12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)


plotList <- list()
index <- 1
for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman"),
                    add.params = list(fill = "lightgray")
  )
  
  x.lab <- paste0(mo[col1], "-month weight loss (%)")
  y.lab <- paste0(mo[col2], "-month weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  
}


file.path <- paste0(outputDir, "weight_loss_spearman_plots.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

##### Plots for %EWL correlations (Spearman) #####
df <- data.frame(PatientID, Timepoint, PEWL_kg)
df <- na.omit(df)

df <- spread(df, Timepoint, PEWL_kg)

MONTHS <- c("ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(1,6,12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)


plotList <- list()
index <- 1
for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman"),
                    add.params = list(fill = "lightgray")
  )
  
  x.lab <- paste0(mo[col1], "-month excess weight loss (%)")
  y.lab <- paste0(mo[col2], "-month excess weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  
}


file.path <- paste0(outputDir, "excess_weight_loss_spearman_plots.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

