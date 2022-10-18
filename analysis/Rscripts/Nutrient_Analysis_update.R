#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Analysis of nutreint intake between responders for 24-month patients at 12 and 18 months.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("Nutrient_Analysis_update")
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

##### Set up for correlation plots #####
MONTHS <- c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(0,1,6,12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg


tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")

##### Kendall plots (without outliers) #####
message("\n##### Kendall plots (without outliers) #####\n")
plotList <- list()
index <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'kendall')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_kendall_by_timepoint_no_outliers.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()





##### Kendall plots (with outliers) #####
message("\n##### Kendall plots (with outliers) #####\n")
plotList <- list()
index <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    # df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'kendall')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_kendall_by_timepoint_with_outliers.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

##### Spearman plots (without outliers) #####
message("\n##### Spearman plots (without outliers) #####\n")
plotList <- list()
index <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'spearman')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
                        )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_spearman_by_timepoint_no_outliers.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()





##### Spearman plots (with outliers) #####
message("\n##### Spearman plots (with outliers) #####\n")
plotList <- list()
index <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    # df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'spearman')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_spearman_by_timepoint_with_outliers.pdf")
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()

##### Set up for final tables #####
message("\n##### Set up for final tables #####\n")

macros <- c("KCAL", "CARB", "PROT", "TFAT")
months <- c("6", "12", "18", "24")
included <- c("0", "1")

dFrame <- data.frame()
t.pValue <- vector()
lm.pValue <- vector()
outcomeMonth <- vector()
monthTested <- vector()
r.squared <- vector()
nNon <- vector()
nRes <- vector()
macronutrient <- vector()
non_avg <- vector()
non_sd <- vector()
res_avg <- vector()
res_sd <- vector()
wilcox.pValue <- vector()

for (month in months) {
  
  included <- c(included, month)
  
  for (macro in macros) {
    
    PatientID <- myTable$PatientID
    Timepoint <- myTable$time
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    PEWL_kg <- myTable$PEWL_kg
    
    df <- data.frame(PatientID, Timepoint, PEWL_kg)
    df <- na.omit(df)
    df$Timepoint <- paste0("M", df$Timepoint)
    
    df <- spread(df, Timepoint, PEWL_kg)
    m <- which(colnames(df) == paste0("M", month))
    
    df$ResponderStatus <- ifelse(df[,m] >= 50, "Responder",
                                 ifelse(df[,m] < 50, "Non-responder", 
                                        NA))
    n_Responders <- length(which(df$ResponderStatus == "Responder"))
    n_Nonresponders <- length(which(df$ResponderStatus == "Non-responder"))
    
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
    
    names(MACRO) = myTable[, "SampleID"]
    
    # Add info to each row
    df$MACRO = MACRO[df$SampleID]
    df <- na.omit(df)
    
    n <- length(unique(df$PatientID))
    index <- 1
    
    for (mo in included) {
      
      label <-  paste0(macro, " at month ",mo, "\nPatient outcome status when classified ", month.labs[which(included == month)], " (n = ", n, ")")
      monthTested[index] <- mo
      outcomeMonth[index] <- month
      df2 <- df[df$Timepoint == mo,]
      
      lm.sum <- summary(lm(MACRO ~ ResponderStatus, data = df2))
      lm.pValue[index] <- lm.sum$coefficients[2,4]
      r.squared[index] <- lm.sum$r.squared
      # message(paste0("Linear model: R^2 = ", r2, "; p = ", p))
      
      t <- t.test(MACRO ~ ResponderStatus, data = df2)
      t.pValue[index] <- t$p.value
      # message(paste0("t-test: p = ", p))
      
      wilcox <- wilcox.test(MACRO ~ ResponderStatus, data = df2)
      wilcox.pValue[index] <- wilcox$p.value
      
      table <- table(df2$ResponderStatus)
      nNon[index] <- table[1]
      nRes[index] <- table[2]
      macronutrient[index] <- macro
      
      averages <- df2 %>%
        group_by(ResponderStatus) %>%
        get_summary_stats(MACRO, type = "mean_sd")
      
      non_avg[index] <- averages$mean[1]
      non_sd[index] <- averages$sd[1]
      res_avg[index] <- averages$mean[2]
      res_sd[index] <- averages$sd[2]
      
      index <- index + 1
      
    } # for (mo in included)
    
    lm.adj <- p.adjust(lm.pValue, method = "BH")
    t.adj <- p.adjust(t.pValue, method = "BH")
    wilcox.adj <- p.adjust(wilcox.pValue, method = "BH")
    
    dFrame.temp <- data.frame(macronutrient, outcomeMonth, monthTested, lm.pValue, lm.adj, r.squared, 
                              t.pValue, t.adj,
                              wilcox.pValue, wilcox.adj,
                              nNon, non_avg, non_sd,
                              nRes, res_avg, res_sd)
    dFrame.temp$n_Responders <- n_Responders
    dFrame.temp$n_Nonresponders <- n_Nonresponders
    
    dFrame <- rbind(dFrame, dFrame.temp)
  } # for (macro in macros)
  
} # for (month in months)


##### Make final tables #####
message("\n##### Make final tables #####\n")
for (oMonth in unique(dFrame$outcomeMonth)) {
  
  df <- dFrame[dFrame$outcomeMonth == oMonth,]
  
  for (macro in macros) {
    
    df2 <- df[df$macronutrient == macro,]
    
    TP <- vector()
    Responder <- vector()
    Nonresponder <- vector()
    pValue <- vector()
    pAdj <- vector()
    r2 <- vector()
    index <- 1
    
    TP[index] <- NA
    Responder[index] <- paste0("(n = ", df2$n_Responders[1], ")")
    Nonresponder[index] <- paste0("(n = ", df2$n_Nonresponders[1], ")")
    pValue[index] <- NA
    pAdj[index] <- NA
    r2[index] <- NA
    
    for (i in 1:nrow(df2)) {
      
      index <- index + 1
      
      TP[index] <- df2$monthTested[i]
      Responder[index] <- paste0(round(df2$res_avg[i], 2), " ± ", round(df2$res_sd[i], 2))
      Nonresponder[index] <- paste0(round(df2$non_avg[i], 2), " ± ", round(df2$non_sd[i], 2))
      pValue[index] <- df2$wilcox.pValue[i]
      pAdj[index] <- df2$wilcox.adj[i]
      r2[index] <- df2$r.squared[i]
      
    } # for (i in 1:nrow(df2))
    
    temp <- data.frame(TP, Responder, Nonresponder, pValue, pAdj, r2)
    names(temp)[names(temp) == "Responder"] <- paste0("Responder_", macro)
    names(temp)[names(temp) == "Nonresponder"] <- paste0("Nonresponder_", macro)
    names(temp)[names(temp) == "pValue"] <- paste0("pValue_", macro)
    names(temp)[names(temp) == "pAdj"] <- paste0("pAdj_", macro)
    names(temp)[names(temp) == "r2"] <- paste0("r2_", macro)
    
    if (macro == macros[1]) {
      resultsTable <- temp
    } else {
      resultsTable <- cbind(resultsTable, temp)
    } # if (macro == macros[1])
    
  } # for (macro in macros)
  
  file.path <- paste0(outputDir, "Wilcoxon_Nutrient_Table_", oMonth, "month_status.tsv")
  write.table(resultsTable, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (oMonth in unique(dFrame$outcomeMonth))

finalTableNames <- c("S2", "2a", "S2", "2b")
tableIndex <- 0
for (oMonth in unique(dFrame$outcomeMonth)) {
  
  tableIndex <- tableIndex + 1
  df <- dFrame[dFrame$outcomeMonth == oMonth,]
  
  for (macro in macros) {
    
    df2 <- df[df$macronutrient == macro,]
    
    TP <- vector()
    Responder <- vector()
    Nonresponder <- vector()
    pValue <- vector()
    pAdj <- vector()
    r2 <- vector()
    index <- 1
    
    TP[index] <- NA
    Responder[index] <- paste0("(n = ", df2$n_Responders[1], ")")
    Nonresponder[index] <- paste0("(n = ", df2$n_Nonresponders[1], ")")
    pValue[index] <- NA
    pAdj[index] <- NA
    r2[index] <- NA
    
    for (i in 1:nrow(df2)) {
      
      index <- index + 1
      
      TP[index] <- df2$monthTested[i]
      Responder[index] <- paste0(round(df2$res_avg[i], 2), " ± ", round(df2$res_sd[i], 2))
      Nonresponder[index] <- paste0(round(df2$non_avg[i], 2), " ± ", round(df2$non_sd[i], 2))
      pValue[index] <- df2$wilcox.pValue[i]
      pAdj[index] <- df2$wilcox.adj[i]
      r2[index] <- df2$r.squared[i]
      
    } # for (i in 1:nrow(df2))
    
    temp <- data.frame(TP, Responder, Nonresponder, pValue, pAdj, r2)
    temp$pValue <- ifelse(temp$pValue < 0.001, "<0.001",
                          round(temp$pValue, 3))
    temp$r2 <- ifelse(temp$r2 < 0.001, "<0.001",
                          round(temp$r2, 3))
    names(temp)[names(temp) == "Responder"] <- paste0("Responder_", macro)
    names(temp)[names(temp) == "Nonresponder"] <- paste0("Nonresponder_", macro)
    names(temp)[names(temp) == "pValue"] <- paste0("pValue_", macro)
    names(temp)[names(temp) == "pAdj"] <- paste0("pAdj_", macro)
    names(temp)[names(temp) == "r2"] <- paste0("r2_", macro)
    
    if (macro == macros[1]) {
      resultsTable <- temp
    } else {
      resultsTable <- cbind(resultsTable, temp)
    } # if (macro == macros[1])
    
  } # for (macro in macros)
  
  
  resultsTable2 <- resultsTable[,-c(5,7,11,13,17,19,23)]
  resultsTable2 <- rbind(colnames(resultsTable2), resultsTable2)
  resultsTable2[1,] <- gsub("_KCAL", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("_CARB", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("_PROT", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("_TFAT", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("TP", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("pValue", "p value", resultsTable2[1,])
  resultsTable2[1,] <- gsub("Nonresponder", "Non-responder", resultsTable2[1,])
  
  resultsTable2 <- rbind(colnames(resultsTable2), resultsTable2)
  resultsTable2[1,] <- gsub("Responder_", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("Nonresponder_", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("pValue_", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("r2_", "", resultsTable2[1,])
  resultsTable2[1,] <- gsub("TP", "Timepoint", resultsTable2[1,])
  resultsTable2[1,] <- gsub("KCAL", "Calories (kcal)", resultsTable2[1,])
  resultsTable2[1,] <- gsub("CARB", "Carbohydrates (g)", resultsTable2[1,])
  resultsTable2[1,] <- gsub("PROT", "Protein (g)", resultsTable2[1,])
  resultsTable2[1,] <- gsub("TFAT", "Total fat (g)", resultsTable2[1,])
  
  resultsTable2$TP <- gsub("0", "Baseline", resultsTable2$TP)
  resultsTable2$TP[5] <- paste0("Post-op ", resultsTable2$TP[5], " month")
  resultsTable2$TP[6:nrow(resultsTable2)] <- paste0("Post-op ", resultsTable2$TP[6:nrow(resultsTable2)], " months")
  
  
  file.path <- paste0(outputDir, "Table_", finalTableNames[tableIndex], "_", oMonth, "M.tsv")
  write.table(resultsTable2, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (oMonth in unique(dFrame$outcomeMonth))


##### Kendall & Spearman plots macronutrient v %EWL at each time point #####
message("\n##### Kendall plots macronutrient v %EWL at each time point #####\n")
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg

plotList <- list()
index <- 1

Macro <- vector()
Month <- vector()
Kendall_p <- vector()
Kendall_R <- vector()
Spearman_p <- vector()
Spearman_R <- vector()

tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    df2$Title <- month.labs[i]
    
    # df2$MACRO <- remove_outliers(df2$MACRO)
    # df2 <- na.omit(df2)
    # outliers <- n - nrow(df2)
    # caption.lab <- paste0("*", outliers, " outliers removed")
    
    myKen <- cor.test( df2$MACRO, df2$PEWL, method = "kendall" )
    mySpear <- cor.test( df2$MACRO, df2$PEWL, method = "spearman" )
    
    Macro[index] <- macro
    Month[index] <- month.labs[i]
    Kendall_p[index] <- myKen$p.value
    Kendall_R[index] <- myKen$estimate
    Spearman_p[index] <- mySpear$p.value
    Spearman_R[index] <- mySpear$estimate
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'kendall')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    # plot <- plot + stat_cor(label.x = xMax, label.y = yMax)
    
    plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab, 
                        # title = title.lab, 
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
                        tag = tag
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))

dFrame <- data.frame(Macro, Month, Kendall_p, Kendall_R, Spearman_p, Spearman_R)

file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall_spearman.tsv")
write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()

##### Kendall plots macronutrient v %EWL at each time point (outliers removed) #####
message("\n##### Kendall plots macronutrient v %EWL at each time point (outliers removed) #####\n")
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    df2$Title <- month.labs[i]
    
    df2$MACRO <- remove_outliers(df2$MACRO)
    df2 <- na.omit(df2)
    outliers <- n - nrow(df2)
    caption.lab <- paste0("*", outliers, " outliers removed")
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'kendall')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        # title = title.lab, 
                        tag = tag
                        # ,subtitle = subtitle.lab
                        ,caption = caption.lab
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall_no_outliers.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()


##### Kendall plots macronutrient v %EWL at each time point 12 and 24 months (outliers removed) #####
message("\n##### Kendall plots macronutrient v %EWL at each time point 12 and 24 months (outliers removed) #####\n")
MONTHS <- c("BL", "TWELVE", "TWENTY_FOUR")
mo <- c(0,12,24)
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  # tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    # df2$Title <- month.labs[i]
    
    df2$MACRO <- remove_outliers(df2$MACRO)
    df2 <- na.omit(df2)
    outliers <- n - nrow(df2)
    caption.lab <- paste0("*", outliers, " outliers removed")
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'kendall')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    # plot <- plot + stat_cor(label.x = xMax, label.y = yMax)
    
    # plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        # title = title.lab, 
                        tag = tag
                        # ,subtitle = subtitle.lab
                        ,caption = caption.lab
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall_12_24_no_outliers.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()


##### Kendall plots macronutrient v %EWL at each time point 12 and 24 months #####
message("\n##### Kendall plots macronutrient v %EWL at each time point 12 and 24 months #####\n")
MONTHS <- c("BL", "TWELVE", "TWENTY_FOUR")
mo <- c(0,12,24)
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  # tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    # df2$Title <- month.labs[i]
    
    # df2$MACRO <- remove_outliers(df2$MACRO)
    # df2 <- na.omit(df2)
    # outliers <- n - nrow(df2)
    # caption.lab <- paste0("*", outliers, " outliers removed")
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'kendall')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    # plot <- plot + stat_cor(label.x = xMax, label.y = yMax)
    
    # plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab, 
                        # title = title.lab, 
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
                        tag = tag
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall_12_24_mo_only.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()

##### Kendall plots 12 v. 24 months (with outliers) #####
message("\n##### Kendall plots (with outliers) #####\n")
MONTHS <- c("TWELVE", "TWENTY_FOUR")
mo <- c(12,24)
y <- 2
all_combinations <- combn(MONTHS,y)
month.labs <- c("12 months post-op",
                "24 months post-op")

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    # df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    tag <- tags[tagIndex]
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'kendall')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab, tag = tag,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plot
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_kendall_by_timepoint_with_outliers_12_24_only.pdf")
pdf(file.path, width = 6, height = 6)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], 
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()


##### Kendall plots macronutrient v %EWL at each time point 12, 18, and 24 months (outliers removed) #####
message("\n##### Kendall plots macronutrient v %EWL at each time point 12 and 24 months (outliers removed) #####\n")
MONTHS <- c("BL", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(0,12,18,24)
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  # tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    # df2$Title <- month.labs[i]
    
    df2$MACRO <- remove_outliers(df2$MACRO)
    df2 <- na.omit(df2)
    outliers <- n - nrow(df2)
    caption.lab <- paste0("*", outliers, " outliers removed")
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'kendall')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "kendall", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    # plot <- plot + stat_cor(label.x = xMax, label.y = yMax)
    
    # plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        # title = title.lab, 
                        tag = tag
                        # ,subtitle = subtitle.lab
                        ,caption = caption.lab
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_kendall_12_18_24_no_outliers.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()

##### Spearman plots macronutrient v %EWL at each time point 12, 18, and 24 months (outliers removed) #####
message("\n##### Spearman plots macronutrient v %EWL at each time point 12 and 24 months (outliers removed) #####\n")
MONTHS <- c("BL", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(0,12,18,24)
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  # tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    # df2$Title <- month.labs[i]
    
    df2$MACRO <- remove_outliers(df2$MACRO)
    df2 <- na.omit(df2)
    outliers <- n - nrow(df2)
    caption.lab <- paste0("*", outliers, " outliers removed")
    
    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'spearman')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    # plot <- plot + stat_cor(label.x = xMax, label.y = yMax)
    
    # plot <- plot + facet_grid(. ~ Title)
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        # title = title.lab, 
                        tag = tag
                        # ,subtitle = subtitle.lab
                        ,caption = caption.lab
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_spearman_12_18_24_no_outliers.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()





##### Spearman plots macronutrient v %EWL at each time point 12, 18, and 24 months (outliers included) #####
MONTHS <- c("BL", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(0,12,18,24)
xLabels <- c("Energy (kcal)", "Carbohydrate (g)", "Protein (g)", "Total fat (g)")

PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg



plotList <- list()
index <- 1

tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")
tagIndex <- 1
for (i in 2:length(MONTHS)) {
  
  monthName <- MONTHS[i]
  numericalMonth <- mo[i]
  
  # tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO, PEWL)
    
    df2 <- df[df$Timepoint == monthName,]
    df2 <- na.omit(df2)
    n <- nrow(df2)
    # df2$Title <- month.labs[i]
    
    df2 <- na.omit(df2)

    tag <- tags[tagIndex]
    title.lab <- month.labs[i]
    x.lab <- xLabels[which(macros == macro)]
    y.lab <- "Excess weight loss (%)"
    # caption.lab <- paste0("cor_test(", macro, " ~ %EWL, method = 'spearman')")
    subtitle.lab <- paste0("(n = ", n, ")")
    
    plot <- ggscatter(df2, x = "MACRO", y = "PEWL",
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman", p.accuracy = 0.001, label.x.npc = 0.25, label.y.npc = 1),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab,
                        # title = title.lab, 
                        tag = tag
    )
    plot <- plot + theme(
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text.x = element_text(size = 13)
    )
    plot 
    
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
  
} # for (i in 1:ncol(all_combinations))






file.path <- paste0(outputDir, "macronutrient_v_PEWL_spearman_12_18_24_with_outliers.pdf")
pdf(file.path, width = 12, height = 3)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], plotList[[i+2]], plotList[[i+3]],
               ncol = 4, nrow = 1)
}
dev.off()




##### Spearman plots 12 v. 24 months (with outliers) #####
message("\n##### Spearman plots (with outliers) #####\n")
MONTHS <- c("TWELVE", "TWENTY_FOUR")
mo <- c(12,24)
y <- 2
all_combinations <- combn(MONTHS,y)
month.labs <- c("12 months post-op",
                "24 months post-op")

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    # df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0(diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    tag <- tags[tagIndex]
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    caption.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'spearman')")
    subtitle.lab <- paste0(outlier_message, " (unadj p values)")
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab, tag = tag,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        # ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plot
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_spearman_by_timepoint_with_outliers_12_24_only.pdf")
pdf(file.path, width = 6, height = 6)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], 
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()



##### Spearman plots 12 v. 24 months (with outliers) #####
message("\n##### Spearman plots (with outliers) #####\n")
MONTHS <- c("TWELVE", "TWENTY_FOUR")
mo <- c(12,24)
y <- 2
all_combinations <- combn(MONTHS,y)
month.labs <- c("12 months post-op",
                "24 months post-op")

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  tagIndex <- 1
  
  for (macro in macros) {
    
    MACRO <- myTable[,which(colnames(myTable) == macro)]
    
    df <- data.frame(PatientID, Timepoint, MACRO)
    df <- df[df$Timepoint %in% c(m1, m2),]
    df <- na.omit(df)
    data_w_outliers <- nrow(df)
    
    df$MACRO <- remove_outliers(df$MACRO)
    df <- na.omit(df)
    data_wout_outliers <- nrow(df)
    diff <- data_w_outliers - data_wout_outliers
    
    outlier_message <- paste0("*", diff, " outlier(s) removed")
    
    df <- spread(df, Timepoint, MACRO)
    
    tag <- tags[tagIndex]
    title.lab <- xLabels[which(macros == macro)]
    x.lab <- month.labs[which(MONTHS == m1)]
    y.lab <- month.labs[which(MONTHS == m2)]
    subtitle.lab <- paste0("cor_test(", m1, " ~ ", m2, " method = 'spearman')")
    caption.lab <- paste0(outlier_message)
    
    plot <- ggscatter(df, x = m1, y = m2,
                      shape = 21, size = 2.5, # Points shape and size
                      add = "reg.line",  # Add regression line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman"),
                      add.params = list(fill = "lightgray")
    )
    
    plot <- plot + labs(x=x.lab, y = y.lab, tag = tag,
                        title = title.lab
                        # ,subtitle = subtitle.lab
                        ,caption = caption.lab
    )
    
    plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", lty = 2)
    plot
    plotList[[index]] <- plot
    index <- index + 1
    tagIndex <- tagIndex + 1
    
  } # for (macro in macros)
  
} # for (i in 1:ncol(all_combinations))


file.path <- paste0(outputDir, "macronutrient_spearman_by_timepoint_12_24_only_no_outliers.pdf")
pdf(file.path, width = 6, height = 6)
for (i in seq(1, length(plotList), 4)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], 
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()


