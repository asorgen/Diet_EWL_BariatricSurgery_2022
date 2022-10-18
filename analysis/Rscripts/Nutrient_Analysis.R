#Author: Alicia Sorgen
#Date: 2022 Mar 01
#Description: Model changes in weight by macronutrients

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

rm(list=ls())





##### To edit #####
# ANALYSIS <- "actigraph"
ANALYSIS <- "ASA24"

module <- paste0("Nutrient_Analysis")



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


##### Compare and calculate average (sd) nutrient intake between responders and non-responders #####
# 4.184 kJ in 1 kcal
myTable$CARB_in_kcal <- (myTable$CARB * 16.7) / 4.184 # 16.7 kJ in 1 gram of carbs
myTable$PROT_in_kcal <- (myTable$PROT * 16.7) / 4.184 # 16.7 kJ in 1 gram of protein
myTable$TFAT_in_kcal <- (myTable$TFAT * 37.7) / 4.184 # 37.7 kJ in 1 gram of fat

myTable$CARB_energy_ratio <- (myTable$CARB_in_kcal / myTable$KCAL) * 100
myTable$PROT_energy_ratio <- (myTable$PROT_in_kcal / myTable$KCAL) * 100
myTable$TFAT_energy_ratio <- (myTable$TFAT_in_kcal / myTable$KCAL) * 100


results <- data.frame()
ratio_results <- data.frame()
summary <- data.frame()

macros <- c("KCAL", "CARB", "PROT", "TFAT")
months <- c("6", "12", "18", "24")
included <- c("0", "1")

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
    
    
    # t-test
    stats.summ <- df %>%
      group_by(Timepoint) %>%
      t_test(MACRO ~ ResponderStatus) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    stats.summ$p_value <- roundP(stats.summ$p.adj)
    stats.summ$ResponderMonth <- month
    stats.summ$Model <- "t-test"
    stats.summ <- stats.summ[,-(which(colnames(stats.summ) == "df"))]
    stats.summ$.y. <- macro
    results <- rbind(results, stats.summ)
    
    
    
    # Wilcoxon
    stats.summ <- df %>%
      group_by(Timepoint) %>%
      wilcox_test(MACRO ~ ResponderStatus) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    stats.summ$p_value <- roundP(stats.summ$p.adj)
    stats.summ$ResponderMonth <- month
    stats.summ$Model <- "Wilcoxon"
    stats.summ$.y. <- macro
    results <- rbind(results, stats.summ)
    
    
    # Averages
    averages <- df %>%
      group_by(Timepoint, ResponderStatus) %>%
      get_summary_stats(MACRO, type = "mean_sd")
    averages$variable <- macro
    averages$ResponderMonth <- month
    
    if (macro == "KCAL") {
      averages$Ratio <- NA
    } else {
      Ratio <- myTable[,which(colnames(myTable) == paste0(macro, "_energy_ratio"))]
      names(Ratio) = myTable[, "SampleID"]
      df$Ratio = Ratio[df$SampleID]
      
      ratios <- df %>%
        group_by(Timepoint, ResponderStatus) %>%
        get_summary_stats(Ratio, type = "mean")
      averages$Ratio <- ratios$mean
      
      # t-tests
      stats.summ <- df %>%
        group_by(Timepoint) %>%
        t_test(Ratio ~ ResponderStatus) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stats.summ$p_value <- roundP(stats.summ$p.adj)
      stats.summ$ResponderMonth <- month
      stats.summ$Model <- "t-test"
      stats.summ <- stats.summ[,-(which(colnames(stats.summ) == "df"))]
      stats.summ$.y. <- paste0(macro, "_energy_ratio")
      ratio_results <- rbind(ratio_results, stats.summ)
      
      # Wilcoxon
      stats.summ <- df %>%
        group_by(Timepoint) %>%
        wilcox_test(Ratio ~ ResponderStatus) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stats.summ$p_value <- roundP(stats.summ$p.adj)
      stats.summ$ResponderMonth <- month
      stats.summ$Model <- "Wilcoxon"
      stats.summ$.y. <- paste0(macro, "_energy_ratio")
      ratio_results <- rbind(ratio_results, stats.summ)
      
    }
    
    summary <- rbind(summary, averages)
    
  }
}

summary$final <- paste0(round(summary$mean, 2), " ± ", round(summary$sd, 2), " (", round(summary$Ratio, 1), "%)")

resultsFileName <- "Macros_by_6M_12M_18M_24M_Outcome__wilcox_ttest.tsv"
file.path <- paste0(outputDir, resultsFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

summaryFileName <- "Macros_by_6M_12M_18M_24M_Outcome__Mean_SD_EnergyRatio.tsv"
file.path <- paste0(outputDir, summaryFileName)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)

energyRatioFileName <- "EnergyRatio_by_6M_12M_18M_24M_Outcome__wilcox_ttest.tsv"
file.path <- paste0(outputDir, energyRatioFileName)
write.table(ratio_results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Compare responder/non-responder nutrient intake over time #####
results <- data.frame()
summary <- data.frame()

months <- c("6", "12", "18", "24")
included <- c("0", "1")
macros <- c("KCAL", "CARB", "PROT", "TFAT")
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
  df$TFAT = TFAT[df$SampleID]
  
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
      tTable$Macronutrient <- macro
      tTable$Group <- aStatus
      rownames(tTable) <- NULL
      
      results <- rbind(results, tTable)
      
    }
    
    
  }
  
  
}


macronutrientLMFileName <- "Macros_by_Timepoint_grouped_by_6M_12M_18M_24M_Outcome__MixedLinearModels.tsv"
file.path <- paste0(outputDir, macronutrientLMFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Correlation between %EWL and nutrient intake #####
results <- data.frame()
summary <- data.frame()
plotList <- list()

palette <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
months <- c("6", "12", "18", "24")
included <- c("0", "1")
macros <- c("KCAL", "CARB", "PROT", "TFAT")
index <- 1
monthindex <- 1

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
  df$TFAT = TFAT[df$SampleID]
  
  for (m in 2:length(included)) {
    for (macro in macros) {
      
      macroCol <- df[,which(colnames(df) == macro)]
      pID <- df$PatientID
      time <- df$Timepoint
      status <- df$ResponderStatus
      pewl <- df$PEWL_kg
      
      df2 <- data.frame(pID, time, status, macroCol, pewl)
      df2 <- na.omit(df2)
      sampleSize <- length(unique(df2$pID))
      label <- paste0("Patients classified at ", month, " months (n = ", sampleSize, ")")
      # print(label)
      
      df3 <- df2[df2$time == included[m],]
      
      # Kendall
      stat.test <- df3 %>%
        # group_by(status) %>%
        cor_test(pewl, macroCol, method = "kendall") %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      
      stat.test$var2 <- macro
      stat.test$p_value <- roundP(stat.test$p.adj)
      stat.test$Timepoint <- included[m]
      stat.test$ResponderMonth <- month
      results <- rbind(results, stat.test)
      
      # Make correlation plots
      title.lab <- label
      caption.lab <- "* Unadjusted Kendall Correlation p values"
      # caption.lab <- "* Unadjusted Pearson Correlation p values"
      xLabels <- c("Energy (kcal)", "Carbohydrates (g)", "Protein (g)", "Total fat (g)")
      x.lab <- paste0(xLabels[which(macros == macro)], " at ", included[m], " month(s)")
      y.lab = paste0("%EWL (BL - ", included[m], "M)")
      
      i <- which(macros == macro)
      plot <- ggscatter(df3, x = "macroCol", y = "pewl",
                        color = palette[i],
                        shape = 21, size = 2.5, # Points shape and size
                        add = "reg.line",  # Add regression line
                        # add.params = list(color = "black"), # Customize reg. line
                        conf.int = TRUE # Add confidence interval
                        # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                        # cor.coeff.args = list(method = "kendall"),
                        # cor.coef.coord = c(-1, 400),
                        # facet.by = "Timepoint" #, label = "PatientID", repel = TRUE
      )
      
      
      plot <- plot + scale_x_log10()
      plot <- plot + stat_cor(method = "kendall")
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      
      plot <- plot + labs(caption = caption.lab)
      plot <- plot + labs(title = title.lab); plot
      plotList[[index]] <- plot
      
      
      # print(plot)
      
      
      index <- index + 1
      
    }
  }
  
  monthindex <- monthindex + 1
  
}

plotFileName <- "Macro_by_PEWL_6M_12M_18M_24M_Outcome_kendall_plots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 4)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()

kendallFileName <- "Macro_by_PEWL_6M_12M_18M_24M_Outcome_kendall_results.tsv"
file.path <- paste0(outputDir, kendallFileName)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Correlation between %EWL and nutrient intake (separated by responder group) #####
results <- data.frame()
summary <- data.frame()
plotList <- list()

# palette <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
months <- c("6", "12", "18", "24")
included <- c("0", "1")
macros <- c("KCAL", "CARB", "PROT", "TFAT")
index <- 1
monthindex <- 1

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
  df$TFAT = TFAT[df$SampleID]
  
  for (m in 2:length(included)) {
    for (macro in macros) {
      
      macroCol <- df[,which(colnames(df) == macro)]
      pID <- df$PatientID
      time <- df$Timepoint
      status <- df$ResponderStatus
      pewl <- df$PEWL_kg
      
      df2 <- data.frame(pID, time, status, macroCol, pewl)
      df2 <- na.omit(df2)
      sampleSize <- length(unique(df2$pID))
      label <- paste0("Patients classified at ", month, " months (n = ", sampleSize, ")")
      # print(label)
      
      df3 <- df2[df2$time == included[m],]
      
      # Kendall
      stat.test <- df3 %>%
        group_by(status) %>%
        cor_test(pewl, macroCol, method = "kendall") %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      
      stat.test$var2 <- macro
      stat.test$p_value <- roundP(stat.test$p.adj)
      stat.test$Timepoint <- months[m]
      stat.test$ResponderMonth <- month
      results <- rbind(results, stat.test)
      
      # Make correlation plots
      title.lab <- label
      caption.lab <- "* Unadjusted Kendall Correlation p values"
      # caption.lab <- "* Unadjusted Pearson Correlation p values"
      
      xLabels <- c("Energy (kcal)", "Carbohydrates (g)", "Protein (g)", "Total fat (g)")
      x.lab <- paste0(xLabels[which(macros == macro)], " at ", included[m], " month(s)")
      y.lab = paste0("%EWL (BL - ", included[m], "M)")
      
      plot <- ggscatter(df3, x = "macroCol", y = "pewl", color = "status",
                        palette = c("#C41522", "#2BCE2A"),
                        shape = 21,
                        size = 2.5, # Points shape and size
                        add = "reg.line",  # Add regression line
                        conf.int = TRUE # Add confidence interval
      )
      
      
      plot <- plot + scale_x_log10()
      plot <- plot + stat_cor(aes(color = status), method = "kendall")
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      plot <- plot + theme(legend.title = element_blank())
      plot <- plot + labs(title = title.lab)
      plot <- plot + labs(caption = caption.lab)
      
      plotList[[index]] <- plot
      
      
      # print(plot)
      
      
      index <- index + 1
      
    }
  }
  
  monthindex <- monthindex + 1
  
}

plotFileName2 <- "Macro_by_PEWL_6M_12M_18M_24M_Outcome_grouped_by_outcome_kendall_plots.pdf"
file.path <- paste0(outputDir, plotFileName2)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 4)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()

kendallFileName2 <- "Macro_by_PEWL_6M_12M_18M_24M_Outcome_grouped_by_outcome_kendall_results.tsv"
file.path <- paste0(outputDir, kendallFileName2)
write.table(results, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Make and generate data tables for nutrients, PEWL, and responder status at 6, 12, 18, 24 months #####
months <- c("6", "12", "18", "24")
included <- c("0", "1")

for (month in months) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  KCAL <- myTable$KCAL
  CARB <- myTable$CARB
  PROT <- myTable$PROT
  TFAT <- myTable$TFAT
  PEWL_kg <- myTable$PEWL_kg
  Baseline_kg <- myTable$Baseline_kg
  Percent_Loss_kg <- myTable$Percent_Loss_kg
  
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
  names(Baseline_kg) = myTable[, "SampleID"]
  names(Percent_Loss_kg) = myTable[, "SampleID"]
  
  # Add info to each row
  df$KCAL = KCAL[df$SampleID]
  df$CARB = CARB[df$SampleID]
  df$PROT = PROT[df$SampleID]
  df$TFAT = TFAT[df$SampleID]
  df$Baseline_kg = Baseline_kg[df$SampleID]
  df$Percent_Loss_kg = Percent_Loss_kg[df$SampleID]
  
  # Output resulting dataframe for each responder month
  file.path <- paste0(outputDir, paste0("Nutrient_PEWL_Outcomes_", month, "M.tsv"))
  write.table(df, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (month in months)

##### Plots for Kendall correlations between macronutrient & % weight loss and absolute weight loss (kg) classifying patients by heaviest, middle, and lightest baseline weight #####
plotList2 <- list()
plotList3 <- list()

months <- c("6", "12", "18", "24")
included <- c("0", "1")
macros <- c("KCAL", "CARB", "PROT", "TFAT")
index <- 1
monthindex <- 1

for (month in months) {
  
  # month <- 12
  
  IDS <- myTable$PatientID[which(myTable$ResponderStatus %in% c(paste0(month, "-month responder"), paste0(month, "-month non-responder")))]
  df <- myTable[myTable$PatientID %in% IDS,]
  
  # included <- c(0, 1, 6, 12)
  included <- c(included, month)
  df <- df[df$time %in% included,]
  
  df.BL <- df[df$time == 0,]
  
  kg <- df.BL$Baseline_kg
  n <- length(kg)
  ordered.kg <- sort(kg)
  
  q3 <- (2/3) * (n + 1)
  up1 <- ceiling(q3)
  down1 <- floor(q3)
  q3 <- (ordered.kg[up1] + ordered.kg[down1]) / 2
  
  q2 <- (1/3) * (n + 1)
  up1 <- ceiling(q2)
  down1 <- floor(q2)
  q2 <- (ordered.kg[up1] + ordered.kg[down1]) / 2
  
  df$Baseline.3rd <- ifelse(df$Baseline_kg >= q3, "Heaviest 3rd",
                            ifelse(df$Baseline_kg < q3 & df$Baseline_kg >= q2 , "Middle 3rd",
                                   "Lightest 3rd"))
  
  df$Baseline.3rd <- as.factor(df$Baseline.3rd)
  df$Baseline.3rd <- factor(df$Baseline.3rd, levels = c("Lightest 3rd", "Middle 3rd", "Heaviest 3rd"))
  
  for (m in 2:length(included)) {
    
    for (macro in macros) {
      
      macroCol <- df[,which(colnames(df) == macro)]
      pID <- df$PatientID
      time <- df$time
      status <- df$ResponderStatus
      pewl <- df$PEWL_kg
      thirds <- df$Baseline.3rd
      pwl <- df$Percent_Loss_kg
      awl <- df$Loss_from_BL_kg
      
      df2 <- data.frame(pID, time, status, macroCol, pewl, thirds, pwl, awl)
      df2 <- na.omit(df2)
      sampleSize <- length(unique(df2$pID))
      label <- paste0("Patients classified at ", month, " months (n = ", sampleSize, ")")
      # print(label)
      
      df3 <- df2[df2$time == included[m],]
      
      # Make correlation plots for PWL
      title.lab <- label
      x.lab <- paste0(xLabels[which(macros == macro)], " at ", included[m], " month(s)")
      y.lab = paste0("Weight loss (%): BL - ", included[m], "M")
      
      i <- which(macros == macro)
      plot <- ggscatter(df3, x = "macroCol", y = "pwl",
                        color = palette[i],
                        shape = 1, size = 3 # Points shape and size
      )
      
      # plot <- plot + scale_x_log10()
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      
      plot <- plot + stat_cor(method = "kendall")
      plot <- plot +geom_smooth(color="black", method="lm", se=TRUE, linetype="solid", fullrange=TRUE, fill = "lightgray")
      plot <- plot + theme(legend.title = element_blank())
      # plot <- plot + labs(title = title.lab)
      
      plotList2[[index]] <- plot
      
      
      
      
      # Make correlation plots for WL
      y.lab = paste0("Weight loss (kg): BL - ", included[m], "M")
      
      i <- which(macros == macro)
      plot <- ggscatter(df3, x = "macroCol", y = "awl",
                        color = palette[i],
                        shape = 20, size = 3 # Points shape and size
      )
      
      # plot <- plot + scale_x_log10()
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      
      plot <- plot + stat_cor(method = "kendall")
      plot <- plot +geom_smooth(color="black", method="lm", se=FALSE, linetype="solid", fullrange=TRUE)
      plot <- plot + theme(legend.title = element_blank())
      plot <- plot + labs(title = title.lab)
      
      plotList3[[index]] <- plot
      
      index <- index + 1
      
    } # for (macro in macros)
  } # for (m in 1:monthindex)
  
  monthindex <- monthindex + 1
  
} # for (month in months)


plotFileName3 <- "Macro_by_PWL_6M_12M_18M_24M_Outcome_kendall_plots.pdf"
file.path <- paste0(outputDir, plotFileName3)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 4)) {
  grid.arrange(plotList2[[i]], plotList2[[i+1]],
               plotList2[[i+2]], plotList2[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()


plotFileName4 <- "Macro_by_WL_6M_12M_18M_24M_Outcome_kendall_plots.pdf"
file.path <- paste0(outputDir, plotFileName4)
pdf(file.path, width = 10, height = 10)
for (i in seq(1, length(plotList), 4)) {
  grid.arrange(plotList3[[i]], plotList3[[i+1]],
               plotList3[[i+2]], plotList3[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()

##### Macronutrients over time (linear model) #####
results <- data.frame()

macros <- c("KCAL", "CARB", "PROT", "TFAT")
months <- c("6", "12", "18", "24")
included <- c("0", "1")

dFrame <- data.frame()

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
    
    for (status in c("Responder", "Non-responder")) {
      
      R <- df[df$ResponderStatus == status,]
      uni.lm <- lm(R$MACRO ~ R$Timepoint)
      uni.lm.sum <- summary(uni.lm)
      
      Status <- vector()
      Comparison <- vector()
      Macronutrient <- vector()
      pValue <- vector()
      StatusMonth <- vector()
      index <- 1
      
      for (j in 2:length(included)) {
        
        Status[index] <- status
        StatusMonth[index] <- month
        Comparison[index] <- paste0("BL v ", included[j], "M")
        Macronutrient[index] <- macro
        pValue[index] <- uni.lm.sum$coefficients[j,4]
        index <- index + 1
        
      } # for (j in 2:length(included))
      
      result <- data.frame(Status, StatusMonth, Comparison, Macronutrient, pValue)
      result$Adj_pValue <- p.adjust(result$pValue, method = "BH")
      result$Significance <- sigStars(result$Adj_pValue)
      
      dFrame <- rbind(dFrame, result)
      
    } # for (status in c("Responder", "Non-responder"))
    
  } # for (macro in macros)
  
} # for (month in months)

macroLMFileName <- "Macros_by_Timepoint_grouped_by_6M_12M_18M_24M_Outcome__PairwiseLinearModels.tsv"
file.path <- paste0(outputDir, macroLMFileName)
write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)

