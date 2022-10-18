#Author: Alicia Sorgen
#Date: 2021 Dec 01
#Description: Mixed linear models on metadata variables.

R <- sessionInfo()
message(R$R.version$version.string)

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "ASA24"
module <- paste0("MetaLinearModeling")

included <- c(  0
                , 1
                , 6
                , 12
                , 18
                , 24
)

##### Libraries #####
library(nlme); message("nlme: Version", packageVersion("nlme"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))


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
prevModule <- paste0("WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
metaFile <- "metadata.tsv"
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")
statsOutput <- "MLM_metadata_statistics.tsv"

##### Set up data input #####
metaVariables <- c("KCAL", "CARB", "PROT", "TFAT")
xLabels <- c("Energy (kcal)", "Carbohydrates (g)", "Protein (g)", "Total fat (g)")

tags <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")

file.path <- paste0(inputDir, metaFile)
myTable<-read.table(file.path,sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
# myTable$time <- as.factor(myTable$time)

Timepoint <- myTable$time
PatientID <- myTable$PatientID
Weight <- myTable$Weight_kg
BMI <- myTable$BMI_kgm2
Surgery <- myTable$Surgery
PEWL <- myTable$PEWL_kg
PWL <- myTable$Percent_Loss_kg


##### Prep for modeling #####
VariableName <- vector()
p_Timepoint <- vector()
p_Weight_kg <- vector()
p_Surgery <- vector()
p_PEWL <- vector()
p_PWL <- vector()
index <- 1

##### Multivariate mixed linear model #####
for (variable in metaVariables) {
  
  message("\n***** Starting ", variable, " *****\n")
  VariableName[index] <- variable
  metaCol <- myTable[,which(colnames(myTable) == variable)]
  
  
  df<-data.frame(PatientID, Timepoint, Weight, PEWL, PWL, Surgery, metaCol)
  df <- na.omit(df)
  df$Timepoint <- factor(df$Timepoint)

  fit <- anova(lme(metaCol ~ Weight + Timepoint + Surgery, method="REML", random=~1 | PatientID, data=df))
  p_Weight_kg[index] <- fit$`p-value`[2]
  p_Timepoint[index] <- fit$`p-value`[3]
  p_Surgery[index] <- fit$`p-value`[4]
  
  df2 <- df[!(df$Timepoint == 0),]
  fit <- anova(lme(metaCol ~ PEWL + Timepoint + Surgery, method="REML", random=~1 | PatientID, data=df2))
  p_PEWL[index] <- fit$`p-value`[2]
  
  fit <- anova(lme(metaCol ~ PWL + Timepoint + Surgery, method="REML", random=~1 | PatientID, data=df2))
  p_PWL[index] <- fit$`p-value`[2]
  
  index <- index + 1

} # for (variable in metaVariables)

dFrame <- data.frame(VariableName, p_Weight_kg, p_Timepoint, p_Surgery, p_PEWL, p_PWL)
dFrame$pAdj_Weight_kg <- p.adjust(dFrame$p_Weight_kg, method = "BH")
dFrame$pAdj_Timepoint <- p.adjust(dFrame$p_Timepoint, method = "BH")
dFrame$pAdj_Surgery <- p.adjust(dFrame$p_Surgery, method = "BH")
dFrame$pAdj_PEWL <- p.adjust(dFrame$p_PEWL, method = "BH")
dFrame$pAdj_PWL <- p.adjust(dFrame$p_PWL, method = "BH")

##### Plots from modeling #####
file.path <- paste0(outputDir,"MLM_BLto", included[length(included)], "_plots.pdf")
pdf(file.path)

par(mfrow=c(2,2))
tFrame <- data.frame()
for (row in 1:nrow(dFrame)) {
  
  
  Variable <- dFrame$VariableName[row]
  metaCol <- which(colnames(myTable) == Variable)
  
  if (dFrame$pAdj_Weight_kg[row] < 0.05) {
    color <- "red"
  } else {
    color <- "black"
  }
  
  title.lab <- paste0( Variable, "\nWeight (adj p = ", format(dFrame$pAdj_Weight_kg[row], length = 3), ")")
  plot(  myTable[, metaCol], myTable$Weight_kg,
         xlab=Variable,
         ylab ="Weight (kg)",
         main = title.lab,
         cex.main = 1,
         col.main = color
  )

  myTable2 <- myTable[!(myTable$time == 0),]

  if (dFrame$pAdj_PEWL[row] < 0.05) {
    color <- "red"
  } else {
    color <- "black"
  }
  
  title.lab <- paste0( Variable , "\n%EWL (adj p = ", format(dFrame$pAdj_PEWL[row], length = 3), ")")
  plot(  myTable2[, metaCol], myTable2$PEWL_kg,
         xlab=Variable,
         ylab ="%EWL", 
         main = title.lab,
         cex.main = 1,
         col.main = color
  )
  
  if (dFrame$pAdj_PWL[row] < 0.05) {
    color <- "red"
  } else {
    color <- "black"
  }
  title.lab <- paste0( Variable , "\n%WL (adj p = ", format(dFrame$pAdj_PWL[row], length = 3), ")")
  plot(  myTable2[, metaCol], myTable2$Percent_Loss_kg,
         xlab=Variable,
         ylab ="Weight loss (%)", 
         main = title.lab,
         cex.main = 1,
         col.main = color
  )
  
  if (dFrame$pAdj_Timepoint[row] < 0.05) {
    color <- "red"
  } else {
    color <- "black"
  }
  title.lab <- paste0( Variable , "\nTime (adj p = ", format(dFrame$pAdj_Timepoint[row], length = 3), ")")
  plot(myTable$time, myTable[, metaCol],
         xlab="Time (months)",
         ylab =Variable,
         main = title.lab,
         cex.main = 1,
       col.main = color
  )
  # legend("topright", legend=c("SG", "RYGB"),
  #        col=c("red", "blue"), pch = 21, cex=0.8)
  
  
  if (dFrame$pAdj_Surgery[row] < 0.05) {
    color <- "red"
  } else {
    color <- "black"
  }
  title.lab <- paste0( Variable , "\nSurgery (adj p = ", format(dFrame$pAdj_Surgery[row], length = 3), ")")
  boxplot( myTable[ ,metaCol] ~ myTable$Surgery, 
       xlab="Surgery",
       ylab =Variable,
       main = title.lab,
       cex.main = 1,
       col.main = color
  )
  
  Weight_kg_pValues <- vector()
  timeNumbers <- vector()
  tvarName <- vector()
  tIndex <- 1
  
  for (j in included) {
    
    varCol <- myTable[, metaCol]
    timeCol <- myTable$time
    wghtCol <- myTable$Weight_kg
    id <- myTable$PatientID
    
    lm.df <- data.frame(id, varCol, timeCol, wghtCol)
    lm.df <- na.omit(lm.df)
    lm.df <- lm.df[lm.df$timeCol == j,]
    
    myLm <- lm( lm.df$varCol ~ lm.df$wghtCol )
    p <- anova(myLm)$"Pr(>F)"[1]
    Weight_kg_pValues[tIndex] <- p
    timeNumbers[tIndex] <- j
    tvarName[tIndex] <- Variable
    tIndex <- tIndex + 1
    
    if (p < 0.05) {
      color <- "red"
    } else {
      color <- "black"
    }
    title.lab <- paste0( month.labs[which(included == j)], " (", length(unique(lm.df$id)), " patients)\n", Variable , "\nWeight (p = ", format(p, length = 3), ")")
    plot(  lm.df$wghtCol, lm.df$varCol ,
           ylab =Variable , 
           xlab="Weight (kg)", 
           main=title.lab,
           cex.main = 1,
           col.main = color
    )
    
  }
  tempFrame <- data.frame(tvarName, timeNumbers, Weight_kg_pValues)
  tempFrame$adj_Weight_kg_pValues <- p.adjust(tempFrame$Weight_kg_pValues, method="BH")


  
  PEWL_pValues <- vector()
  tIndex <- 1
  
  for (j in included[2:length(included)]) {
    
    varCol <- myTable2[, metaCol]
    timeCol <- myTable2$time
    wghtCol <- myTable2$PEWL_kg
    id <- myTable2$PatientID
    
    lm.df <- data.frame(id, varCol, timeCol, wghtCol)
    lm.df <- na.omit(lm.df)
    lm.df <- lm.df[lm.df$timeCol == j,]
    
    myLm <- lm( lm.df$varCol ~ lm.df$wghtCol )
    p <- anova(myLm)$"Pr(>F)"[1]
    PEWL_pValues[tIndex] <- p
    tIndex <- tIndex + 1
    
    if (p < 0.05) {
      color <- "red"
    } else {
      color <- "black"
    }
    title.lab <- paste0( month.labs[which(included == j)], " (", length(unique(lm.df$id)), " patients)\n", Variable , "\n%EWL (p = ", format(p, length = 3), ")")
    plot(  lm.df$wghtCol, lm.df$varCol ,
           ylab =Variable , 
           xlab="%EWL", 
           main=title.lab,
           cex.main = 1,
           col.main = color
    )
    
  }
  PEWL_pValues <- c(NA, PEWL_pValues)
  tempFrame$PEWL_pValues <-PEWL_pValues
  tempFrame$adj_PEWL_pValues <- p.adjust(tempFrame$PEWL_pValues, method="BH")
  
  
  PWL_pValues <- vector()
  tIndex <- 1
  
  for (j in included[2:length(included)]) {
    
    varCol <- myTable2[, metaCol]
    timeCol <- myTable2$time
    wghtCol <- myTable2$Percent_Loss_kg
    id <- myTable2$PatientID
    
    lm.df <- data.frame(id, varCol, timeCol, wghtCol)
    lm.df <- na.omit(lm.df)
    lm.df <- lm.df[lm.df$timeCol == j,]
    
    myLm <- lm( lm.df$varCol ~ lm.df$wghtCol )
    p <- anova(myLm)$"Pr(>F)"[1]
    PWL_pValues[tIndex] <- p
    tIndex <- tIndex + 1
    
    if (p < 0.05) {
      color <- "red"
    } else {
      color <- "black"
    }
    title.lab <- paste0( month.labs[which(included == j)], " (", length(unique(lm.df$id)), " patients)\n", Variable , "\n%WL (p = ", format(p, length = 3), ")")
    plot(  lm.df$wghtCol, lm.df$varCol ,
           ylab =Variable , 
           xlab="%WL", 
           main=title.lab,
           cex.main = 1,
           col.main = color
    )
    
  }
  PWL_pValues <- c(NA, PWL_pValues)
  tempFrame$PWL_pValues <-PWL_pValues
  tempFrame$adj_PWL_pValues <- p.adjust(tempFrame$PWL_pValues, method="BH")
  
  tFrame <- rbind(tFrame, tempFrame)
  plot.new()
  
} # for (row in 1:nrow(dFrame))

dev.off()


file.path <- paste0(outputDir, "Weight_MLM_pValues_at_each_timepoint_BLto", included[length(included)], ".tsv")
write.table(tFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)



file.path <- paste0(outputDir, statsOutput)
write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Univariate linear model #####
PatientID <- myTable$PatientID
Timepoint <- myTable$time
metaName <- vector()
pValue_Timepoint.MLM <- vector()
pValue_Timepoint.LM <- vector()
p_BLv1M.LM <- vector()
p_BLv6M.LM <- vector()
p_BLv12M.LM <- vector()
p_BLv18M.LM <- vector()
p_BLv24M.LM <- vector()

index <- 1

for (i in metaVariables){
  
  meta <- myTable[,which(colnames(myTable) == i)]
  df<-data.frame(meta,PatientID, Timepoint)
  df <- df[df$Timepoint %in% included,]
  df <- na.omit(df)
  df$Timepoint <- as.factor(df$Timepoint)
  
  mlm <- lme(meta ~ Timepoint, method="REML", random=~1 | PatientID, data=df)
  smry <- summary(mlm)
  fit<-anova(mlm)
  pValue_Timepoint.MLM[index]<-fit$`p-value`[2]
  
  lm <- lm(df$meta ~ df$Timepoint)
  smry <- summary(lm)
  fit<-anova(lm)
  pValue_Timepoint.LM[index]<-fit$`Pr(>F)`[1]
  p_BLv1M.LM[index] <- smry$coefficients[2,4]
  if (length(included) > 2) {  p_BLv6M.LM[index] <- smry$coefficients[3,4]}
  if (length(included) > 3) {  p_BLv12M.LM[index] <- smry$coefficients[4,4]}
  if (length(included) > 4) {  p_BLv18M.LM[index] <- smry$coefficients[5,4]}
  if (length(included) > 5) {  p_BLv24M.LM[index] <- smry$coefficients[6,4]}
  metaName[index]<-i
  
  index<-index+1
  
  
} #for (i in 1:ncol(myT))


taxa.df<-data.frame(metaName, pValue_Timepoint.MLM, pValue_Timepoint.LM, p_BLv1M.LM)

taxa.df$Adj_pValue_Timepoint.MLM <- p.adjust(taxa.df$pValue_Timepoint.MLM,method = "BH")
taxa.df$Adj_pValue_Timepoint.MLM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.MLM)

taxa.df$Adj_pValue_Timepoint.LM <- p.adjust(taxa.df$pValue_Timepoint.LM,method = "BH")
taxa.df$Adj_pValue_Timepoint.LM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.LM)

taxa.df$Adj_p_BLv1M.LM<- p.adjust(taxa.df$p_BLv1M.LM,method = "BH")
taxa.df$Adj_p_BLv1M.LM_sig <- sigStars(taxa.df$Adj_p_BLv1M.LM)

if (length(included) > 2) {
  taxa.df$p_BLv6M.LM <- p_BLv6M.LM
  taxa.df$Adj_p_BLv6M.LM <- p.adjust(taxa.df$p_BLv6M.LM,method = "BH")
  taxa.df$Adj_p_BLv6M.LM_sig <- sigStars(taxa.df$Adj_p_BLv6M.LM)
}

if (length(included) > 3) {
  taxa.df$p_BLv12M.LM <- p_BLv12M.LM
  taxa.df$Adj_p_BLv12M.LM <- p.adjust(taxa.df$p_BLv12M.LM,method = "BH")
  taxa.df$Adj_p_BLv12M.LM_sig <- sigStars(taxa.df$Adj_p_BLv12M.LM)
}

if (length(included) > 4) {
  taxa.df$p_BLv18M.LM <- p_BLv18M.LM
  taxa.df$Adj_p_BLv18M.LM <- p.adjust(taxa.df$p_BLv18M.LM,method = "BH")
  taxa.df$Adj_p_BLv18M.LM_sig <- sigStars(taxa.df$Adj_p_BLv18M.LM)
}

if (length(included) > 5) {
  taxa.df$p_BLv24M.LM <- p_BLv24M.LM
  taxa.df$Adj_p_BLv24M.LM <- p.adjust(taxa.df$p_BLv24M.LM,method = "BH")
  taxa.df$Adj_p_BLv24M.LM_sig <- sigStars(taxa.df$Adj_p_BLv24M.LM)
}

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months.tsv")
write.table(taxa.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### Univariate linear model plots + pairwise t-test: Nutrient intake over time #####
myTable2 <- myTable[myTable$time %in% included,]

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months.tsv")
pVal.df <- read.table(file.path, sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
pVal.df$pFinal <- roundP(pVal.df$Adj_pValue_Timepoint.MLM)

plotList <- list()
finalStats <- data.frame()

for (row in 1:nrow(pVal.df)) {
  
  metaName <- rownames(pVal.df[row,])
  
  metaCol <- myTable2[,which(colnames(myTable2) == metaName)]
  month <- myTable2$time
  df <- data.frame(metaCol, month)
  df <- na.omit(df)
  
  df$metaCol <- log10(df$metaCol + 1)
  
  stat.test <- df %>%
    t_test(metaCol ~ month, ref.group = "0") %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  stat.test$.y. <- metaName
  finalStats <- rbind(finalStats, stat.test)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "month", fun = "max")
  
  caption.lab <- paste0("pairwise adj p: t_test(", metaName, " ~ Timepoint)")
  subtitle.lab <- paste0("anova(lme(metaName ~ Timepoint, method=REML, random=~1 | PatientID))")
  
  plot <- ggboxplot(
    df, x = "month", y = "metaCol",
    # fill = "#2EDF52",
    # facet.by = c("Surgery"),
    scales = "free", add = "jitter"
  )
  
  plot <- plot +
    labs(x = "Post-op Time (months)", y = xLabels[which(metaVariables == metaName)])
  
  plot <- plot +
    labs(title = paste0(metaName, " (adj ", pVal.df$pFinal[row], ")"),
         caption = caption.lab,
         subtitle = subtitle.lab)
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 6,
      hide.ns = TRUE,
      label = "p.adj.signif",
      step.increase = 0.07
    )
  # plot <- plot + scale_y_log10()
  plot
  
  plotList[[row]] <- plot
  
} # for (row in 1:nrow(pVal.df))

file.path <- paste0(outputDir, "ttestResults_BLto", included[length(included)], "months.tsv")
write.table(finalStats, file.path,sep="\t",quote = FALSE, row.names = FALSE)


# Output plot
file.path <- paste0(outputDir, "UnivariateMLM_Timepoint_BLto", included[length(included)], "months_Boxplots_ttest.pdf")
pdf(file.path, width = 7, height = 4)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()




##### Univariate linear model plots + pairwise Wilcoxon: Nutrient intake over time #####
myTable2 <- myTable[myTable$time %in% included,]

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months.tsv")
pVal.df <- read.table(file.path, sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
pVal.df$pFinal <- roundP(pVal.df$Adj_pValue_Timepoint.MLM)

plotList <- list()
finalStats <- data.frame()

for (row in 1:nrow(pVal.df)) {
  
  metaName <- rownames(pVal.df[row,])
  
  metaCol <- myTable2[,which(colnames(myTable2) == metaName)]
  month <- myTable2$time
  df <- data.frame(metaCol, month)
  df <- na.omit(df)
  
  df$metaCol <- log10(df$metaCol + 1)
  
  stat.test <- df %>%
    wilcox_test(metaCol ~ month, ref.group = "0") %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  stat.test$.y. <- metaName
  finalStats <- rbind(finalStats, stat.test)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "month", fun = "max")
  
  caption.lab <- paste0("pairwise adj p: wilcox_test(", metaName, " ~ Timepoint)")
  subtitle.lab <- paste0("anova(lme(metaName ~ Timepoint, method=REML, random=~1 | PatientID))")
  
  plot <- ggboxplot(
    df, x = "month", y = "metaCol",
    # fill = "#2EDF52",
    # facet.by = c("Surgery"),
    scales = "free", add = "jitter"
  )
  
  plot <- plot +
    labs(x = "Post-op Time (months)", y = xLabels[which(metaVariables == metaName)])
  
  plot <- plot +
    labs(title = paste0(metaName, " (adj ", pVal.df$pFinal[row], ")"),
         caption = caption.lab,
         subtitle = subtitle.lab)
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 6,
      hide.ns = TRUE,
      label = "p.adj.signif",
      step.increase = 0.07
    )
  # plot <- plot + scale_y_log10()
  plot
  
  plotList[[row]] <- plot
  
} # for (row in 1:nrow(pVal.df))

file.path <- paste0(outputDir, "WilcoxonResults_BLto", included[length(included)], "months.tsv")
write.table(finalStats, file.path,sep="\t",quote = FALSE, row.names = FALSE)


# Output plot
file.path <- paste0(outputDir, "UnivariateMLM_Timepoint_BLto", included[length(included)], "months_Boxplots_wilcox.pdf")
pdf(file.path, width = 7, height = 4)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()










##### Univariate linear model: Nutrient intake over time by surgery type #####
PatientID <- myTable$PatientID
Timepoint <- myTable$time
Surgery <- myTable$Surgery

metaName <- vector()
pValue_Timepoint.MLM <- vector()
pValue_Timepoint.LM <- vector()
p_BLv1M.LM <- vector()
p_BLv6M.LM <- vector()
p_BLv12M.LM <- vector()
p_BLv18M.LM <- vector()
p_BLv24M.LM <- vector()
SurgeryType <- vector()

index <- 1

for (i in metaVariables){
  
  meta <- myTable[,which(colnames(myTable) == i)]
  df<-data.frame(meta,PatientID, Timepoint, Surgery)
  df <- df[df$Timepoint %in% included,]
  df <- na.omit(df)
  df$Timepoint <- as.factor(df$Timepoint)
  
  for (type in unique(df$Surgery)) {
    
    df2 <- df[df$Surgery == type,]
    SurgeryType[index]<-type
    
    mlm <- lme(meta ~ Timepoint, method="REML", random=~1 | PatientID, data=df2)
    smry <- summary(mlm)
    fit<-anova(mlm)
    pValue_Timepoint.MLM[index]<-fit$`p-value`[2]
    
    lm <- lm(df2$meta ~ df2$Timepoint)
    smry <- summary(lm)
    fit<-anova(lm)
    pValue_Timepoint.LM[index]<-fit$`Pr(>F)`[1]
    p_BLv1M.LM[index] <- smry$coefficients[2,4]
    if (length(included) > 2) {  p_BLv6M.LM[index] <- smry$coefficients[3,4]}
    if (length(included) > 3) {  p_BLv12M.LM[index] <- smry$coefficients[4,4]}
    if (length(included) > 4) {  p_BLv18M.LM[index] <- smry$coefficients[5,4]}
    if (length(included) > 5) {  p_BLv24M.LM[index] <- smry$coefficients[6,4]}
    metaName[index]<-i
    
    index<-index+1
    
  }
  
} #for (i in 1:ncol(myT))


taxa.df<-data.frame(metaName, SurgeryType, pValue_Timepoint.MLM, pValue_Timepoint.LM, p_BLv1M.LM)

taxa.df$Adj_pValue_Timepoint.MLM <- p.adjust(taxa.df$pValue_Timepoint.MLM,method = "BH")
taxa.df$Adj_pValue_Timepoint.MLM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.MLM)

taxa.df$Adj_pValue_Timepoint.LM <- p.adjust(taxa.df$pValue_Timepoint.LM,method = "BH")
taxa.df$Adj_pValue_Timepoint.LM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.LM)

taxa.df$Adj_p_BLv1M.LM<- p.adjust(taxa.df$p_BLv1M.LM,method = "BH")
taxa.df$Adj_p_BLv1M.LM_sig <- sigStars(taxa.df$Adj_p_BLv1M.LM)

if (length(included) > 2) {
  taxa.df$p_BLv6M.LM <- p_BLv6M.LM
  taxa.df$Adj_p_BLv6M.LM <- p.adjust(taxa.df$p_BLv6M.LM,method = "BH")
  taxa.df$Adj_p_BLv6M.LM_sig <- sigStars(taxa.df$Adj_p_BLv6M.LM)
}

if (length(included) > 3) {
  taxa.df$p_BLv12M.LM <- p_BLv12M.LM
  taxa.df$Adj_p_BLv12M.LM <- p.adjust(taxa.df$p_BLv12M.LM,method = "BH")
  taxa.df$Adj_p_BLv12M.LM_sig <- sigStars(taxa.df$Adj_p_BLv12M.LM)
}

if (length(included) > 4) {
  taxa.df$p_BLv18M.LM <- p_BLv18M.LM
  taxa.df$Adj_p_BLv18M.LM <- p.adjust(taxa.df$p_BLv18M.LM,method = "BH")
  taxa.df$Adj_p_BLv18M.LM_sig <- sigStars(taxa.df$Adj_p_BLv18M.LM)
}

if (length(included) > 5) {
  taxa.df$p_BLv24M.LM <- p_BLv24M.LM
  taxa.df$Adj_p_BLv24M.LM <- p.adjust(taxa.df$p_BLv24M.LM,method = "BH")
  taxa.df$Adj_p_BLv24M.LM_sig <- sigStars(taxa.df$Adj_p_BLv24M.LM)
}

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months_by_SurgeryType.tsv")
write.table(taxa.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### Univariate linear model plots + pairwise t-test: Nutrient intake over time by surgery type #####
myTable2 <- myTable[myTable$time %in% included,]

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months_by_SurgeryType.tsv")
pVal.df <- read.table(file.path, sep="\t",header = TRUE,check.names = FALSE)
pVal.df$pFinal <- roundP(pVal.df$Adj_pValue_Timepoint.MLM)

plotList <- list()
finalStats <- data.frame()
index <- 1

# Filter by variable 
for (metaVariable in metaVariables) {
  
  pVal.df2 <- pVal.df[pVal.df$metaName == metaVariable,]
  
  # Filter by surgery type
  for (surgery in c("RYGB", "SG")) {
    
    pVal.df3 <- pVal.df2[pVal.df2$SurgeryType == surgery,]
    
    for (row in 1:nrow(pVal.df3)) {
      
      metaName <- pVal.df3$metaName[row]
      
      metaCol <- myTable2[,which(colnames(myTable2) == metaName)]
      month <- myTable2$time
      df <- data.frame(metaCol, month)
      df <- na.omit(df)
      
      df$metaCol <- log10(df$metaCol + 1)
      
      stat.test <- df %>%
        t_test(metaCol ~ month, ref.group = "0") %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stat.test$.y. <- metaName
      finalStats <- rbind(finalStats, stat.test)
      
      stat.test <- stat.test %>%
        add_xy_position(x = "month", fun = "max")
      
      caption.lab <- paste0("pairwise adj p: t_test(", metaName, " ~ Timepoint)")
      subtitle.lab <- paste0("lme(metaName ~ Timepoint, method=REML, random=~1 | PatientID)")
      
      plot <- ggboxplot(
        df, x = "month", y = "metaCol",
        # fill = "#2EDF52",
        # facet.by = c("Surgery"),
        scales = "free", add = "jitter"
      )
      
      plot <- plot +
        labs(x = "Post-op Time (months)", y = xLabels[which(metaVariables == metaName)])
      
      plot <- plot +
        labs(title = paste0(surgery, " (adj ", pVal.df3$pFinal[row], ")"),
             caption = caption.lab,
             subtitle = subtitle.lab)
      
      plot <- plot +
        theme(plot.subtitle = element_text(size = 8))
      
      plot <- plot +
        stat_pvalue_manual(
          stat.test,
          bracket.nudge.y = 0.01,
          size = 6,
          hide.ns = TRUE,
          label = "p.adj.signif",
          step.increase = 0.07
        )
      # plot <- plot + scale_y_log10()
      plot
      
      plotList[[index]] <- plot
      index <- index + 1
    } # for (row in 1:nrow(pVal.df3))
    
  } # for (surgery in c("RYGB", "SG"))
  
} # for (metaVariable in metaVariables)
file.path <- paste0(outputDir, "ttestResults_BLto", included[length(included)], "months_by_Surgery.tsv")
write.table(finalStats, file.path,sep="\t",quote = FALSE, row.names = FALSE)


# Output plot
file.path <- paste0(outputDir, "UnivariateMLM_Timepoint_BLto", included[length(included)], "months_Boxplots_ttest_by_Surgery.pdf")
pdf(file.path, width = 10, height = 4)
for (i in seq(1, length(plotList), 2)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               ncol = 2, nrow = 1)
}
dev.off()






##### Univariate linear model: Nutrient intake over time by responder status #####
SampleID <- rownames(myTable)
myTable <- cbind(SampleID, myTable)
months <- c("12", "18", "24")
included <- c("0", "1", "6")
# index <- 1

for (month in months) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  KCAL <- myTable$KCAL
  CARB <- myTable$CARB
  PROT <- myTable$PROT
  TFAT <- myTable$TFAT
  PEWL_kg <- myTable$PEWL_kg
  
  outcome.df <- data.frame(PatientID, Timepoint, PEWL_kg)
  outcome.df <- na.omit(outcome.df)
  outcome.df$Timepoint <- paste0("M", outcome.df$Timepoint)
  
  outcome.df <- spread(outcome.df, Timepoint, PEWL_kg)
  m <- which(colnames(outcome.df) == paste0("M", month))
  
  outcome.df$ResponderStatus <- ifelse(outcome.df[,m] >= 50, "Responder",
                                       ifelse(outcome.df[,m] < 50, "Non-responder",
                                              NA))
  outcome.df <- outcome.df[!is.na(outcome.df$ResponderStatus),]
  bl <- which(colnames(outcome.df) == "M0")
  end <- which(colnames(outcome.df) == "ResponderStatus")-1
  outcome.df <- outcome.df %>%
    gather("Timepoint", "PEWL_kg", bl:end)
  outcome.df$Timepoint <- gsub("M", "", outcome.df$Timepoint)
  
  included <- c(included, month)
  
  outcome.df <- outcome.df[outcome.df$Timepoint %in% included,]
  outcome.df$Timepoint <- as.factor(outcome.df$Timepoint)
  outcome.df$Timepoint <- factor(outcome.df$Timepoint, levels = included)
  outcome.df$SampleID <- paste0(outcome.df$PatientID, "-", str_pad(outcome.df$Timepoint, width=2, pad="0"))
  
  names(KCAL) = myTable[, "SampleID"]
  names(CARB) = myTable[, "SampleID"]
  names(PROT) = myTable[, "SampleID"]
  names(TFAT) = myTable[, "SampleID"]
  
  # Add info to each row
  outcome.df$KCAL = KCAL[outcome.df$SampleID]
  outcome.df$CARB = CARB[outcome.df$SampleID]
  outcome.df$PROT = PROT[outcome.df$SampleID]
  outcome.df$TFAT = PROT[outcome.df$SampleID]
  
  # Set up for analysis
  PatientID <- outcome.df$PatientID
  Timepoint <- outcome.df$Timepoint
  ResponderStatus <- outcome.df$ResponderStatus
  
  metaName <- vector()
  pValue_Timepoint.MLM <- vector()
  pValue_Timepoint.LM <- vector()
  p_BLv1M.LM <- vector()
  p_BLv6M.LM <- vector()
  p_BLv12M.LM <- vector()
  p_BLv18M.LM <- vector()
  p_BLv24M.LM <- vector()
  Outcome <- vector()
  
  index <- 1
  
  for (i in metaVariables){
    
    meta <- outcome.df[,which(colnames(outcome.df) == i)]
    df<-data.frame(meta,PatientID, Timepoint, ResponderStatus)
    df <- df[df$Timepoint %in% included,]
    df <- na.omit(df)
    df$Timepoint <- as.factor(df$Timepoint)
    
    for (type in unique(df$ResponderStatus)) {
      
      df2 <- df[df$ResponderStatus == type,]
      Outcome[index]<-type
      
      mlm <- lme(meta ~ Timepoint, method="REML", random=~1 | PatientID, data=df2)
      smry <- summary(mlm)
      fit<-anova(mlm)
      pValue_Timepoint.MLM[index]<-fit$`p-value`[2]
      
      lm <- lm(df2$meta ~ df2$Timepoint)
      smry <- summary(lm)
      fit<-anova(lm)
      pValue_Timepoint.LM[index]<-fit$`Pr(>F)`[1]
      p_BLv1M.LM[index] <- smry$coefficients[2,4]
      if (length(included) > 2) {  p_BLv6M.LM[index] <- smry$coefficients[3,4]}
      if (length(included) > 3) {  p_BLv12M.LM[index] <- smry$coefficients[4,4]}
      if (length(included) > 4) {  p_BLv18M.LM[index] <- smry$coefficients[5,4]}
      if (length(included) > 5) {  p_BLv24M.LM[index] <- smry$coefficients[6,4]}
      metaName[index]<-i
      
      index<-index+1
      
    }
    
  } #for (i in 1:ncol(myT))
  
  
  taxa.df<-data.frame(metaName, Outcome, pValue_Timepoint.MLM, pValue_Timepoint.LM, p_BLv1M.LM)
  
  taxa.df$Adj_pValue_Timepoint.MLM <- p.adjust(taxa.df$pValue_Timepoint.MLM,method = "BH")
  taxa.df$Adj_pValue_Timepoint.MLM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.MLM)
  
  taxa.df$Adj_pValue_Timepoint.LM <- p.adjust(taxa.df$pValue_Timepoint.LM,method = "BH")
  taxa.df$Adj_pValue_Timepoint.LM_sig <- sigStars(taxa.df$Adj_pValue_Timepoint.LM)
  
  taxa.df$Adj_p_BLv1M.LM<- p.adjust(taxa.df$p_BLv1M.LM,method = "BH")
  taxa.df$Adj_p_BLv1M.LM_sig <- sigStars(taxa.df$Adj_p_BLv1M.LM)
  
  if (length(included) > 2) {
    taxa.df$p_BLv6M.LM <- p_BLv6M.LM
    taxa.df$Adj_p_BLv6M.LM <- p.adjust(taxa.df$p_BLv6M.LM,method = "BH")
    taxa.df$Adj_p_BLv6M.LM_sig <- sigStars(taxa.df$Adj_p_BLv6M.LM)
  }
  
  if (length(included) > 3) {
    taxa.df$p_BLv12M.LM <- p_BLv12M.LM
    taxa.df$Adj_p_BLv12M.LM <- p.adjust(taxa.df$p_BLv12M.LM,method = "BH")
    taxa.df$Adj_p_BLv12M.LM_sig <- sigStars(taxa.df$Adj_p_BLv12M.LM)
  }
  
  if (length(included) > 4) {
    taxa.df$p_BLv18M.LM <- p_BLv18M.LM
    taxa.df$Adj_p_BLv18M.LM <- p.adjust(taxa.df$p_BLv18M.LM,method = "BH")
    taxa.df$Adj_p_BLv18M.LM_sig <- sigStars(taxa.df$Adj_p_BLv18M.LM)
  }
  
  if (length(included) > 5) {
    taxa.df$p_BLv24M.LM <- p_BLv24M.LM
    taxa.df$Adj_p_BLv24M.LM <- p.adjust(taxa.df$p_BLv24M.LM,method = "BH")
    taxa.df$Adj_p_BLv24M.LM_sig <- sigStars(taxa.df$Adj_p_BLv24M.LM)
  }
  
  file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months_by_", month, "M_Outcome.tsv")
  write.table(taxa.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (month in months)

##### Univariate linear model plots + pairwise t-test: Nutrient intake over time by responder status #####
months <- c("12", "18", "24")
included <- c("0", "1", "6")
# index <- 1

for (month in months) {
  
  included <- c(included, month)
  
  myTable2 <- myTable[myTable$time %in% included,]
  
  file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months_by_", month, "M_Outcome.tsv")
  pVal.df <- read.table(file.path, sep="\t",header = TRUE,check.names = FALSE)
  pVal.df$pFinal <- roundP(pVal.df$Adj_pValue_Timepoint.MLM)
  
  plotList <- list()
  finalStats <- data.frame()
  index <- 1
  
  # Filter by variable 
  for (metaVariable in metaVariables) {
    
    pVal.df2 <- pVal.df[pVal.df$metaName == metaVariable,]
    
    # Filter by responder status
    for (outcome in c("Responder", "Non-responder")) {
      
      pVal.df3 <- pVal.df2[pVal.df2$Outcome == outcome,]
      
      for (row in 1:nrow(pVal.df3)) {
        
        metaName <- pVal.df3$metaName[row]
        
        metaCol <- myTable2[,which(colnames(myTable2) == metaName)]
        Month <- myTable2$time
        df <- data.frame(metaCol, Month)
        df <- na.omit(df)
        
        df$metaCol <- log10(df$metaCol + 1)
        
        stat.test <- df %>%
          t_test(metaCol ~ Month, ref.group = "0") %>%
          adjust_pvalue(method = "BH") %>%
          add_significance()
        stat.test$.y. <- metaName
        finalStats <- rbind(finalStats, stat.test)
        
        stat.test <- stat.test %>%
          add_xy_position(x = "Month", fun = "max")
        
        caption.lab <- paste0("pairwise adj p: t_test(", metaName, " ~ Timepoint)")
        subtitle.lab <- paste0("lme(metaName ~ Timepoint, method=REML, random=~1 | PatientID)")
        
        plot <- ggboxplot(
          df, x = "Month", y = "metaCol",
          # fill = "#2EDF52",
          # facet.by = c("Outcome"),
          scales = "free", add = "jitter"
        )
        
        plot <- plot +
          labs(x = "Post-op Time (months)", y = xLabels[which(metaVariables == metaName)])
        
        plot <- plot +
          labs(title = paste0(outcome, "s (adj ", pVal.df3$pFinal[row], ")"),
               caption = caption.lab,
               subtitle = subtitle.lab)
        
        plot <- plot +
          theme(plot.subtitle = element_text(size = 8))
        
        plot <- plot +
          stat_pvalue_manual(
            stat.test,
            bracket.nudge.y = 0.01,
            size = 6,
            hide.ns = TRUE,
            label = "p.adj.signif",
            step.increase = 0.07
          )
        # plot <- plot + scale_y_log10()
        plot
        
        plotList[[index]] <- plot
        index <- index + 1
      } # for (row in 1:nrow(pVal.df3))
      
    } # for (outcome in c("RYGB", "SG"))
    
  } # for (metaVariable in metaVariables)
  file.path <- paste0(outputDir, "ttestResults_BLto", included[length(included)], "months_by_", month, "M_Outcome.tsv")
  write.table(finalStats, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
  
  # Output plot
  file.path <- paste0(outputDir, "UnivariateMLM_Timepoint_BLto", included[length(included)], "months_Boxplots_ttest_by_", month, "M_Outcome.pdf")
  pdf(file.path, width = 10, height = 4)
  for (i in seq(1, length(plotList), 2)) {
    grid.arrange(plotList[[i]], plotList[[i+1]],
                 ncol = 2, nrow = 1)
  }
  dev.off()
  
} # for (month in months)


##### Univariate linear model plots + pairwise Wilcoxon: Nutrient intake over time PUBLICATION READY #####
myTable2 <- myTable[myTable$time %in% included,]

file.path <- paste0(outputDir, "UnivariateMLMResults_BLto", included[length(included)], "months.tsv")
pVal.df <- read.table(file.path, sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
pVal.df$pFinal <- roundP(pVal.df$Adj_pValue_Timepoint.MLM)

plotList <- list()
finalStats <- data.frame()

for (row in 1:nrow(pVal.df)) {
  
  metaName <- rownames(pVal.df[row,])
  
  metaCol <- myTable2[,which(colnames(myTable2) == metaName)]
  month <- myTable2$time
  df <- data.frame(metaCol, month)
  df <- na.omit(df)
  
  df$metaCol <- log10(df$metaCol + 1)
  
  stat.test <- df %>%
    wilcox_test(metaCol ~ month, ref.group = "0") %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  stat.test$.y. <- metaName
  finalStats <- rbind(finalStats, stat.test)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "month", fun = "max")
  
  caption.lab <- paste0("pairwise adj p: wilcox_test(", metaName, " ~ Timepoint)")
  subtitle.lab <- paste0("anova(lme(metaName ~ Timepoint, method=REML, random=~1 | PatientID))")
  
  plot <- ggboxplot(
    df, x = "month", y = "metaCol",
    # fill = "#2EDF52",
    # facet.by = c("Surgery"),
    scales = "free", add = "jitter"
  )
  
  tag <- tags[row]
  plot <- plot +
    labs(x = "Post-op Time (months)", y = xLabels[which(metaVariables == metaName)])
  
  plot <- plot +
    labs(tag = tag)
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 6,
      hide.ns = TRUE,
      label = "p.adj.signif",
      step.increase = 0.07
    )
  # plot <- plot + scale_y_log10()
  plot
  
  plotList[[row]] <- plot
  
} # for (row in 1:nrow(pVal.df))

file.path <- paste0(outputDir, "WilcoxonResults_BLto", included[length(included)], "months.tsv")
write.table(finalStats, file.path,sep="\t",quote = FALSE, row.names = FALSE)


# Output plot
file.path <- paste0(outputDir, "UnivariateMLM_Timepoint_BLto", included[length(included)], "months_Boxplots_wilcox_PUB.pdf")
pdf(file.path, width = 7.5, height = 7.5)
for (i in seq(1, length(plotList), 4)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], plotList[[i+3]],
               ncol = 2, nrow = 2)
}
dev.off()
