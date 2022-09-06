
#This function parses out concatenated taxa names.
parseQIIME2Taxa<-function(table){
  
  bugNameSplit=gsub(pattern =  "d__", x = table$bugName,replacement = "")
  bugNameSplit=gsub(pattern =  "[.].__", x = bugNameSplit,replacement = "_/_")
  other <- paste(replicate(count,"Other"), collapse = "_/_")
  bugNameSplit=gsub(pattern ="^Other$",x = bugNameSplit,replacement = other)
  string <- strsplit(as.character(bugNameSplit),split = "_/_")
  
  temp_string=do.call(rbind,string)
  table <- cbind(temp_string,table)
  
  colnames(table)[colnames(table)=="1"] <- "Domain"
  colnames(table)[colnames(table)=="2"] <- "Phylum"
  colnames(table)[colnames(table)=="3"] <- "Class"
  colnames(table)[colnames(table)=="4"] <- "Order"
  colnames(table)[colnames(table)=="5"] <- "Family"
  colnames(table)[colnames(table)=="6"] <- "Genus"
  
  return(table)
  
  }

##### This function rounds p values based on value #####
roundP <- function(pvalue) {
  
  p <- ifelse(pvalue == 0, "p < 2e-16", 
              ifelse(pvalue < 0.001, "p < 0.001",
                     ifelse(pvalue >= 0.001, paste0("p = ", round(pvalue, digits = 3)),
                            "NS")))
  
  

}


sigStars <- function(pvalue) {
  
  pStar <- ifelse(pvalue < 0.001, "***",
                  ifelse(pvalue < 0.01, "**",
                         ifelse(pvalue < 0.05, "*",
                                "NS")))
  

}


getRelAbun <- function(table) {
  
  relabun <- table
  
  for (x in 1:nrow(table)){
    relabun[x,1:ncol(table)] = ((table[x,1:ncol(table)]) / (rowSums(table[,1:ncol(table)])[x]))
  }
  
  return(relabun)
  
}


getQuintileGroup <- function( val, quintile )
{
  for( i in 1:4)
  {
    if( val >= as.numeric(quintile[i]) & val <= as.numeric(quintile[i+1]) ) 
      return (i) 
  }
  
  return (5)
  
}

getQuartileGroup <- function( val, quartile )
{
  for( i in 1:3)
  {
    if( val >= as.numeric(quartile[i]) & val <= as.numeric(quartile[i+1]) ) 
      return (i) 
  }
  
  return (4)
  
}

getHalfGroup <- function( val, half )
{
  for( i in 1:1)
  {
    if( val >= as.numeric(half[i]) & val <= as.numeric(half[i+1]) ) 
      return (i) 
  }
  
  return (2)
  
}

getThirdGroup <- function( val, third )
{
  for( i in 1:2)
  {
    if( val >= as.numeric(third[i]) & val <= as.numeric(third[i+1]) ) 
      return (i) 
  }
  
  return (3)
  
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
