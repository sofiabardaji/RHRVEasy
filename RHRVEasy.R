#install.packages("RHRV")
library(RHRV)

# Post hoc Dunn test
library(dunn.test)
library(FSA)
library(PMCMR)

file_validation<-function(path){
  # 1. Check if path really exists
  if (dir.exists(path) != TRUE){
    stop("\nThe path ", path, " does not exist")
  }else{
    cat("\nThe path ", path, " exists ")
  }
  
  # 2. The path contains files:
  if ((length(list.files(path))>0) != TRUE){
    stop("but there are no files in it")
  }else{
    cat("and there are files in it\n\n")
  }
}

preparing_analysis<-function(file, rrs, format){
  hrv.data = CreateHRVData()
  hrv.data = SetVerbose(hrv.data, FALSE)
  
  hrv.data = LoadBeat(fileType = format, HRVData = hrv.data,  Recordname = file, RecordPath = rrs)
  
  hrv.data=BuildNIHR(hrv.data)
  hrv.data=FilterNIHR(hrv.data)
  hrv.data$Beat = hrv.data$Beat[2: nrow(hrv.data$Beat),]
  hrv.data
}

#Calls an RHRV function with hrv.data after cleaning the parameters
easy_call <- function(hrv.data, mf, ...) {
  args.list = plotrix::clean.args(list(...), mf)
  args.list$HRVData = hrv.data
  do.call(mf, args.list)
}

# Creating time analysis data frames
time_analysis<-function(format, files, class, rrs2, ...){
  dataFrame = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, CreateTimeAnalysis, ...)
    results=hrv.data$TimeAnalysis[]
    name_file = list ("filename" = file)
    group = list ("group" = class)
    # group_name = list("group" = group)
    row_list = c (name_file, results, group)
    df=as.data.frame(row_list)
    dataFrame=rbind(dataFrame, df)
  }
  dataFrame
}

# Frequency analysis
freq_analysis<-function(format, files, class, rrs2, ...){
  dataFrame = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data = easy_call(hrv.data, CalculatePSD, doPlot = F, ...)
    name_file = list ("filename" = file)
    x1 = easy_call(hrv.data, CalculateEnergyInPSDBands, ...)
    names(x1) = c("ULF", "VLF", "LF", "HF")
    group = list ("group" = class)
    row_list = c (name_file, x1, group)
    df = data.frame()
    df = rbind(df, as.data.frame(row_list))
    dataFrame=rbind(dataFrame, df)
  }
  dataFrame
}

#  Wavelet analysis
wavelet_analysis<-function(format, files, class, rrs2, ...){
  dataFrameMWavelet = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=SetVerbose(hrv.data, FALSE)
    
    hrv.data = easy_call(hrv.data, CalculatePowerBand, ...)
    
    index = length (hrv.data$FreqAnalysis)
    resultsWavelet = hrv.data$FreqAnalysis[[index]]
    resultsWavelet$File = file
    resultsWavelet$HRV = NULL
    resultsWavelet$ULF = sum(hrv.data$FreqAnalysis[[index]]$ULF)
    resultsWavelet$VLF = sum(hrv.data$FreqAnalysis[[index]]$VLF)
    resultsWavelet$LF = sum(hrv.data$FreqAnalysis[[index]]$LF)
    resultsWavelet$HF = sum(hrv.data$FreqAnalysis[[index]]$HF)
    resultsWavelet$LFHF = NULL
    resultsWavelet$Time = NULL
    name_file = list ("filename" = file)
    x1 = as.list(resultsWavelet)
    group = list ("group" = class)
    row_list = c (name_file, x1, group)
    dataFrameMWavelet = rbind(dataFrameMWavelet, as.data.frame(row_list))
    
  }
  dataFrameMWavelet
}

# PostHoc Dunn
dunnfreq<-function(dfM){
  dfM$group = factor(dfM$group)
  list (ULF = posthoc.kruskal.dunn.test(ULF ~ group, data=dfM, p.adjust="bonf"),
        VLF = posthoc.kruskal.dunn.test(VLF ~ group, data=dfM, p.adjust="bonf"),
        LF = posthoc.kruskal.dunn.test(LF ~ group, data=dfM, p.adjust="bonf"),
        HF = posthoc.kruskal.dunn.test(HF ~ group, data=dfM, p.adjust="bonf") )
}

dunntime<-function(dfM){
  dfM$group = factor(dfM$group)
  list (SDNN= posthoc.kruskal.dunn.test(SDNN ~ group, data = dfM, p.adjust="bonf"),
        SDANN = posthoc.kruskal.dunn.test(SDANN ~ group, data = dfM, p.adjust="bonf"),
        SDNNIDX = posthoc.kruskal.dunn.test(SDNNIDX ~ group, data = dfM, p.adjust="bonf"),
        pNN50 = posthoc.kruskal.dunn.test(pNN50 ~ group, data = dfM, p.adjust="bonf"),
        SDSD = posthoc.kruskal.dunn.test(SDSD ~ group, data = dfM, p.adjust="bonf"),
        rMSSD = posthoc.kruskal.dunn.test(rMSSD ~ group, data = dfM, p.adjust="bonf"),
        IRRR = posthoc.kruskal.dunn.test(IRRR ~ group, data = dfM, p.adjust="bonf"),
        MADRR = posthoc.kruskal.dunn.test(MADRR ~ group, data = dfM, p.adjust="bonf"),
        TINN = posthoc.kruskal.dunn.test(TINN ~ group, data = dfM, p.adjust="bonf"),
        HRVi = posthoc.kruskal.dunn.test(HRVi ~ group, data = dfM, p.adjust="bonf"))
}

statistical_analysisFreq<-function(dfM, verbose, numberOfExperimentalGroups){
  anova = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  kruskal = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  dunn = NA
  list = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  vec = c("group" = NA, "p-value ULF" = NA, "p-value VLF" = NA, "p-value LF" = NA, 
          "p-value HF" = NA)
  
  
  for(objeto in names(listDF)){
    vec$group = objeto
    
    for (column in c('ULF', 'VLF', 'LF', 'HF')){
      destino = paste0('p-value ', column)
      vec[[destino]] = (shapiro.test(listDF[[objeto]][[column]])$p.value)
    }

    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
  }
  
  # for (column in c('ULF', 'VLF', 'LF', 'HF')){
  #   if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.column > signif_level)) {
  #     if (verbose == TRUE){
  #       cat(column, "Normal: Anova. P-values = ", dataFramePvalues$p.value.ULF, "\n")
  #     }
  #     list$anova[[column]] = aov(column ~ group, data = dfM)
  #   } else {
  #     if (verbose == TRUE){
  #       cat(column, " NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.column, "\n")
  #     }
  #     list$kruskal$column = kruskal.test(column ~ group, data = dfM)
  #   }
  # }
  #### NO FUNCIONA EL AOV ~ !!!!!!!
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.ULF > signif_level)) {
    if (verbose == TRUE){
      cat(ULF, "Normal: Anova. P-values = ", dataFramePvalues$p.value.ULF, "\n")
    }
    list$anova[[ULF]] = aov(ULF ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat(ULF, " NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.ULF, "\n")
    }
    list$kruskal$ULF = kruskal.test(ULF ~ group, data = dfM)
  }
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.VLF > signif_level)) {
    if (verbose == TRUE){
      cat("VLF Normal: Anova. P-values = ", dataFramePvalues$p.value.VLF, "\n")
    }
    aov(VLF ~ group, data = dfM)
    list$anova$VLF = aov(VLF ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("VLF NOT normal: Kruskal. P-values = ",   dataFramePvalues$p.value.VLF,  "\n")
    }
    list$kruskal$VLF = kruskal.test(VLF ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.LF > signif_level)) {
    if (verbose == TRUE){
      cat("LF Normal: Anova. P-values = ",  dataFramePvalues$p.value.LF, "\n")
    }
    list$anova$LF = aov(LF ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("LF NOT normal: Kruskal. P-values = ",  dataFramePvalues$p.value.LF, "\n")
    }
    list$kruskal$LF = kruskal.test(LF ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.HF > signif_level)) {
    if (verbose == TRUE){
      cat("HF Normal: Anova. P-values = ",  dataFramePvalues$p.value.HF,  "\n")
    }
    list$anova$HF = aov(HF ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("HF NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HF,  "\n")
    }
    list$kruskal$HF = kruskal.test(HF ~ group, data = dfM)
  }
  list$dunn = dunnfreq(dfM)
  list
  
}

statistical_analysisTime<-function(dfM, verbose, numberOfExperimentalGroups){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  list = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  # Can't be a vector, must be a list to use the $ notation
  vec = list("group" = NA, "p-value SDNN" = NA, "p-value SDANN" = NA,
             "p-value SDNNIDX" = NA, "p-value pNN50" = NA, "p-value SDSD" = NA, 
             "p-value rMSSD" = NA, "p-value IRRR" = NA, "p-value MADRR" = NA, 
             "p-value TINN" = NA, "p-value HRVi" = NA)
  
  for(objeto in names(listDF)){
    vec$group = objeto
    
    for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                     'MADRR', 'TINN', 'HRVi')){
      destino = paste0('p-value ', column)
      vec[[destino]] = (shapiro.test(listDF[[objeto]][[column]])$p.value)
    }
    
    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
  }

  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDNN > signif_level)) {
    if (verbose == TRUE){
      cat("SDNN Normal: Anova. P-values = ", dataFramePvalues$p.value.SDNN, "\n")
    }
    list$anova$SDNN = aov(SDNN ~ group, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDNN, "\n")
    }
    list$kruskal$SDNN = kruskal.test(SDNN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDANN > signif_level)) {
    if (verbose == TRUE){
      cat("SDANN Normal: Anova. P-values = ", dataFramePvalues$p.value.SDANN, "\n")
    }
    list$anova$SDANN = aov(SDANN ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDANN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDANN, "\n")
    }
    list$kruskal$SDANN = kruskal.test(SDANN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDNNIDX > signif_level)) {
    if (verbose == TRUE){
      cat("SDNNIDX Normal: Anova. P-values = ", dataFramePvalues$p.value.SDNNIDX, "\n")
    }
    list$anova$SDNNIDX = aov(SDNNIDX ~ group, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNNIDX NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDNNIDX, "\n")
    }
    list$kruskal$SDNNIDX = kruskal.test(SDNNIDX ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.pNN50 > signif_level)) {
    if (verbose == TRUE){
      cat("pNN50 Normal: Anova. P-values = ", dataFramePvalues$p.value.pNN50, "\n")
    }
    list$anova$pNN50 = aov(pNN50 ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("pNN50 NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.pNN50, "\n")
    }
    list$kruskal$pNN50 = kruskal.test(pNN50 ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDSD > signif_level)) {
    if (verbose == TRUE){
      cat("SDSD Normal: Anova. P-values = ", dataFramePvalues$p.value.SDSD, "\n")
    }
    list$anova$SDSD = aov(SDSD ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDSD NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDSD, "\n")
    }
    list$kruskal$SDSD = kruskal.test(SDSD ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.rMSSD > signif_level)) {
    if (verbose == TRUE){
      cat("rMSSD Normal: Anova. P-values = ", dataFramePvalues$p.value.rMSSD, "\n")
    }
    list$anova$rMSSD = aov(rMSSD ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("rMSSD NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.rMSSD, "\n")
    }
    list$kruskal$rMSSD = kruskal.test(rMSSD ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.IRRR> signif_level)){
    if (verbose == TRUE){
      cat("IRRR Normal: Anova. P-values = ", dataFramePvalues$p.value.IRRR, "\n")
    }
    list$anova$IRRR = aov(IRRR ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("IRRR NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.IRRR, "\n")
    }
    list$kruskal$IRRR = kruskal.test(IRRR ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.MADRR > signif_level)){
    if (verbose == TRUE){
      cat("MADRR Normal: Anova. P-values = ", dataFramePvalues$p.value.MADRR, "\n")
    }
    list$anova$MADRR = aov(MADRR ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("MADRR NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.MADRR, "\n")
    }
    list$kruskal$MADRR = kruskal.test(MADRR ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.TINN > signif_level)){
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.TINN, "\n")
    }
    list$anova$TINN = aov(TINN ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.TINN, "\n")
    }
    list$kruskal$TINN = kruskal.test(TINN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.HRVi > signif_level)){
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HRVi, "\n")
    }
    list$anova$HRVi = aov(HRVi ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HRVi, "\n")
    }
    list$kruskal$HRVi = kruskal.test(HRVi ~ group, data = dfM)
  }
  
  list$dunn = dunntime(dfM)
  list
}

correctpValues <- function(listTime, listFreq, correction, method){
  
  listpValues = list(ULF = NA, VLF = NA, LF = NA, HF = NA,
                     SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, 
                     rMSSD = NA, IRRR = NA, MADRR = NA, TINN = NA, HRVi = NA)
  
  listpValuesCorrected = list(ULF = NA, VLF = NA, LF = NA, HF = NA, SDNN = NA, SDANN = NA, 
                              SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                              MADRR = NA, TINN = NA, HRVi = NA)
  
  for (column in c('ULF', 'VLF', 'LF', 'HF')){
    if(is.na(listFreq$anova[[column]])){
      listpValues[[column]] = listFreq$kruskal[[column]]$p.value
    }else{
      listpValues[[column]] = extract_ANOVA_pvalue(listFreq$anova[[column]])
    }
  }
  
  for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    if(is.na(listTime$anova[[column]])){
      listpValues[[column]] = listTime$kruskal[[column]]$p.value
    }else{
      listpValues[[column]] = extract_ANOVA_pvalue(listTime$anova[[column]])
    }
  }
  
  
  if (correction == TRUE){
    listpValuesCorrected = p.adjust(listpValues, method)
    listpValuesCorrected <- as.list(listpValuesCorrected)
    
  }else if (correction == FALSE){
    listpValuesCorrected = listpValues
  }
  listpValuesCorrected
}

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

extract_ANOVA_pvalue<-function(anovaObject){
  pvalue = summary(anovaObject)[[1]][1, 5]
  pvalue
}

print.RHRVEasyResult <- function(results){
  
  listDF = split(results$TimeAnalysis, results$TimeAnalysis$group)
  
  differencesFound = FALSE
  
  cat("\n\nResult of the analysis of the variability of the heart rate of the group",
      levels(results$TimeAnalysis$group)[1], 
      "versus the group", levels(results$TimeAnalysis$group)[2], ":\n\n")
  
  for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    if(is.na(results$StatysticalAnalysisTime$anova[[column]])){
      #report kruskal
      if(results$pValues[[column]]<signif_level){#error pvalue 1
        differencesFound = TRUE
        cat("There is a statistically significant difference in", column,  "; pvalue: ", 
            results$StatysticalAnalysisTime$kruskal[[column]]$p.value, "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          column = ""
          cat(column, " for the group", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
      }
    }
    #report anova
    else{
      if(results$pValues[[column]]<signif_level){
        differencesFound = TRUE
        cat("There is a statistically significant difference in", column, "; pvalue: ", 
            extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova[[column]]), "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group ", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
        
      }
    }
  }
  
  # 2 Conditions to report Dunn:
  # 1. We have more than 2 groups, we check that by looking at the length of listDF
  if(length(listDF)>2){
    for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                     'MADRR', 'TINN', 'HRVi')){
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
      var = which(a4$StatysticalAnalysisTime$dunn[[column]][["p.value"]]<signif_level, arr.ind = TRUE)
      if(length(var)>0){
      for (i in 1:nrow(var)){
        cat(a4$StatysticalAnalysisTime$dunn[[column]][["p.value"]][var][i], 'is the p value of',
            row.names(a4$StatysticalAnalysisTime$dunn[[column]][["p.value"]])[var[i,1]], 'vs',
            colnames(a4$StatysticalAnalysisTime$dunn[[column]][["p.value"]])[var[i,2]], '\n')
      }
      cat('for', column, '\n')
      }
    }
  }
  
  
  
  
  
  listDF = split(results$FrequencyAnalysis, results$FrequencyAnalysis$group)
  
  for (column in c('ULF', 'VLF', 'LF', 'HF')){
    if(is.na(results$StatysticalAnalysisFrequency$anova[[column]])){
      #report kruskal
      if(results$pValues[[column]]<signif_level){#error pvalue 1
        differencesFound = TRUE
        cat("There is a statistically significant difference in", column,  "; pvalue: ", 
            results$StatysticalAnalysisFrequency$kruskal[[column]]$p.value, "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          column = ""
          cat(column, " for the group", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
      }
    }
    #report anova
    else{
      if(results$pValues[[column]]<signif_level){
        differencesFound = TRUE
        cat("There is a statistically significant difference in", column, "; pvalue: ", 
            extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova[[column]]), "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group ", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
        
      }
    }
  }
  
  if(length(listDF)>2){
    for (column in c('ULF', 'VLF', 'LF', 'HF')){
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
      
      var = which(a4$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]]<signif_level, arr.ind = TRUE)
      if(length(var)>0){
        for (i in 1:nrow(var)){
          cat(a4$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]][var][i], 'is the p value of',
              row.names(a4$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]])[var[i,1]], 'vs',
              colnames(a4$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]])[var[i,2]], '\n')
        }
        cat('for', column, '\n')
      }
    }
    
  }
  
  
  if(!differencesFound){
    cat("No statistically significant difference were found\n")
  }
}


RHRVEasy<-function(folders, correction = FALSE, method = "bonferroni", verbose=FALSE, 
                   format = "RR", typeAnalysis = 'fourier', significance_level = 0.05, ...) {
  dataFrameMWavelet = data.frame()
  dataFrameMTime = data.frame()
  dataFrameMFreq = data.frame()
  #We create a global variable signif_level with the level significance
  signif_level <<- significance_level
  
  files = list()
  
  for (folder in folders){
    file_validation(folder)
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format, 
                                                                          list.files(folder), split_path(folder)[1], folder, ...))
  }
  
  numberOfExperimentalGroups = length(folders)
  # Statistical analysis of both
  
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime, 
                                                         verbose, numberOfExperimentalGroups)
  
  # FREQUENCY:
  if(typeAnalysis == "fourier"){
    for (folder in folders){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format, 
                                                                            list.files(folder), split_path(folder)[1], folder, ...))
    }
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, 
                                                             verbose, numberOfExperimentalGroups)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, 
                                                             verbose, numberOfExperimentalGroups)
    }
  }
  
  # WAVELET
  if(typeAnalysis == "wavelet"){
    for (folder in folders){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format,
                                                                                        list.files(folder), split_path(folder)[1], folder, 
                                                                                        type = typeAnalysis, ...))
    }
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, 
                                                             verbose, numberOfExperimentalGroups)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, 
                                                             verbose, numberOfExperimentalGroups)
    }
    dataFrameMFreq = dataFrameMWavelet
  }
  
  listpValues = correctpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis, 
                               correction, method)
  
  results = list("TimeAnalysis" = dataFrameMTime, "StatysticalAnalysisTime" = listTimeStatysticalAnalysis,
                 "FrequencyAnalysis" = dataFrameMFreq, 
                 "StatysticalAnalysisFrequency" = listFreqStatysticalAnalysis, "pValues" = listpValues)
  
  class(results) = "RHRVEasyResult"
  results
}
