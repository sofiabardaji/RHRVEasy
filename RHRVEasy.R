#install.packages("RHRV")
library(RHRV)


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
  
  hrv.data = LoadBeat(fileType = format, HRVData = hrv.data,  Recordname = file,RecordPath = rrs)
  
  hrv.data=BuildNIHR(hrv.data)
  hrv.data=FilterNIHR(hrv.data)
  hrv.data$Beat = hrv.data$Beat[2: nrow(hrv.data$Beat),]
  hrv.data
}
easy_call <- function(hrv.data, mf, ...) {
  args.list = plotrix::clean.args(list(...), mf)
  args.list$HRVData = hrv.data
  do.call(mf, args.list)
}

# CREATING TIME ANALYSIS DATA FRAMES
time_analysis<-function(format, files, class, rrs2, dataFrame3, ...){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, CreateTimeAnalysis, ...)
    results=hrv.data$TimeAnalysis[]
    name_file = list ("filename" = file)
    group = list ("group" = class)
    # group_name = list("group" = group)
    row_list = c (name_file, results, group)
    dataFrame=as.data.frame(row_list)
    dataFrame3=rbind(dataFrame3, dataFrame)
  }
  dataFrame3
}

# FREQUENCY ANALYSIS
freq_analysis<-function(format, files, class, rrs2, dataFrame2, ...){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data = easy_call(hrv.data, CalculatePSD, ...)
    name_file = list ("filename" = file)
    x1 = easy_call(hrv.data, CalculateEnergyInPSDBands, ...)
    names(x1) = c("ULF", "VLF", "LF", "HF")
    group = list ("group" = class)
    row_list = c (name_file, x1, group)
    df = data.frame()
    df = rbind(df, as.data.frame(row_list))
    dataFrame2=rbind(dataFrame2, df)
  }
  dataFrame2
}

#  WAVELET ANALYSIS
wavelet_analysis<-function(format, files, class, rrs2, dataFrameMWavelet, ...){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=SetVerbose(hrv.data, FALSE)
    
    
    
    hrv.data = CalculatePowerBand(hrv.data, indexFreqAnalysis = 1, size = 60, shift = 30, 
                                  sizesp = NULL, scale = "linear", ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
                                  LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4,
                                  type = "fourier", wavelet = "d4", bandtolerance = 0.01, relative = FALSE, verbose = FALSE)
    
    
    # hrv.data = easy_call(hrv.data, CalculatePowerBand, ...)
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

# POST HOC DUNN TEST
library(dunn.test)
library(FSA)
library(PMCMR)

# POSTHOC DUNN
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
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  vec = c("group" = NA, "p-value ULF" = NA, "p-value VLF" = NA, 
          "p-value LF" = NA, "p-value HF" = NA)
  
  
  
  for(objeto in sapply(dfM, levels)$group){
    
    vec$group = objeto
    vec$`p-value ULF` = shapiro.test(listaDF[[objeto]][["ULF"]])$p.value
    vec$`p-value VLF` = shapiro.test(listaDF[[objeto]][["VLF"]])$p.value
    vec$`p-value LF` = shapiro.test(listaDF[[objeto]][["LF"]])$p.value
    vec$`p-value HF` = shapiro.test(listaDF[[objeto]][["HF"]])$p.value

    
    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
    
  }

  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.ULF > 0.05)) {
    if (verbose == TRUE){
      cat("ULF Normal: Anova. P-values = ", dataFramePvalues$p.value.ULF, "\n")
    }
    lista$anova$ULF = aov(ULF ~ group, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("ULF NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.ULF, "\n")
    }
    lista$kruskal$ULF = kruskal.test(ULF ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.VLF > 0.05)) {
    if (verbose == TRUE){
      cat("VLF Normal: Anova. P-values = ", dataFramePvalues$p.value.VLF, "\n")
    }
    aov(VLF ~ group, data = dfM)
    lista$anova$VLF = aov(VLF ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("VLF NOT normal: Kruskal. P-values = ",   dataFramePvalues$p.value.VLF,  "\n")
    }
    lista$kruskal$VLF = kruskal.test(VLF ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.LF > 0.05)) {
    if (verbose == TRUE){
      cat("LF Normal: Anova. P-values = ",  dataFramePvalues$p.value.LF, "\n")
    }
    lista$anova$LF = aov(LF ~ group, data = dfM)  
  } else {
    if (verbose == TRUE){
      cat("LF NOT normal: Kruskal. P-values = ",  dataFramePvalues$p.value.LF, "\n")
    }
    lista$kruskal$LF = kruskal.test(LF ~ group, data = dfM)
  }
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.HF > 0.05)) {
    if (verbose == TRUE){
      cat("HF Normal: Anova. P-values = ",  dataFramePvalues$p.value.HF,  "\n")
    }
    lista$anova$HF = aov(HF ~ group, data = dfM) 
  } else {
    if (verbose == TRUE){
      cat("HF NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HF,  "\n")
    }
    lista$kruskal$HF = kruskal.test(HF ~ group, data = dfM)
  }
  lista$dunn = dunnfreq(dfM)
  lista
  
}
statistical_analysisTime<-function(dfM, verbose, numberOfExperimentalGroups){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  vec = c("group" = NA, "p-value SDNN" = NA, "p-value SDANN" = NA, 
          "p-value SDNNIDX" = NA, "p-value pNN50" = NA, "p-value SDSD" = NA, "p-value rMSSD" = NA, "p-value IRRR" = NA,
          "p-value MADRR" = NA, "p-value TINN" = NA, "p-value HRVi" = NA)
 
  for(objeto in sapply(dfM, levels)$group){

    vec$group = objeto
    vec$`p-value SDNN` = shapiro.test(listaDF[[objeto]][["SDNN"]])$p.value
    vec$`p-value SDANN` = shapiro.test(listaDF[[objeto]][["SDANN"]])$p.value
    vec$`p-value SDNNIDX` = shapiro.test(listaDF[[objeto]][["SDNNIDX"]])$p.value
    vec$`p-value pNN50` = shapiro.test(listaDF[[objeto]][["pNN50"]])$p.value
    vec$`p-value SDSD` = shapiro.test(listaDF[[objeto]][["SDSD"]])$p.value
    vec$`p-value rMSSD` = shapiro.test(listaDF[[objeto]][["rMSSD"]])$p.value
    vec$`p-value IRRR` = shapiro.test(listaDF[[objeto]][["IRRR"]])$p.value
    vec$`p-value MADRR` = shapiro.test(listaDF[[objeto]][["MADRR"]])$p.value
    vec$`p-value TINN` = shapiro.test(listaDF[[objeto]][["TINN"]])$p.value
    vec$`p-value HRVi` = shapiro.test(listaDF[[objeto]][["HRVi"]])$p.value

    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)

  }
  
  
    if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDNN > 0.05)) { 
    if (verbose == TRUE){
      cat("SDNN Normal: Anova. P-values = ", dataFramePvalues$p.value.SDNN, "\n")
    }
    lista$anova$SDNN = aov(SDNN ~ group, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDNN, "\n")
    }
    lista$kruskal$SDNN = kruskal.test(SDNN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDANN > 0.05)) { 
    if (verbose == TRUE){
      cat("SDANN Normal: Anova. P-values = ", dataFramePvalues$p.value.SDANN, "\n")
    }
    lista$anova$SDANN = aov(SDANN ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDANN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDANN, "\n")
    }
    lista$kruskal$SDANN = kruskal.test(SDANN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDNNIDX > 0.05)) { 
    if (verbose == TRUE){
      cat("SDNNIDX Normal: Anova. P-values = ", dataFramePvalues$p.value.SDNNIDX, "\n")
    }
    lista$anova$SDNNIDX = aov(SDNNIDX ~ group, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNNIDX NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDNNIDX, "\n")
    }
    lista$kruskal$SDNNIDX = kruskal.test(SDNNIDX ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.pNN50 > 0.05)) { 
    if (verbose == TRUE){
      cat("pNN50 Normal: Anova. P-values = ", dataFramePvalues$p.value.pNN50, "\n")
    }
    lista$anova$pNN50 = aov(pNN50 ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("pNN50 NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.pNN50, "\n")
    }
    lista$kruskal$pNN50 = kruskal.test(pNN50 ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.SDSD > 0.05)) {
    if (verbose == TRUE){
      cat("SDSD Normal: Anova. P-values = ", dataFramePvalues$p.value.SDSD, "\n")
    }
    lista$anova$SDSD = aov(SDSD ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDSD NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.SDSD, "\n")
    }
    lista$kruskal$SDSD = kruskal.test(SDSD ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.rMSSD > 0.05)) { 
    if (verbose == TRUE){
      cat("rMSSD Normal: Anova. P-values = ", dataFramePvalues$p.value.rMSSD, "\n")
    }
    lista$anova$rMSSD = aov(rMSSD ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("rMSSD NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.rMSSD, "\n")
    }
    lista$kruskal$rMSSD = kruskal.test(rMSSD ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.IRRR> 0.05)){
    if (verbose == TRUE){
      cat("IRRR Normal: Anova. P-values = ", dataFramePvalues$p.value.IRRR, "\n")
    }
    lista$anova$IRRR = aov(IRRR ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("IRRR NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.IRRR, "\n")
    }
    lista$kruskal$IRRR = kruskal.test(IRRR ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.MADRR > 0.05)){ 
    if (verbose == TRUE){
      cat("MADRR Normal: Anova. P-values = ", dataFramePvalues$p.value.MADRR, "\n")
    }
    lista$anova$MADRR = aov(MADRR ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("MADRR NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.MADRR, "\n")
    }
    lista$kruskal$MADRR = kruskal.test(MADRR ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.TINN > 0.05)){ 
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.TINN, "\n")
    }
    lista$anova$TINN = aov(TINN ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.TINN, "\n")
    }
    lista$kruskal$TINN = kruskal.test(TINN ~ group, data = dfM)
  }
  
  
  if (numberOfExperimentalGroups > 2 || all(dataFramePvalues$p.value.HRVi > 0.05)){ 
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HRVi, "\n")
    }
    lista$anova$HRVi = aov(HRVi ~ group, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", dataFramePvalues$p.value.HRVi, "\n")
    }
    lista$kruskal$HRVi = kruskal.test(HRVi ~ group, data = dfM)
  }
  
  lista$dunn = dunntime(dfM)
  lista
  
}


correctpValues <- function(listTime, listFreq, correction, method){
  
  listapValues = list(ULF = NA, VLF = NA, LF = NA, HF = NA, 
                      SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                      MADRR = NA, TINN = NA, HRVi = NA)
  
  listapValuesCorrected = list(ULF = NA, VLF = NA, LF = NA, HF = NA, 
                               SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                               MADRR = NA, TINN = NA, HRVi = NA)
  
  
  if(is.na(listFreq$anova$ULF)){
    listapValues$ULF = listFreq$kruskal$ULF$p.value
  }else{
    listapValues$ULF = extract_ANOVA_pvalue(listFreq$anova$ULF)
  }
  
  if(is.na(listFreq$anova$VLF)){
    listapValues$VLF = listFreq$kruskal$VLF$p.value
  }else{
    listapValues$VLF = extract_ANOVA_pvalue(listFreq$anova$VLF)
  }
  
  if(is.na(listFreq$anova$LF)){
    listapValues$LF = listFreq$kruskal$LF$p.value
  }else{
    listapValues$LF = extract_ANOVA_pvalue(listFreq$anova$LF)
  }
  
  if(is.na(listFreq$anova$HF)){
    listapValues$HF = listFreq$kruskal$HF$p.value
  }else{
    listapValues$HF = extract_ANOVA_pvalue(listFreq$anova$HF)
  }
  
  
  if(is.na(listTime$anova$SDNN)){
    listapValues$SDNN = listTime$kruskal$SDNN$p.value
  }else{
    listapValues$SDNN = extract_ANOVA_pvalue(listTime$anova$SDNN)
  }
  
  if(is.na(listTime$anova$SDANN)){
    listapValues$SDANN = listTime$kruskal$SDANN$p.value
  }else{
    listapValues$SDANN = extract_ANOVA_pvalue(listTime$anova$SDANN)
  }
  
  if(is.na(listTime$anova$SDNNIDX)){
    listapValues$SDNNIDX = listTime$kruskal$SDNNIDX$p.value
  }else{
    listapValues$SDNNIDX = extract_ANOVA_pvalue(listTime$anova$SDNNIDX)
  }
  
  if(is.na(listTime$anova$pNN50)){
    listapValues$pNN50 = listTime$kruskal$pNN50$p.value
  }else{
    listapValues$pNN50 = extract_ANOVA_pvalue(listTime$anova$pNN50)
  }
  
  if(is.na(listTime$anova$SDSD)){
    listapValues$SDSD = listTime$kruskal$SDSD$p.value
  }else{
    listapValues$SDSD = extract_ANOVA_pvalue(listTime$anova$SDSD)
  }
  
  if(is.na(listTime$anova$rMSSD)){
    listapValues$rMSSD = listTime$kruskal$rMSSD$p.value
  }else{
    listapValues$rMSSD = extract_ANOVA_pvalue(listTime$anova$rMSSD)
  }
  
  if(is.na(listTime$anova$IRRR)){
    listapValues$IRRR = listTime$kruskal$IRRR$p.value
  }else{
    listapValues$IRRR = extract_ANOVA_pvalue(listTime$anova$IRRR)
  }
  
  if(is.na(listTime$anova$MADRR)){
    listapValues$MADRR = listTime$kruskal$MADRR$p.value
  }else{
    listapValues$MADRR = extract_ANOVA_pvalue(listTime$anova$MADRR)
  }
  
  if(is.na(listTime$anova$TINN)){
    listapValues$TINN = listTime$kruskal$TINN$p.value
  }else{
    listapValues$TINN = extract_ANOVA_pvalue(listTime$anova$TINN)
  }
  
  if(is.na(listTime$anova$HRVi)){
    listapValues$HRVi = listTime$kruskal$HRVi$p.value
  }else{
    listapValues$HRVi = extract_ANOVA_pvalue(listTime$anova$HRVi)
  }
  
  if (correction == TRUE){
    listapValuesCorrected = p.adjust(listapValues, method)
    listapValuesCorrected <- as.list(listapValuesCorrected) 
    
  }else if (correction == FALSE){
    listapValuesCorrected = listapValues
  }
  listapValuesCorrected
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
  
  listaDF = split(results$TimeAnalysis, results$TimeAnalysis$group)
  
  differencesFound = FALSE
  
  cat("\n\nResult of the analysis of the variability of the heart rate of the group",
      levels(results$TimeAnalysis$group)[1], "versus the group", levels(results$TimeAnalysis$group)[2], ":\n\n")
  
  if(is.na(results$StatysticalAnalysisTime$anova$SDNN)){
    #report kruskal
    if(results$pValues$SDNN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNN; pvalue: ", results$StatysticalAnalysisTime$kruskal$SDNN$p.value, "\n")
      
      cat("SDNN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDNN), "+-", sd(listaDF$normal$SDNN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDNN), "+-", sd(listaDF$chf$SDNN), "\n\n")
    }
  }
  #report anova
  
  
  else{
    if(results$pValues$SDNN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNN; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$SDNN), "\n")
      
      cat("SDNN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDNN), "+-", sd(listaDF$normal$SDNN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDNN), "+-", sd(listaDF$chf$SDNN), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$SDANN)){
    #report kruskal
    if(results$pValues$SDANN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDANN; pvalue: ", results$StatysticalAnalysisTime$kruskal$SDANN$p.value, "\n")
      cat("SDANN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDANN), "+-", sd(listaDF$normal$SDANN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDANN), "+-", sd(listaDF$chf$SDANN), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$SDANN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDANN; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$SDANN), "ºn")
      cat("SDANN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDANN), "+-", sd(listaDF$normal$SDANN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDANN), "+-", sd(listaDF$chf$SDANN), "\n\n")
    }
  }
  
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$SDNNIDX)){
    #report kruskal
    if(results$pValues$SDNNIDX<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", results$StatysticalAnalysisTime$kruskal$SDNNIDX$p.value, "\n")
      cat("SDNNIDX for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDNNIDX), "+-", sd(listaDF$normal$SDNNIDX))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDNNIDX), "+-", sd(listaDF$chf$SDNNIDX), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$SDNNIDX<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$SDNNIDX), "ºn")
      cat("SDNNIDX for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDNNIDX), "+-", sd(listaDF$normal$SDNNIDX))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDNNIDX), "+-", sd(listaDF$chf$SDNNIDX), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$pNN50)){
    #report kruskal
    if(results$pValues$pNN50<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in pNN50; pvalue: ", results$StatysticalAnalysisTime$kruskal$pNN50$p.value, "\n")
      cat("pNN50 for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$pNN50), "+-", sd(listaDF$normal$pNN50))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$pNN50), "+-", sd(listaDF$chf$pNN50), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$pNN50<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in pNN50; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$pNN50), "ºn")
      cat("pNN50 for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$pNN50), "+-", sd(listaDF$normal$pNN50))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$pNN50), "+-", sd(listaDF$chf$pNN50), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$SDSD)){
    #report kruskal
    if(results$pValues$SDSD<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDSD; pvalue: ", results$StatysticalAnalysisTime$kruskal$SDSD$p.value, "\n")
      cat("SDSD for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDSD), "+-", sd(listaDF$normal$SDSD))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDSD), "+-", sd(listaDF$chf$SDSD), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$SDSD<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDSD; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$SDSD), "ºn")
      cat("SDSD for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$SDSD), "+-", sd(listaDF$normal$SDSD))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$SDSD), "+-", sd(listaDF$chf$SDSD), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$rMSSD)){
    #report kruskal
    if(results$pValues$rMSSD<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in rMSSD; pvalue: ", results$StatysticalAnalysisTime$kruskal$rMSSD$p.value, "\n")
      cat("rMSSD for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$rMSSD), "+-", sd(listaDF$normal$rMSSD))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$rMSSD), "+-", sd(listaDF$chf$rMSSD), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$rMSSD<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in rMSSD; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$rMSSD), "ºn")
      cat("rMSSD for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$rMSSD), "+-", sd(listaDF$normal$rMSSD))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$rMSSD), "+-", sd(listaDF$chf$rMSSD), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$IRRR)){
    #report kruskal
    if(results$pValues$IRRR<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in IRRR; pvalue: ", results$StatysticalAnalysisTime$kruskal$IRRR$p.value, "\n")
      cat("IRRR for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$IRRR), "+-", sd(listaDF$normal$IRRR))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$IRRR), "+-", sd(listaDF$chf$IRRR), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$IRRR<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in IRRR; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$SDNN), "ºn")
      cat("IRRR for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$IRRR), "+-", sd(listaDF$normal$IRRR))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$IRRR), "+-", sd(listaDF$chf$IRRR), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$MADRR)){
    #report kruskal
    if(results$pValues$MADRR<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in MADRR; pvalue: ", results$StatysticalAnalysisTime$kruskal$MADRR$p.value, "\n")
      cat("MADRR for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$MADRR), "+-", sd(listaDF$normal$MADRR))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$MADRR), "+-", sd(listaDF$chf$MADRR), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$MADRR<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in MADRR; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$MADRR), "ºn")
      cat("MADRR for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$MADRR), "+-", sd(listaDF$normal$MADRR))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$MADRR), "+-", sd(listaDF$chf$MADRR), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$TINN)){
    #report kruskal
    if(results$pValues$TINN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in TINN; pvalue: ", results$StatysticalAnalysisTime$kruskal$TINN$p.value, "\n")
      cat("TINN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$TINN), "+-", sd(listaDF$normal$TINN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$TINN), "+-", sd(listaDF$chf$TINN), "\n\n")
      
    }
  }
  #report anova
  else{
    if(results$pValues$TINN<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in TINN; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$TINN), "ºn")
      cat("TINN for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$TINN), "+-", sd(listaDF$normal$TINN))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$TINN), "+-", sd(listaDF$chf$TINN), "\n\n")
      
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisTime$anova$HRVi)){
    #report kruskal
    if(results$pValues$HRVi<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HRVi; pvalue: ", results$StatysticalAnalysisTime$kruskal$HRVi$p.value, "\n")
      cat("HRVi for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$HRVi), "+-", sd(listaDF$normal$HRVi))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$HRVi), "+-", sd(listaDF$chf$HRVi), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$HRVi<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HRVi; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova$HRVi), "ºn")
      cat("HRVi for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF$normal$HRVi), "+-", sd(listaDF$normal$HRVi))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF$chf$HRVi), "+-", sd(listaDF$chf$HRVi), "\n\n")
    }
  }
  
  
  listaDF1 = split(results$FrequencyAnalysis, results$FrequencyAnalysis$group)
  
  if(is.na(results$StatysticalAnalysisFrequency$anova$ULF)){
    #report kruskal
    if(results$pValues$ULF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in ULF; pvalue: ", results$StatysticalAnalysisFrequency$kruskal$ULF$p.value, "\n")
      cat("ULF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$ULF), "+-", sd(listaDF1$normal$ULF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$ULF), "+-", sd(listaDF1$chf$ULF), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$ULF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in ULF; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova$ULF), "ºn")
      cat("ULF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$ULF), "+-", sd(listaDF1$normal$ULF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$ULF), "+-", sd(listaDF1$chf$ULF), "\n\n")
    }
  }
  
  
  if(is.na(results$StatysticalAnalysisFrequency$anova$VLF)){
    #report kruskal
    if(results$pValues$VLF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in VLF; pvalue: ", results$StatysticalAnalysisFrequency$kruskal$VLF$p.value, "\n")
      cat("VLF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$VLF), "+-", sd(listaDF1$normal$VLF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$VLF), "+-", sd(listaDF1$chf$VLF), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$VLF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in VLF; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova$VLF), "ºn")
      cat("VLF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$VLF), "+-", sd(listaDF1$normal$VLF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$VLF), "+-", sd(listaDF1$chf$VLF), "\n\n")
    }
  }
  if(is.na(results$StatysticalAnalysisFrequency$anova$LF)){
    #report kruskal
    if(results$pValues$LF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in LF; pvalue: ", results$StatysticalAnalysisFrequency$kruskal$LF$p.value, "\n")
      cat("LF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$LF), "+-", sd(listaDF1$normal$LF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$LF), "+-", sd(listaDF1$chf$LF), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$LF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in LF; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova$LF), "ºn")
      cat("LF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$LF), "+-", sd(listaDF1$normal$LF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$LF), "+-", sd(listaDF1$chf$LF), "\n\n")
    }
  }
  if(is.na(results$StatysticalAnalysisFrequency$anova$HF)){
    #report kruskal
    if(results$pValues$HF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HF; pvalue: ", results$StatysticalAnalysisFrequency$kruskal$HF$p.value, "\n")
      cat("HF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$HF), "+-", sd(listaDF1$normal$HF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$HF), "+-", sd(listaDF1$chf$HF), "\n\n")
    }
  }
  #report anova
  else{
    if(results$pValues$HF<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HF; pvalue: ", extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova$HF), "ºn")
      cat("HF for the group ", levels(results$TimeAnalysis$group)[1], "is", 
          mean(listaDF1$normal$HF), "+-", sd(listaDF1$normal$HF))
      cat(" and for the group", levels(results$TimeAnalysis$group)[2], " is", 
          mean(listaDF1$chf$HF), "+-", sd(listaDF1$chf$HF), "\n\n")
    }
  }
  
  if(!differencesFound){
    cat("No statistically significant difference were found\n")
  }
}



RHRVEasy<-function(correction = FALSE, method = "bonferroni", verbose=FALSE, format = "RR", type  = "fourier", directorios, ...) {
  dataFrame = data.frame()
  dataFrame2 = data.frame()
  dataFrameMWavelet = data.frame()
  dataFrameMTime = data.frame()
  dataFrameMFreq = data.frame()
  
  files = list()
  
  for (directorio in directorios){
    file_validation(directorio)
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format, list.files(directorio), split_path(directorio)[1], directorio, dataFrame2, ...))
  }
  
  numberOfExperimentalGroups = length(directorios)
    
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime, verbose, numberOfExperimentalGroups)
  
  # FREQUENCY:
  if(type == "fourier"){

    for (directorio in directorios){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format, list.files(directorio), split_path(directorio)[1], directorio, dataFrame, ...))
    }

    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, verbose, numberOfExperimentalGroups)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, verbose, numberOfExperimentalGroups)
    }

  }

  # WAVELET
  if(type == "wavelet"){
    for (directorio in directorios){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format, list.files(directorio), split_path(directorio)[1], directorio, dataFrame, ...))
    }

    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, verbose, numberOfExperimentalGroups)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, verbose, numberOfExperimentalGroups)
    }
    dataFrameMFreq = dataFrameMWavelet
  }

 listapValues = correctpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis, correction, method)

  
  
  results =  list("TimeAnalysis" = dataFrameMTime, "StatysticalAnalysisTime" = listTimeStatysticalAnalysis,
                  "FrequencyAnalysis" = dataFrameMFreq, "StatysticalAnalysisFrequency" = listFreqStatysticalAnalysis,
                  "pValues" = listapValues)
    
  class(results) = "RHRVEasyResult"
  results
}
