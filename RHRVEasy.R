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
    cat("and there are files in it\n")
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

# CREATING TIME ANALYSIS DATA FRAMES
time_analysis<-function(format, files, class, rrs2, dataFrame3, size, numofbins, interval, verbose){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data=CreateTimeAnalysis(hrv.data, size=size, numofbins=numofbins, interval=interval, verbose=verbose)
    results=hrv.data$TimeAnalysis[]
    name_file = list ("filename" = file)
    clase = list ("clase" = class)
    row_list = c (name_file, results, clase)
    dataFrame=as.data.frame(row_list)
    dataFrame3=rbind(dataFrame3, dataFrame)
  }
  dataFrame3
}

# FREQUENCY ANALYSIS
freq_analysis<-function(format, files, class, rrs2, dataFrame2, freqhr, methodInterpolation, verbose, 
                        methodCalculationPSD, doPlot, ULFmin, ULFmax, VLFmin, VLFmax, 
                        LFmin,LFmax, HFmin, HFmax){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    # PlotNIHR(hrv.data)
    hrv.data=InterpolateNIHR(hrv.data, freqhr = freqhr, method = methodInterpolation, verbose = verbose)
    # Find the zeros in HR
    zero_indexes = which(hrv.data$HR == 0)
    # Compute the median of HR after removing 0s
    hr_median = median(hrv.data$HR[-zero_indexes])
    # Substitute the 0s in HR by the median
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=CalculatePSD(hrv.data, indexFreqAnalysis = 1, method = methodCalculationPSD, doPlot = doPlot)
    
    name_file = list ("filename" = file)
    x1 = as.list(CalculateEnergyInPSDBands(hrv.data, indexFreqAnalysis = 1, ULFmin = ULFmin, ULFmax = ULFmax, VLFmin = VLFmin, VLFmax = VLFmax, 
                                           LFmin = LFmin, LFmax = LFmax, HFmin = HFmin, HFmax = HFmax))
    
    names(x1) = c("ULF", "VLF", "LF", "HF")
    clase = list ("clase" = class)
    row_list = c (name_file, x1, clase)
    df = data.frame()
    df = rbind(df, as.data.frame(row_list))
    dataFrame2=rbind(dataFrame2, df)
  }
  dataFrame2
}

#  WAVELET ANALYSIS
wavelet_analysis<-function(format, files, class, rrs2, dataFrameMWavelet, freqhr, method, verbose, 
                           sizesp, scale, ULFmin, ULFmax, VLFmin, VLFmax, 
                           LFmin,LFmax, HFmin, HFmax, 
                           type, mother, bandtolerance, relative){
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    # PlotNIHR(hrv.data)
    hrv.data=InterpolateNIHR(hrv.data, freqhr = freqhr, method = method, verbose = verbose)
    # Find the zeros in HR
    zero_indexes = which(hrv.data$HR == 0)
    # Compute the median of HR after removing 0s
    hr_median = median(hrv.data$HR[-zero_indexes])
    # Substitute the 0s in HR by the median
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=SetVerbose(hrv.data, FALSE)
    hrv.data = CalculatePowerBand(hrv.data, indexFreqAnalysis = 1, 
                                  sizesp = sizesp, scale = scale, ULFmin = ULFmin, ULFmax = ULFmax, VLFmin = VLFmin, VLFmax = VLFmax, 
                                  LFmin = LFmin, LFmax = LFmax, HFmin = HFmin, HFmax = HFmax, 
                                  type = type, wavelet = mother, bandtolerance = bandtolerance, relative = relative, verbose = verbose)
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
    clase = list ("clase" = class)
    row_list = c (name_file, x1, clase)
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
  dfM$clase = factor(dfM$clase)  
  list (ULF = posthoc.kruskal.dunn.test(ULF ~ clase, data=dfM, p.adjust="bonf"),
        VLF = posthoc.kruskal.dunn.test(VLF ~ clase, data=dfM, p.adjust="bonf"),
        LF = posthoc.kruskal.dunn.test(LF ~ clase, data=dfM, p.adjust="bonf"),
        HF = posthoc.kruskal.dunn.test(HF ~ clase, data=dfM, p.adjust="bonf") )
}
dunntime<-function(dfM){
  dfM$clase = factor(dfM$clase)  
  list (SDNN= posthoc.kruskal.dunn.test(SDNN ~ clase, data = dfM, p.adjust="bonf"), 
        SDANN = posthoc.kruskal.dunn.test(SDANN ~ clase, data = dfM, p.adjust="bonf"), 
        SDNNIDX = posthoc.kruskal.dunn.test(SDNNIDX ~ clase, data = dfM, p.adjust="bonf"), 
        pNN50 = posthoc.kruskal.dunn.test(pNN50 ~ clase, data = dfM, p.adjust="bonf"),
        SDSD = posthoc.kruskal.dunn.test(SDSD ~ clase, data = dfM, p.adjust="bonf"), 
        rMSSD = posthoc.kruskal.dunn.test(rMSSD ~ clase, data = dfM, p.adjust="bonf"),
        IRRR = posthoc.kruskal.dunn.test(IRRR ~ clase, data = dfM, p.adjust="bonf"),
        MADRR = posthoc.kruskal.dunn.test(MADRR ~ clase, data = dfM, p.adjust="bonf"),
        TINN = posthoc.kruskal.dunn.test(TINN ~ clase, data = dfM, p.adjust="bonf"), 
        HRVi = posthoc.kruskal.dunn.test(HRVi ~ clase, data = dfM, p.adjust="bonf"))
}

statistical_analysisFreq<-function(dfM, correctSigLevel, verbose){
  anova = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  kruskal = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)

  listaDF = split(dfM, dfM$clase)
  
  shapiroFreqULFCase = shapiro.test(listaDF[[1]]$ULF)
  shapiroFreqULFControl = shapiro.test(listaDF[[2]]$ULF)
  pvaluesULF = c(shapiroFreqULFCase$p.value,shapiroFreqULFControl$p.value)
  
  shapiroFreqVLFCase = shapiro.test(listaDF[[1]]$VLF)
  shapiroFreqVLFControl = shapiro.test(listaDF[[2]]$VLF)
  pvaluesVLF = c(shapiroFreqVLFCase$p.value,shapiroFreqVLFControl$p.value)
  
  shapiroFreqLFCase = shapiro.test(listaDF[[1]]$LF)
  shapiroFreqLFControl = shapiro.test(listaDF[[2]]$LF)
  pvaluesLF = c(shapiroFreqLFCase$p.value,shapiroFreqLFControl$p.value)
  
  shapiroFreqHFCase = shapiro.test(listaDF[[1]]$HF)
  shapiroFreqHFControl = shapiro.test(listaDF[[2]]$HF)
  pvaluesHF = c(shapiroFreqHFCase$p.value,shapiroFreqHFControl$p.value)
  

  if (all(pvaluesULF > 0.05)) {
    if (verbose == TRUE){
      cat("ULF Normal: Anova. P-values = ", pvaluesULF, "\n")
    }
    lista$anova$ULF = aov(ULF ~ clase, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("ULF NOT normal: Kruskal. P-values = ", pvaluesULF, "\n")
    }
    lista$kruskal$ULF = kruskal.test(ULF ~ clase, data = dfM)
  }

  
  if (all(pvaluesVLF > 0.05)) {
    if (verbose == TRUE){
      cat("VLF Normal: Anova. P-values = ", pvaluesVLF, "\n")
    }
    aov(VLF ~ clase, data = dfM)
    lista$anova$VLF = aov(VLF ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("VLF NOT normal: Kruskal. P-values = ",   pvaluesVLF,  "\n")
    }
    lista$kruskal$VLF = kruskal.test(VLF ~ clase, data = dfM)
  }
 

  if (all(pvaluesLF > 0.05)) {
    if (verbose == TRUE){
      cat("LF Normal: Anova. P-values = ",  pvaluesLF, "\n")
    }
    lista$anova$LF = aov(LF ~ clase, data = dfM)  
  } else {
    if (verbose == TRUE){
      cat("LF NOT normal: Kruskal. P-values = ",  pvaluesLF, "\n")
    }
    lista$kruskal$LF = kruskal.test(LF ~ clase, data = dfM)
  }
  
  if (all(pvaluesHF > 0.05)) {
    if (verbose == TRUE){
      cat("HF Normal: Anova. P-values = ",  pvaluesHF,  "\n")
    }
    lista$anova$HF = aov(HF ~ clase, data = dfM) 
  } else {
    if (verbose == TRUE){
      cat("HF NOT normal: Kruskal. P-values = ", pvaluesHF,  "\n")
    }
    lista$kruskal$HF = kruskal.test(HF ~ clase, data = dfM)
  }
  lista$dunn = dunnfreq(dfM)
  lista
  
}

statistical_analysisTime<-function(dfM, correctSigLevel, verbose){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$clase)
  
  shapiroTimeSDNNCase = shapiro.test(listaDF[[1]]$SDNN)
  shapiroTimeSDNNControl = shapiro.test(listaDF[[2]]$SDNN)
  pvaluesSDNN = c(shapiroTimeSDNNCase$p.value,shapiroTimeSDNNControl$p.value)
  
  shapiroTimeSDANNCase = shapiro.test(listaDF[[1]]$SDANN)
  shapiroTimeSDANNControl = shapiro.test(listaDF[[2]]$SDANN)
  pvaluesSDANN = c(shapiroTimeSDANNCase$p.value,shapiroTimeSDANNControl$p.value)
  
  shapiroTimeSDNNIDXCase = shapiro.test(listaDF[[1]]$SDNNIDX)
  shapiroTimeSDNNIDXControl = shapiro.test(listaDF[[2]]$SDNNIDX)
  pvaluesSDNNIDX = c(shapiroTimeSDNNIDXCase$p.value,shapiroTimeSDNNIDXControl$p.value)
  
  shapiroTimepNN50Case = shapiro.test(listaDF[[1]]$pNN50)
  shapiroTimepNN50Control = shapiro.test(listaDF[[2]]$pNN50)
  pvaluespNN50 = c(shapiroTimepNN50Case$p.value, shapiroTimepNN50Control$p.value)
  
  shapiroTimeSDSDCase = shapiro.test(listaDF[[1]]$SDSD)
  shapiroTimeSDSDControl = shapiro.test(listaDF[[2]]$SDSD)
  pvaluesSDSD = c(shapiroTimeSDSDCase$p.value,shapiroTimeSDSDControl$p.value)
  
  shapiroTimerMSSDCase = shapiro.test(listaDF[[1]]$rMSSD)
  shapiroTimerMSSDControl = shapiro.test(listaDF[[2]]$rMSSD)
  pvaluesrMSSD = c(shapiroTimerMSSDCase$p.value,shapiroTimerMSSDControl$p.value)
  
  shapiroTimeIRRRCase = shapiro.test(listaDF[[1]]$IRRR)
  shapiroTimeIRRRControl = shapiro.test(listaDF[[2]]$IRRR)
  pvaluesIRRR = c(shapiroTimeIRRRCase$p.value,shapiroTimeIRRRControl$p.value)
  
  shapiroTimeMADRRCase = shapiro.test(listaDF[[1]]$MADRR)
  shapiroTimeMADRRControl = shapiro.test(listaDF[[2]]$MADRR)
  pvaluesMADRR = c(shapiroTimeMADRRCase$p.value,shapiroTimeMADRRControl$p.value)
  
  shapiroTimeTINNCase = shapiro.test(listaDF[[1]]$TINN)
  shapiroTimeTINNControl = shapiro.test(listaDF[[2]]$TINN)
  pvaluesTINN = c(shapiroTimeTINNCase$p.value,shapiroTimeTINNControl$p.value)
  
  shapiroTimeHRViCase = shapiro.test(listaDF[[1]]$HRVi)
  shapiroTimeHRViControl = shapiro.test(listaDF[[2]]$HRVi)
  pvaluesHRVi = c(shapiroTimeHRViCase$p.value,shapiroTimeHRViControl$p.value)
  
  
  
  if (all(pvaluesSDNN > 0.05)) { 
    if (verbose == TRUE){
      cat("SDNN Normal: Anova. P-values = ", pvaluesSDNN, "\n")
    }
    lista$anova$SDNN = aov(SDNN ~ clase, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNN NOT normal: Kruskal. P-values = ", pvaluesSDNN, "\n")
    }
    lista$kruskal$SDNN = kruskal.test(SDNN ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesSDANN > 0.05)) { 
    if (verbose == TRUE){
      cat("SDANN Normal: Anova. P-values = ", pvaluesSDANN, "\n")
    }
    lista$anova$SDANN = aov(SDANN ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDANN NOT normal: Kruskal. P-values = ", pvaluesSDANN, "\n")
    }
    lista$kruskal$SDANN = kruskal.test(SDANN ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesSDNNIDX > 0.05)) { 
    if (verbose == TRUE){
      cat("SDNNIDX Normal: Anova. P-values = ", pvaluesSDNNIDX, "\n")
    }
    lista$anova$SDNNIDX = aov(SDNNIDX ~ clase, data = dfM)
  }else {
    if (verbose == TRUE){
      cat("SDNNIDX NOT normal: Kruskal. P-values = ", pvaluesSDNNIDX, "\n")
    }
    lista$kruskal$SDNNIDX = kruskal.test(SDNNIDX ~ clase, data = dfM)
  }
  
  
  if (all(pvaluespNN50 > 0.05)) { 
    if (verbose == TRUE){
      cat("pNN50 Normal: Anova. P-values = ", pvaluespNN50, "\n")
    }
    lista$anova$pNN50 = aov(pNN50 ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("pNN50 NOT normal: Kruskal. P-values = ", pvaluespNN50, "\n")
    }
    lista$kruskal$pNN50 = kruskal.test(pNN50 ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesSDSD > 0.05)) {
    if (verbose == TRUE){
      cat("SDSD Normal: Anova. P-values = ", pvaluesSDSD, "\n")
    }
    lista$anova$SDSD = aov(SDSD ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("SDSD NOT normal: Kruskal. P-values = ", pvaluesSDSD, "\n")
    }
    lista$kruskal$SDSD = kruskal.test(SDSD ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesrMSSD > 0.05)) { 
    if (verbose == TRUE){
      cat("rMSSD Normal: Anova. P-values = ", pvaluesrMSSD, "\n")
    }
    lista$anova$rMSSD = aov(rMSSD ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("rMSSD NOT normal: Kruskal. P-values = ", pvaluesrMSSD, "\n")
    }
    lista$kruskal$rMSSD = kruskal.test(rMSSD ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesIRRR> 0.05)){
    if (verbose == TRUE){
      cat("IRRR Normal: Anova. P-values = ", pvaluesIRRR, "\n")
    }
    lista$anova$IRRR = aov(IRRR ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("IRRR NOT normal: Kruskal. P-values = ", pvaluesIRRR, "\n")
    }
    lista$kruskal$IRRR = kruskal.test(IRRR ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesMADRR > 0.05)){ 
    if (verbose == TRUE){
      cat("MADRR Normal: Anova. P-values = ", pvaluesMADRR, "\n")
    }
    lista$anova$MADRR = aov(MADRR ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("IRRR NOT normal: Kruskal. P-values = ", pvaluesMADRR, "\n")
    }
    lista$kruskal$MADRR = kruskal.test(MADRR ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesTINN > 0.05)){ 
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", pvaluesTINN, "\n")
    }
    lista$anova$TINN = aov(TINN ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("TINN NOT normal: Kruskal. P-values = ", pvaluesTINN, "\n")
    }
    lista$kruskal$TINN = kruskal.test(TINN ~ clase, data = dfM)
  }
  
  
  if (all(pvaluesHRVi > 0.05)){ 
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", pvaluesHRVi, "\n")
    }
    lista$anova$HRVi = aov(HRVi ~ clase, data = dfM)
  } else {
    if (verbose == TRUE){
      cat("HRVi NOT normal: Kruskal. P-values = ", pvaluesHRVi, "\n")
    }
    lista$kruskal$HRVi = kruskal.test(HRVi ~ clase, data = dfM)
  }
  
  lista$dunn = dunntime(dfM)
  lista
  
}


correctpValues <- function(listTime, listFreq){

  listapValues = list(ULF = NA, VLF = NA, LF = NA, HF = NA, 
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
  
  listapValuesCorrected = p.adjust(listapValues, "bonferroni")
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


print.RHRVEasyResult <- function(result){
  
  listaDF = split(result[[1]], result[[1]]$clase)
  
  differencesFound = FALSE
  
  cat("\n\nResult of the analysis of the variability of the heart rate of the group",
      levels(result[[1]]$clase)[1], "versus the group", levels(result[[1]]$clase)[2], ":\n\n")
 
   if(is.na(result[[2]]$anova$SDNN)){
    #report kruskal
    if(result[[2]]$kruskal$SDNN$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNN; pvalue: ", result[[2]]$kruskal$SDNN$p.value, "\n")
      
      cat("SDNN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDNN), "+-", sd(listaDF[[1]]$SDNN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDNN), "+-", sd(listaDF[[2]]$SDNN), "\n\n")
    }
  }
  #report anova
 
  
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$SDNN)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNN; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$SDNN), "\n")
      
      cat("SDNN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDNN), "+-", sd(listaDF[[1]]$SDNN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDNN), "+-", sd(listaDF[[2]]$SDNN), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$SDANN)){
    #report kruskal
    if(result[[2]]$kruskal$SDANN$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDANN; pvalue: ", result[[2]]$kruskal$SDANN$p.value, "\n")
      cat("SDANN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDANN), "+-", sd(listaDF[[1]]$SDANN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDANN), "+-", sd(listaDF[[2]]$SDANN), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$SDANN)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDANN; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$SDANN), "ºn")
      cat("SDANN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDANN), "+-", sd(listaDF[[1]]$SDANN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDANN), "+-", sd(listaDF[[2]]$SDANN), "\n\n")
    }
  }
  
  
  
  if(is.na(result[[2]]$anova$SDNNIDX)){
    #report kruskal
    if(result[[2]]$kruskal$SDNNIDX$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", result[[2]]$kruskal$SDNNIDX$p.value, "\n")
      cat("SDNNIDX for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDNNIDX), "+-", sd(listaDF[[1]]$SDNNIDX))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDNNIDX), "+-", sd(listaDF[[2]]$SDNNIDX), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$SDNNIDX)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$SDNNIDX), "ºn")
      cat("SDNNIDX for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDNNIDX), "+-", sd(listaDF[[1]]$SDNNIDX))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDNNIDX), "+-", sd(listaDF[[2]]$SDNNIDX), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$pNN50)){
    #report kruskal
    if(result[[2]]$kruskal$pNN50$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in pNN50; pvalue: ", result[[2]]$kruskal$pNN50$p.value, "\n")
      cat("pNN50 for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$pNN50), "+-", sd(listaDF[[1]]$pNN50))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$pNN50), "+-", sd(listaDF[[2]]$pNN50), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$pNN50)>0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in pNN50; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$pNN50), "ºn")
      cat("pNN50 for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$pNN50), "+-", sd(listaDF[[1]]$pNN50))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$pNN50), "+-", sd(listaDF[[2]]$pNN50), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$SDSD)){
    #report kruskal
    if(result[[2]]$kruskal$SDSD$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDSD; pvalue: ", result[[2]]$kruskal$SDSD$p.value, "\n")
      cat("SDSD for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDSD), "+-", sd(listaDF[[1]]$SDSD))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDSD), "+-", sd(listaDF[[2]]$SDSD), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$SDSD)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in SDSD; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$SDSD), "ºn")
      cat("SDSD for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$SDSD), "+-", sd(listaDF[[1]]$SDSD))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$SDSD), "+-", sd(listaDF[[2]]$SDSD), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$rMSSD)){
    #report kruskal
    if(result[[2]]$kruskal$rMSSD$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in rMSSD; pvalue: ", result[[2]]$kruskal$rMSSD$p.value, "\n")
      cat("rMSSD for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$rMSSD), "+-", sd(listaDF[[1]]$rMSSD))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$rMSSD), "+-", sd(listaDF[[2]]$rMSSD), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$rMSSD)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in rMSSD; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$rMSSD), "ºn")
      cat("rMSSD for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$rMSSD), "+-", sd(listaDF[[1]]$rMSSD))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$rMSSD), "+-", sd(listaDF[[2]]$rMSSD), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$IRRR)){
    #report kruskal
    if(result[[2]]$kruskal$IRRR$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in IRRR; pvalue: ", result[[2]]$kruskal$IRRR$p.value, "\n")
      cat("IRRR for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$IRRR), "+-", sd(listaDF[[1]]$IRRR))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$IRRR), "+-", sd(listaDF[[2]]$IRRR), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$IRRR)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in IRRR; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$SDNN), "ºn")
      cat("IRRR for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$IRRR), "+-", sd(listaDF[[1]]$IRRR))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$IRRR), "+-", sd(listaDF[[2]]$IRRR), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$MADRR)){
    #report kruskal
    if(result[[2]]$kruskal$MADRR$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in MADRR; pvalue: ", result[[2]]$kruskal$MADRR$p.value, "\n")
      cat("MADRR for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$MADRR), "+-", sd(listaDF[[1]]$MADRR))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$MADRR), "+-", sd(listaDF[[2]]$MADRR), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$MADRR)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in MADRR; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$MADRR), "ºn")
      cat("MADRR for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$MADRR), "+-", sd(listaDF[[1]]$MADRR))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$MADRR), "+-", sd(listaDF[[2]]$MADRR), "\n\n")
    }
  }
  
  
  if(is.na(result[[2]]$anova$TINN)){
    #report kruskal
    if(result[[2]]$kruskal$TINN$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in TINN; pvalue: ", result[[2]]$kruskal$TINN$p.value, "\n")
      cat("TINN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$TINN), "+-", sd(listaDF[[1]]$TINN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$TINN), "+-", sd(listaDF[[2]]$TINN), "\n\n")
      
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$TINN)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in TINN; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$TINN), "ºn")
      cat("TINN for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$TINN), "+-", sd(listaDF[[1]]$TINN))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$TINN), "+-", sd(listaDF[[2]]$TINN), "\n\n")
      
    }
  }
  
  
  if(is.na(result[[2]]$anova$HRVi)){
    #report kruskal
    if(result[[2]]$kruskal$HRVi$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HRVi; pvalue: ", result[[2]]$kruskal$HRVi$p.value, "\n")
      cat("HRVi for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$HRVi), "+-", sd(listaDF[[1]]$HRVi))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$HRVi), "+-", sd(listaDF[[2]]$HRVi), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[2]]$anova$HRVi)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HRVi; pvalue: ", extract_ANOVA_pvalue(result[[2]]$anova$HRVi), "ºn")
      cat("HRVi for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF[[1]]$HRVi), "+-", sd(listaDF[[1]]$HRVi))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF[[2]]$HRVi), "+-", sd(listaDF[[2]]$HRVi), "\n\n")
    }
  }
  
  
  listaDF1 = split(result[[3]], result[[3]]$clase)
  
  if(is.na(result[[4]]$anova$ULF)){
    #report kruskal
    if(result[[4]]$kruskal$ULF$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in ULF; pvalue: ", result[[4]]$kruskal$ULF$p.value, "\n")
      cat("ULF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$ULF), "+-", sd(listaDF1[[1]]$ULF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$ULF), "+-", sd(listaDF1[[2]]$ULF), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[4]]$anova$ULF)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in ULF; pvalue: ", extract_ANOVA_pvalue(result[[4]]$anova$ULF), "ºn")
      cat("ULF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$ULF), "+-", sd(listaDF1[[1]]$ULF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$ULF), "+-", sd(listaDF1[[2]]$ULF), "\n\n")
    }
  }
  
  
  if(is.na(result[[4]]$anova$VLF)){
    #report kruskal
    if(result[[4]]$kruskal$VLF$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in VLF; pvalue: ", result[[4]]$kruskal$VLF$p.value, "\n")
      cat("VLF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$VLF), "+-", sd(listaDF1[[1]]$VLF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$VLF), "+-", sd(listaDF1[[2]]$VLF), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[4]]$anova$VLF)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in VLF; pvalue: ", extract_ANOVA_pvalue(result[[4]]$anova$VLF), "ºn")
      cat("VLF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$VLF), "+-", sd(listaDF1[[1]]$VLF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$VLF), "+-", sd(listaDF1[[2]]$VLF), "\n\n")
    }
  }
  if(is.na(result[[4]]$anova$LF)){
    #report kruskal
    if(result[[4]]$kruskal$LF$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in LF; pvalue: ", result[[4]]$kruskal$LF$p.value, "\n")
      cat("LF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$LF), "+-", sd(listaDF1[[1]]$LF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$LF), "+-", sd(listaDF1[[2]]$LF), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[4]]$anova$LF)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in LF; pvalue: ", extract_ANOVA_pvalue(result[[4]]$anova$LF), "ºn")
      cat("LF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$LF), "+-", sd(listaDF1[[1]]$LF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$LF), "+-", sd(listaDF1[[2]]$LF), "\n\n")
    }
  }
  if(is.na(result[[4]]$anova$HF)){
    #report kruskal
    if(result[[4]]$kruskal$HF$p.value<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HF; pvalue: ", result[[4]]$kruskal$HF$p.value, "\n")
      cat("HF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$HF), "+-", sd(listaDF1[[1]]$HF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$HF), "+-", sd(listaDF1[[2]]$HF), "\n\n")
    }
  }
  #report anova
  else{
    if(extract_ANOVA_pvalue(result[[4]]$anova$HF)<0.05){
      differencesFound = TRUE
      cat("There is a statistically significant difference in HF; pvalue: ", extract_ANOVA_pvalue(result[[4]]$anova$HF), "ºn")
      cat("HF for the group ", levels(result[[1]]$clase)[1], "is", 
          mean(listaDF1[[1]]$HF), "+-", sd(listaDF1[[1]]$HF))
      cat(" and for the group", levels(result[[1]]$clase)[2], " is", 
          mean(listaDF1[[2]]$HF), "+-", sd(listaDF1[[2]]$HF), "\n\n")
    }
  }
  
  if(!differencesFound){
    cat("No statistically significant difference were found\n")
  }
}


RHRVEasy<-function(control, case, correctSigLevel = TRUE, verbose=FALSE, format = "RR",
                   size = 300, numofbins = NULL, interval = 7.8125, verboseTime = NULL,
                   freqhr = 4, methodInterpolation = c("linear", "spline"), verboseFreq = NULL,
                   methodCalculationPSD = c("pgram", "ar", "lomb"), doPlot = F,
                   ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
                   LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4,
                   sizesp = NULL, scale = "linear", 
                   type = "fourier", mother = "d4", bandtolerance = 0.01, relative = FALSE, verboseWavelet = NULL) {
  dataFrame3 = data.frame()
  dataFrame2 = data.frame()
  dataFrameMWavelet = data.frame()
  
  file_validation(control)
  file_validation(case)
  
  
  filesControl = list.files(control)
  classControl = split_path(control)[1]
  
  filesCase = list.files(case)
  classCase = split_path(case)[1]
  
  listFreqStatysticalAnalysis = list()
  listTimeStatysticalAnalysis = list()
  


  dataFrameMTimeControl = time_analysis(format, filesControl, classControl, control, dataFrame3, size, numofbins, interval, verboseTime)
  dataFrameMTimeCase = time_analysis(format, filesCase, classCase, case, dataFrame3, size, numofbins, interval, verboseTime)
  # Creating a DF with both in Time
  dataFrameMTime=rbind(dataFrameMTimeControl, dataFrameMTimeCase)

  # Statistical analysis of both

  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime, correctSigLevel, verbose)

  # FREQUENCY:
  if(type == "fourier"){
    
    dataFrameMFreqControl = freq_analysis(format, filesControl, classControl, control, dataFrame2, freqhr, methodInterpolation, verboseFreq,
                                          methodCalculationPSD, doPlot,ULFmin, ULFmax, VLFmin, VLFmax, LFmin, LFmax, HFmin, HFmax)

    dataFrameMFreqCase = freq_analysis(format, filesCase, classCase, case, dataFrame2, freqhr, methodInterpolation, verboseFreq,
                                       methodCalculationPSD, doPlot,ULFmin, ULFmax, VLFmin, VLFmax, LFmin, LFmax, HFmin, HFmax)
    
    dataFrameMFreq=rbind(dataFrameMFreqControl, dataFrameMFreqCase)
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, correctSigLevel, verbose)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq, correctSigLevel, verbose)
    }
    
  }
  
  # WAVELET
  if(type == "wavelet"){
    dataFrameMWaveletControl = wavelet_analysis(format, filesControl, classControl, control, dataFrameMWavelet, freqhr, methodInterpolation, verboseWavelet,
                                                sizesp, scale, 
                                                ULFmin, ULFmax, VLFmin, VLFmax, LFmin, LFmax, HFmin, HFmax, 
                                                type, mother, bandtolerance, relative)
    dataFrameMWaveletCase = wavelet_analysis(format, filesCase, classCase, case, dataFrameMWavelet, freqhr, methodInterpolation, verboseWavelet,
                                             sizesp, scale, 
                                             ULFmin, ULFmax, VLFmin, VLFmax, LFmin, LFmax, HFmin, HFmax, 
                                             type, mother, bandtolerance, relative)
    dataFrameMWavelet=rbind(dataFrameMWaveletControl, dataFrameMWaveletCase)
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, correctSigLevel, verbose)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet, correctSigLevel, verbose)
    }    
    dataFrameMFreq = dataFrameMWavelet
  }
  results = list(dataFrameMTime, listTimeStatysticalAnalysis, dataFrameMFreq, listFreqStatysticalAnalysis)
  class(results) = "RHRVEasyResult"
  results
}
