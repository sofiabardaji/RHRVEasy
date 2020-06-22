#install.packages("RHRV")
library(RHRV)

# CREATING TIME ANALYSIS DATA FRAMES
preparing_analysis<-function(directory_files,file, rrs){
  hrv.data = CreateHRVData()
  hrv.data = SetVerbose(hrv.data, FALSE)
  hrv.data = LoadBeatRR(hrv.data,
                        RecordName = file,
                        RecordPath = rrs,
                        scale = 1
  )
  hrv.data=BuildNIHR(hrv.data)
  hrv.data=FilterNIHR(hrv.data)
  hrv.data$Beat = hrv.data$Beat[2: nrow(hrv.data$Beat),]
  hrv.data
}
time_analysis<-function(files, class, rrs2, dataFrame3){
  for (file in files) {
    hrv.data = preparing_analysis(directory_files = files, file = file, rrs = rrs2)
    hrv.data=CreateTimeAnalysis(hrv.data)
    results=hrv.data$TimeAnalysis[]
    name_file = list ("filename" = file)
    clase = list ("clase" = class)
    row_list = c (name_file, results, clase)
    dataFrame=as.data.frame(row_list)
    dataFrame3=rbind(dataFrame3, dataFrame)
  }
  dataFrame3
}


freq_analysis<-function(files, class, rrs2, dataFrame2){
  for (file in files) {
    hrv.data = preparing_analysis(directory_files = files, file = file, rrs = rrs2)
    # PlotNIHR(hrv.data)
    hrv.data=InterpolateNIHR(hrv.data)
    # Find the zeros in HR
    zero_indexes = which(hrv.data$HR == 0)
    # Compute the median of HR after removing 0s
    hr_median = median(hrv.data$HR[-zero_indexes])
    # Substitute the 0s in HR by the median
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=CalculatePSD(hrv.data,1,"pgram",doPlot = F)
    #print(CalculateEnergyInPSDBands(hrv.data, 1)) # 
    name_file = list ("filename" = file)
    x1 = as.list(CalculateEnergyInPSDBands(hrv.data, 1))
    names(x1) = c("ULF", "VLF", "LF", "HF")
    clase = list ("clase" = class)
    row_list = c (name_file, x1, clase)
    #unique(row_list)
    df = data.frame()
    df = rbind(df, as.data.frame(row_list))
    dataFrame2=rbind(dataFrame2, df)
  }
  dataFrame2
}
#  WAVELET ANALYSIS
wavelet_analysis<-function(files, class, rrs2, dataFrameMWavelet){
  for (file in files) {
    hrv.data = preparing_analysis(directory_files = files, file = file, rrs = rrs2)
    # PlotNIHR(hrv.data)
    hrv.data=InterpolateNIHR(hrv.data)
    # Find the zeros in HR
    zero_indexes = which(hrv.data$HR == 0)
    # Compute the median of HR after removing 0s
    hr_median = median(hrv.data$HR[-zero_indexes])
    # Substitute the 0s in HR by the median
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data=CreateFreqAnalysis(hrv.data)
    hrv.data=SetVerbose(hrv.data, FALSE)
    hrv.data = CalculatePowerBand(hrv.data,indexFreqAnalysis = length(hrv.data$FreqAnalysis),
                                  type = "wavelet", wavelet = "la8", bandtolerance = 0.01 )
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

# ANALYZING NORMALITY
shapiroFreqULF<-function(directorio){
  shapiroFreqULF = shapiro.test(directorio$ULF)
  shapiroFreqULF
}
shapiroFreqVLF<-function(directorio){
  shapiroFreqVLF = shapiro.test(directorio$VLF)
  shapiroFreqVLF
}
shapiroFreqLF<-function(directorio){
  shapiroFreqLF = shapiro.test(directorio$LF)
  shapiroFreqLF
}
shapiroFreqHF<-function(directorio){
  shapiroFreqHF = shapiro.test(directorio$HF)
  shapiroFreqHF
}
shapiroTimeSDNN<-function(directorio){
  shapiroTimeSDNN = shapiro.test(directorio$SDNN)
  shapiroTimeSDNN
}

shapiroTimeSDANN<-function(directorio){
  shapiroTimeSDANN = shapiro.test(directorio$SDANN)
  shapiroTimeSDANN
}
shapiroTimeSDNNIDX<-function(directorio){
  shapiroTimeSDNNIDX = shapiro.test(directorio$SDNNIDX)
  shapiroTimeSDNNIDX
}
shapiroTimepNN50<-function(directorio){
  shapiroTimepNN50 = shapiro.test(directorio$pNN50)
  shapiroTimepNN50
}
shapiroTimeSDSD<-function(directorio){
  shapiroTimeSDSD = shapiro.test(directorio$SDSD)
  shapiroTimeSDSD
}
shapiroTimerMSSD<-function(directorio){
  shapiroTimerMSSD = shapiro.test(directorio$rMSSD)
  shapiroTimerMSSD
}
shapiroTimeIRRR<-function(directorio){
  shapiroTimeIRRR = shapiro.test(directorio$IRRR)
  shapiroTimeIRRR
}
shapiroTimeMADRR<-function(directorio){
  shapiroTimeMADRR = shapiro.test(directorio$MADRR)
  shapiroTimeMADRR
}
shapiroTimeTINN<-function(directorio){
  shapiroTimeTINN = shapiro.test(directorio$TINN)
  shapiroTimeTINN
}
shapiroTimeHRVi<-function(directorio){
  shapiroTimeHRVi = shapiro.test(directorio$HRVi)
  shapiroTimeHRVi
}


# KRUSKAL WALLIS TEST
kruskalFreqULF<-function(dfM){
  kruskal.test(ULF ~ clase, data = dfM) 
}
kruskalFreqVLF<-function(dfM){
  kruskal.test(VLF ~ clase, data = dfM) 
}
kruskalFreqLF<-function(dfM){
  kruskal.test(LF ~ clase, data = dfM) 
}
kruskalFreqHF<-function(dfM){
  kruskal.test(HF ~ clase, data = dfM) 
}



kruskalTimeSDNN<-function(dfM){
  kruskal.test(SDNN ~ clase, data = dfM)
}
kruskalTimeSDANN<-function(dfM){
  kruskal.test(SDANN ~ clase, data = dfM)
} 
kruskalTimeSDNNIDX<-function(dfM){
  kruskal.test(SDNNIDX ~ clase, data = dfM)
}  
kruskalTimepNN50<-function(dfM){
  kruskal.test(pNN50 ~ clase, data = dfM)
}  
kruskalTimeSDSD<-function(dfM){
  kruskal.test(SDSD ~ clase, data = dfM)
}  

kruskalTimerMSSD<-function(dfM){
  kruskal.test(rMSSD ~ clase, data = dfM)
}  
kruskalTimeIRRR<-function(dfM){
  kruskal.test(IRRR ~ clase, data = dfM)
} 
kruskalTimeMADRR<-function(dfM){
  kruskal.test(MADRR ~ clase, data = dfM)
}  
kruskalTimeTINN<-function(dfM){
  kruskal.test(TINN ~ clase, data = dfM)
}  
kruskalTimeHRVi<-function(dfM){
  kruskal.test(HRVi ~ clase, data = dfM)
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

statistical_analysisFreqVerboseTRUE<-function(dfM, correctSigLevel){
  anova = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  kruskal = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$clase)
  shapiroFreqULFCase = shapiroFreqULF(listaDF[[1]])
  shapiroFreqULFControl = shapiroFreqULF(listaDF[[2]])
  pvaluesULF = c(shapiroFreqULFCase$p.value,shapiroFreqULFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesULF = p.adjust(pvaluesULF)
  }
  if (all(pvaluesULF > 0.05)) {
    cat("ULF Normal: Anova. P-values = ", pvaluesULF, "\n")
    lista$anova$ULF = aov(ULF ~ clase, data = dfM)
  }else {
    cat("ULF NOT normal: Kruskal. P-values = ", pvaluesULF, "\n")
    lista$kruskal$ULF = kruskalFreqULF(dfM)
  }
  shapiroFreqVLFCase = shapiroFreqVLF(listaDF[[1]])
  shapiroFreqVLFControl = shapiroFreqVLF(listaDF[[2]])
  pvaluesVLF = c(shapiroFreqVLFCase$p.value,shapiroFreqVLFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesVLF = p.adjust(pvaluesVLF)
  }
  if (all(pvaluesVLF > 0.05)) {
    cat("VLF Normal: Anova. P-values = ", pvaluesVLF, "\n")
    aov(VLF ~ clase, data = dfM)
    lista$anova$VLF = aov(VLF ~ clase, data = dfM)
  } else {
    cat("VLF NOT normal: Kruskal. P-values = ", pvaluesVLF, "\n")
    lista$kruskal$VLF = kruskalFreqVLF(dfM)
  }
  shapiroFreqLFCase = shapiroFreqLF(listaDF[[1]])
  shapiroFreqLFControl = shapiroFreqLF(listaDF[[2]])
  pvaluesLF = c(shapiroFreqLFCase$p.value,shapiroFreqLFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesLF = p.adjust(pvaluesLF)
  }
  if (all(pvaluesLF > 0.05)) {
    cat("LF Normal: Anova. P-values = ", pvaluesLF, "\n")
    lista$anova$LF = aov(LF ~ clase, data = dfM)  
  } else {
    cat("LF NOT normal: Kruskal. P-values = ", pvaluesLF, "\n")
    lista$kruskal$LF = kruskalFreqLF(dfM)
  }
  shapiroFreqHFCase = shapiroFreqHF(listaDF[[1]])
  shapiroFreqHFControl = shapiroFreqHF(listaDF[[2]])
  pvaluesHF = c(shapiroFreqHFCase$p.value,shapiroFreqHFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesHF = p.adjust(pvaluesHF)
  }
  if (all(pvaluesHF > 0.05)) {
    cat("HF Normal: Anova. P-values = ", pvaluesHF, "\n")
    lista$anova$HF = aov(HF ~ clase, data = dfM) 
  } else {
    cat("HF NOT normal: Kruskal. P-values = ", pvaluesHF, "\n")
    lista$kruskal$HF = kruskalFreqHF(dfM)
  }
  lista$dunn = dunnfreq(dfM)
  lista
  
}
statistical_analysisFreqVerboseFALSE<-function(dfM, correctSigLevel){
  anova = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  kruskal = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$clase)
  shapiroFreqULFCase = shapiroFreqULF(listaDF[[1]])
  shapiroFreqULFControl = shapiroFreqULF(listaDF[[2]])
  pvaluesULF = c(shapiroFreqULFCase$p.value,shapiroFreqULFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesULF = p.adjust(pvaluesULF)
  }
  if (all(pvaluesULF > 0.05)) {
    lista$anova$ULF = aov(ULF ~ clase, data = dfM)
  }else {
    lista$kruskal$ULF = kruskalFreqULF(dfM)
  }
  shapiroFreqVLFCase = shapiroFreqVLF(listaDF[[1]])
  shapiroFreqVLFControl = shapiroFreqVLF(listaDF[[2]])
  pvaluesVLF = c(shapiroFreqVLFCase$p.value,shapiroFreqVLFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesVLF = p.adjust(pvaluesVLF)
  }
  if (all(pvaluesVLF > 0.05)) {
    aov(VLF ~ clase, data = dfM)
    lista$anova$VLF = aov(VLF ~ clase, data = dfM)
  } else {
    lista$kruskal$VLF = kruskalFreqVLF(dfM)
  }
  shapiroFreqLFCase = shapiroFreqLF(listaDF[[1]])
  shapiroFreqLFControl = shapiroFreqLF(listaDF[[2]])
  pvaluesLF = c(shapiroFreqLFCase$p.value,shapiroFreqLFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesLF = p.adjust(pvaluesLF)
  }
  if (all(pvaluesLF > 0.05)) {
    lista$anova$LF = aov(LF ~ clase, data = dfM)  
  } else {
    lista$kruskal$LF = kruskalFreqLF(dfM)
  }
  shapiroFreqHFCase = shapiroFreqHF(listaDF[[1]])
  shapiroFreqHFControl = shapiroFreqHF(listaDF[[2]])
  pvaluesHF = c(shapiroFreqHFCase$p.value,shapiroFreqHFControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesHF = p.adjust(pvaluesHF)
  }
  if (all(pvaluesHF > 0.05)) {
    lista$anova$HF = aov(HF ~ clase, data = dfM) 
  } else {
    lista$kruskal$HF = kruskalFreqHF(dfM)
  }
  lista$dunn = dunnfreq(dfM)
  lista
  
}
statistical_analysisTimeVerboseTRUE<-function(dfM, correctSigLevel){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$clase)
  
  shapiroTimeSDNNCase = shapiroTimeSDNN(listaDF[[1]])
  shapiroTimeSDNNControl = shapiroTimeSDNN(listaDF[[2]])
  pvaluesSDNN = c(shapiroTimeSDNNCase$p.value,shapiroTimeSDNNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDNN = p.adjust(pvaluesSDNN)
  }
  if (all(pvaluesSDNN > 0.05)) { 
    cat("SDNN Normal: Anova. P-values = ", pvaluesSDNN, "\n")
    lista$anova$SDNN = aov(SDNN ~ clase, data = dfM)
  }else {
    cat("SDNN NOT normal: Kruskal. P-values = ", pvaluesSDNN, "\n")
    lista$kruskal$SDNN = kruskalTimeSDNN(dfM)
  }
  shapiroTimeSDANNCase = shapiroTimeSDANN(listaDF[[1]])
  shapiroTimeSDANNControl = shapiroTimeSDANN(listaDF[[2]])
  pvaluesSDANN = c(shapiroTimeSDANNCase$p.value,shapiroTimeSDANNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDANN = p.adjust(pvaluesSDANN)
  }
  if (all(pvaluesSDANN > 0.05)) { 
    cat("SDANN Normal: Anova. P-values = ", pvaluesSDANN, "\n")
    lista$anova$SDANN = aov(SDANN ~ clase, data = dfM)
  } else {
    cat("SDANN NOT normal: Kruskal. P-values = ", pvaluesSDANN, "\n")
    lista$kruskal$SDANN = kruskalTimeSDANN(dfM)
  }
  
  shapiroTimeSDNNIDXCase = shapiroTimeSDNNIDX(listaDF[[1]])
  shapiroTimeSDNNIDXControl = shapiroTimeSDNNIDX(listaDF[[2]])
  pvaluesSDNNIDX = c(shapiroTimeSDNNIDXCase$p.value,shapiroTimeSDNNIDXControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDNNIDX = p.adjust(pvaluesSDNNIDX)
  }
  if (all(pvaluesSDNNIDX > 0.05)) { 
    cat("SDNNIDX Normal: Anova. P-values = ", pvaluesSDNNIDX, "\n")
    lista$anova$SDNNIDX = aov(SDNNIDX ~ clase, data = dfM)
  }else {
    cat("SDNNIDX NOT normal: Kruskal. P-values = ", pvaluesSDNNIDX, "\n")
    lista$kruskal$SDNNIDX = kruskalTimeSDNNIDX(dfM)
  }
  
  shapiroTimepNN50Case = shapiroTimepNN50(listaDF[[1]])
  shapiroTimepNN50Control = shapiroTimepNN50(listaDF[[2]])
  pvaluespNN50 = c(shapiroTimepNN50Case$p.values,shapiroTimepNN50Control$p.value)
  if (correctSigLevel == TRUE){
    pvaluespNN50 = p.adjust(pvaluespNN50)
  }
  if (all(pvaluespNN50 > 0.05)) { 
    cat("pNN50 Normal: Anova. P-values = ", pvaluespNN50, "\n")
    lista$anova$pNN50 = aov(pNN50 ~ clase, data = dfM)
  } else {
    cat("pNN50 NOT normal: Kruskal. P-values = ", pvaluespNN50, "\n")
    lista$kruskal$pNN50 = kruskalTimepNN50(dfM)
  }
  
  shapiroTimeSDSDCase = shapiroTimeSDSD(listaDF[[1]])
  shapiroTimeSDSDControl = shapiroTimeSDSD(listaDF[[2]])
  pvaluesSDSD = c(shapiroTimeSDSDCase$p.value,shapiroTimeSDSDControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDSD = p.adjust(pvaluesSDSD)
  }
  if (all(pvaluesSDSD > 0.05)) {
    cat("SDSD Normal: Anova. P-values = ", pvaluesSDSD, "\n")
    lista$anova$SDSD = aov(SDSD ~ clase, data = dfM)
  } else {
    cat("SDSD NOT normal: Kruskal. P-values = ", pvaluesSDSD, "\n")
    lista$kruskal$SDSD = kruskalTimeSDSD(dfM)
  }
  
  
  shapiroTimerMSSDCase = shapiroTimerMSSD(listaDF[[1]])
  shapiroTimerMSSDControl = shapiroTimerMSSD(listaDF[[2]])
  pvaluesrMSSD = c(shapiroTimerMSSDCase$p.value,shapiroTimerMSSDControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesrMSSD = p.adjust(pvaluesrMSSD)
  }
  if (all(pvaluesrMSSD > 0.05)) { 
    cat("rMSSD Normal: Anova. P-values = ", pvaluesrMSSD, "\n")
    lista$anova$rMSSD = aov(rMSSD ~ clase, data = dfM)
  } else {
    cat("rMSSD NOT normal: Kruskal. P-values = ", pvaluesrMSSD, "\n")
    lista$kruskal$rMSSD = kruskalTimerMSSD(dfM)
  }
  
  shapiroTimeIRRRCase = shapiroTimeIRRR(listaDF[[1]])
  shapiroTimeIRRRControl = shapiroTimeIRRR(listaDF[[2]])
  pvaluesIRRR = c(shapiroTimeIRRRCase$p.value,shapiroTimeIRRRControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesIRRR = p.adjust(pvaluesIRRR)
  }
  if (all(pvaluesIRRR > 0.05)){
    cat("IRRR Normal: Anova. P-values = ", pvaluesIRRR, "\n")
    lista$anova$IRRR = aov(IRRR ~ clase, data = dfM)
  } else {
    cat("IRRR NOT normal: Kruskal. P-values = ", pvaluesIRRR, "\n")
    lista$kruskal$IRRR = kruskalTimeIRRR(dfM)
  }
  
  shapiroTimeMADRRCase = shapiroTimeMADRR(listaDF[[1]])
  shapiroTimeMADRRControl = shapiroTimeMADRR(listaDF[[2]])
  pvaluesMADRR = c(shapiroTimeMADRRCase$p.value,shapiroTimeMADRRControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesMADRR = p.adjust(pvaluesMADRR)
  }
  if (all(pvaluesMADRR > 0.05)){ 
    cat("MADRR Normal: Anova. P-values = ", pvaluesMADRR, "\n")
    lista$anova$MADRR = aov(MADRR ~ clase, data = dfM)
  } else {
    cat("IRRR NOT normal: Kruskal. P-values = ", pvaluesMADRR, "\n")
    lista$kruskal$MADRR = kruskalTimeMADRR(dfM)
  }
  
  shapiroTimeTINNCase = shapiroTimeTINN(listaDF[[1]])
  shapiroTimeTINNControl = shapiroTimeTINN(listaDF[[2]])
  pvaluesTINN = c(shapiroTimeTINNCase$p.value,shapiroTimeTINNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesTINN = p.adjust(pvaluesTINN)
  }
  if (all(pvaluesTINN > 0.05)){ 
    cat("TINN NOT normal: Kruskal. P-values = ", pvaluesTINN, "\n")
    lista$anova$TINN = aov(TINN ~ clase, data = dfM)
  } else {
    cat("TINN NOT normal: Kruskal. P-values = ", pvaluesTINN, "\n")
    lista$kruskal$TINN = kruskalTimeTINN(dfM)
  }
  
  shapiroTimeHRViCase = shapiroTimeHRVi(listaDF[[1]])
  shapiroTimeHRViControl = shapiroTimeHRVi(listaDF[[2]])
  pvaluesHRVi = c(shapiroTimeHRViCase$p.value,shapiroTimeHRViControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesHRVi = p.adjust(pvaluesHRVi)
  }
  if (all(pvaluesHRVi > 0.05)){ 
    cat("HRVi NOT normal: Kruskal. P-values = ", pvaluesHRVi, "\n")
    lista$anova$HRVi = aov(HRVi ~ clase, data = dfM)
  } else {
    cat("HRVi NOT normal: Kruskal. P-values = ", pvaluesHRVi, "\n")
    lista$kruskal$HRVi = kruskalTimeHRVi(dfM)
  }
  
  lista$dunn = dunntime(dfM)
  lista
  
}
statistical_analysisTimeVerboseFALSE<-function(dfM, correctSigLevel){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listaDF = split(dfM, dfM$clase)
  
  shapiroTimeSDNNCase = shapiroTimeSDNN(listaDF[[1]])
  shapiroTimeSDNNControl = shapiroTimeSDNN(listaDF[[2]])
  pvaluesSDNN = c(shapiroTimeSDNNCase$p.value,shapiroTimeSDNNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDNN = p.adjust(pvaluesSDNN)
  }
  if (all(pvaluesSDNN > 0.05)) { 
    lista$anova$SDNN = aov(SDNN ~ clase, data = dfM)
  }else {
    lista$kruskal$SDNN = kruskalTimeSDNN(dfM)
  }
  shapiroTimeSDANNCase = shapiroTimeSDANN(listaDF[[1]])
  shapiroTimeSDANNControl = shapiroTimeSDANN(listaDF[[2]])
  pvaluesSDANN = c(shapiroTimeSDANNCase$p.value,shapiroTimeSDANNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDANN = p.adjust(pvaluesSDANN)
  }
  if (all(pvaluesSDANN > 0.05)) { 
    lista$anova$SDANN = aov(SDANN ~ clase, data = dfM)
  } else {
    lista$kruskal$SDANN = kruskalTimeSDANN(dfM)
  }
  
  shapiroTimeSDNNIDXCase = shapiroTimeSDNNIDX(listaDF[[1]])
  shapiroTimeSDNNIDXControl = shapiroTimeSDNNIDX(listaDF[[2]])
  pvaluesSDNNIDX = c(shapiroTimeSDNNIDXCase$p.value,shapiroTimeSDNNIDXControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDNNIDX = p.adjust(pvaluesSDNNIDX)
  }
  if (all(pvaluesSDNNIDX > 0.05)) { 
    lista$anova$SDNNIDX = aov(SDNNIDX ~ clase, data = dfM)
  }else {
    lista$kruskal$SDNNIDX = kruskalTimeSDNNIDX(dfM)
  }
  
  shapiroTimepNN50Case = shapiroTimepNN50(listaDF[[1]])
  shapiroTimepNN50Control = shapiroTimepNN50(listaDF[[2]])
  pvaluespNN50 = c(shapiroTimepNN50Case$p.values,shapiroTimepNN50Control$p.value)
  if (correctSigLevel == TRUE){
    pvaluespNN50 = p.adjust(pvaluespNN50)
  }
  if (all(pvaluespNN50 > 0.05)) { 
    lista$anova$pNN50 = aov(pNN50 ~ clase, data = dfM)
  } else {
    lista$kruskal$pNN50 = kruskalTimepNN50(dfM)
  }
  
  shapiroTimeSDSDCase = shapiroTimeSDSD(listaDF[[1]])
  shapiroTimeSDSDControl = shapiroTimeSDSD(listaDF[[2]])
  pvaluesSDSD = c(shapiroTimeSDSDCase$p.value,shapiroTimeSDSDControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesSDSD = p.adjust(pvaluesSDSD)
  }
  if (all(pvaluesSDSD > 0.05)) {
    lista$anova$SDSD = aov(SDSD ~ clase, data = dfM)
  } else {
    lista$kruskal$SDSD = kruskalTimeSDSD(dfM)
  }
  
  
  shapiroTimerMSSDCase = shapiroTimerMSSD(listaDF[[1]])
  shapiroTimerMSSDControl = shapiroTimerMSSD(listaDF[[2]])
  pvaluesrMSSD = c(shapiroTimerMSSDCase$p.value,shapiroTimerMSSDControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesrMSSD = p.adjust(pvaluesrMSSD)
  }
  if (all(pvaluesrMSSD > 0.05)) { 
    lista$anova$rMSSD = aov(rMSSD ~ clase, data = dfM)
  } else {
    lista$kruskal$rMSSD = kruskalTimerMSSD(dfM)
  }
  
  shapiroTimeIRRRCase = shapiroTimeIRRR(listaDF[[1]])
  shapiroTimeIRRRControl = shapiroTimeIRRR(listaDF[[2]])
  pvaluesIRRR = c(shapiroTimeIRRRCase$p.value,shapiroTimeIRRRControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesIRRR = p.adjust(pvaluesIRRR)
  }
  if (all(pvaluesIRRR > 0.05)){
    lista$anova$IRRR = aov(IRRR ~ clase, data = dfM)
  } else {
    lista$kruskal$IRRR = kruskalTimeIRRR(dfM)
  }
  
  shapiroTimeMADRRCase = shapiroTimeMADRR(listaDF[[1]])
  shapiroTimeMADRRControl = shapiroTimeMADRR(listaDF[[2]])
  pvaluesMADRR = c(shapiroTimeMADRRCase$p.value,shapiroTimeMADRRControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesMADRR = p.adjust(pvaluesMADRR)
  }
  if (all(pvaluesMADRR > 0.05)){ 
    lista$anova$MADRR = aov(MADRR ~ clase, data = dfM)
  } else {
    lista$kruskal$MADRR = kruskalTimeMADRR(dfM)
  }
  
  shapiroTimeTINNCase = shapiroTimeTINN(listaDF[[1]])
  shapiroTimeTINNControl = shapiroTimeTINN(listaDF[[2]])
  pvaluesTINN = c(shapiroTimeTINNCase$p.value,shapiroTimeTINNControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesTINN = p.adjust(pvaluesTINN)
  }
  if (all(pvaluesTINN > 0.05)){ 
    lista$anova$TINN = aov(TINN ~ clase, data = dfM)
  } else {
    lista$kruskal$TINN = kruskalTimeTINN(dfM)
  }
  
  shapiroTimeHRViCase = shapiroTimeHRVi(listaDF[[1]])
  shapiroTimeHRViControl = shapiroTimeHRVi(listaDF[[2]])
  pvaluesHRVi = c(shapiroTimeHRViCase$p.value,shapiroTimeHRViControl$p.value)
  if (correctSigLevel == TRUE){
    pvaluesHRVi = p.adjust(pvaluesHRVi)
  }
  if (all(pvaluesHRVi > 0.05)){ 
    lista$anova$HRVi = aov(HRVi ~ clase, data = dfM)
  } else {
    lista$kruskal$HRVi = kruskalTimeHRVi(dfM)
  }
  
  lista$dunn = dunntime(dfM)
  lista
  
}

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}


print.RHRVEasyResult <- function(result){
  cat("Result of the analysis of the variability of the heart rate of the group",
      levels(result[[1]]$clase)[1], "versus the group", levels(result[[1]]$clase)[2], "\n\n")
  if(is.na(result[[2]]$anova$SDNN)){
    #report krustal
    if(result[[2]]$kruskal$SDNN$p.value<0.05){
      cat("There is a statistically significant difference in SDNN; pvalue: ", result[[2]]$kruskal$SDNN$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$SDNN$p.value<0.05){
      cat("There is a statistically significant difference in SDNN; pvalue: ", result[[2]]$anova$SDNN$p.value, "ºn")
    }
  }

  if(is.na(result[[2]]$anova$SDANN)){
    #report krustal
    if(result[[2]]$kruskal$SDANN$p.value<0.05){
      cat("There is a statistically significant difference in SDANN; pvalue: ", result[[2]]$kruskal$SDANN$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$SDANN$p.value<0.05){
      cat("There is a statistically significant difference in SDANN; pvalue: ", result[[2]]$anova$SDANN$p.value, "ºn")
    }
  }

  if(is.na(result[[2]]$anova$SDNNIDX)){
    #report krustal
    if(result[[2]]$kruskal$SDNNIDX$p.value<0.05){
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", result[[2]]$kruskal$SDNNIDX$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$SDNNIDX$p.value<0.05){
      cat("There is a statistically significant difference in SDNNIDX; pvalue: ", result[[2]]$anova$SDNNIDX$p.value, "ºn")
    }
  }

  if(is.na(result[[2]]$anova$pNN50)){
    #report krustal
    if(result[[2]]$kruskal$pNN50$p.value<0.05){
      cat("There is a statistically significant difference in pNN50; pvalue: ", result[[2]]$kruskal$pNN50$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$pNN50$p.value>0.05){
      cat("There is a statistically significant difference in pNN50; pvalue: ", result[[2]]$anova$pNN50$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$SDSD)){
    #report krustal
    if(result[[2]]$kruskal$SDSD$p.value<0.05){
      cat("There is a statistically significant difference in SDSD; pvalue: ", result[[2]]$kruskal$SDSD$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$SDSD$p.value<0.05){
      cat("There is a statistically significant difference in SDSD; pvalue: ", result[[2]]$anova$SDSD$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$rMSSD)){
    #report krustal
    if(result[[2]]$kruskal$rMSSD$p.value<0.05){
      cat("There is a statistically significant difference in rMSSD; pvalue: ", result[[2]]$kruskal$rMSSD$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$rMSSD$p.value<0.05){
      cat("There is a statistically significant difference in rMSSD; pvalue: ", result[[2]]$anova$rMSSD$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$IRRR)){
    #report krustal
    if(result[[2]]$kruskal$IRRR$p.value<0.05){
      cat("There is a statistically significant difference in IRRR; pvalue: ", result[[2]]$kruskal$IRRR$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$IRRR$p.value<0.05){
      cat("There is a statistically significant difference in IRRR; pvalue: ", result[[2]]$anova$IRRR$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$MADRR)){
    #report krustal
    if(result[[2]]$kruskal$MADRR$p.value<0.05){
      cat("There is a statistically significant difference in MADRR; pvalue: ", result[[2]]$kruskal$MADRR$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$MADRR$p.value<0.05){
      cat("There is a statistically significant difference in MADRR; pvalue: ", result[[2]]$anova$MADRR$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$TINN)){
    #report krustal
    if(result[[2]]$kruskal$TINN$p.value<0.05){
      cat("There is a statistically significant difference in TINN; pvalue: ", result[[2]]$kruskal$TINN$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$TINN$p.value<0.05){
      cat("There is a statistically significant difference in TINN; pvalue: ", result[[2]]$anova$TINN$p.value, "ºn")
    }
  }
  if(is.na(result[[2]]$anova$HRVi)){
    #report krustal
    if(result[[2]]$kruskal$HRVi$p.value<0.05){
      cat("There is a statistically significant difference in HRVi; pvalue: ", result[[2]]$kruskal$HRVi$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[2]]$anova$HRVi$p.value<0.05){
      cat("There is a statistically significant difference in HRVi; pvalue: ", result[[2]]$anova$HRVi$p.value, "ºn")
    }
  }
  if(is.na(result[[4]]$anova$ULF)){
    #report krustal
    if(result[[4]]$kruskal$ULF$p.value<0.05){
      cat("There is a statistically significant difference in ULF; pvalue: ", result[[4]]$kruskal$ULF$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[4]]$anova$ULF$p.value<0.05){
      cat("There is a statistically significant difference in ULF; pvalue: ", result[[4]]$anova$ULF$p.value, "ºn")
    }
  }
  if(is.na(result[[4]]$anova$VLF)){
    #report krustal
    if(result[[4]]$kruskal$VLF$p.value<0.05){
      cat("There is a statistically significant difference in VLF; pvalue: ", result[[4]]$kruskal$VLF$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[4]]$anova$VLF$p.value<0.05){
      cat("There is a statistically significant difference in VLF; pvalue: ", result[[4]]$anova$VLF$p.value, "ºn")
    }
  }
  if(is.na(result[[4]]$anova$LF)){
    #report krustal
    if(result[[4]]$kruskal$LF$p.value<0.05){
      cat("There is a statistically significant difference in LF; pvalue: ", result[[4]]$kruskal$LF$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[4]]$anova$LF$p.value<0.05){
      cat("There is a statistically significant difference in LF; pvalue: ", result[[4]]$anova$LF$p.value, "ºn")
    }
  }
  if(is.na(result[[4]]$anova$HF)){
    #report krustal
    if(result[[4]]$kruskal$HF$p.value<0.05){
    cat("There is a statistically significant difference in HF; pvalue: ", result[[4]]$kruskal$HF$p.value, "\n")
    }
  }
  #report anova
  else{
    if(result[[4]]$anova$HF$p.value<0.05){
      cat("There is a statistically significant difference in HF; pvalue: ", result[[4]]$anova$HF$p.value, "ºn")
    }
  }
}


RHRVEasy<-function(control, case, useWavelet = FALSE, correctSigLevel = TRUE, verbose=FALSE){
  dataFrame3 = data.frame()
  dataFrame2 = data.frame()
  dataFrameMWavelet = data.frame()
  
  filesControl = list.files(control)
  classControl = split_path(control)[1]
  
  filesCase = list.files(case)
  classCase = split_path(case)[1]
  
  listFreqStatysticalAnalysis = list()
  listTimeStatysticalAnalysis = list()

  dataFrameMTimeControl = time_analysis(filesControl, classControl, control, dataFrame3)
  dataFrameMTimeCase = time_analysis(filesCase, classCase, case, dataFrame3)
  # Creating a DF with both in Time
  dataFrameMTime=rbind(dataFrameMTimeControl, dataFrameMTimeCase)
  # Statistical analysis of both
  if(verbose == TRUE){
    listTimeStatysticalAnalysis = statistical_analysisTimeVerboseTRUE(dataFrameMTime, correctSigLevel)
  }else{
    listTimeStatysticalAnalysis = statistical_analysisTimeVerboseFALSE(dataFrameMTime, correctSigLevel)
  }
  # FREQUENCY:
  if(useWavelet == FALSE){
    dataFrameMFreqControl = freq_analysis(filesControl, classControl, control, dataFrame2)
    dataFrameMFreqCase = freq_analysis(filesCase, classCase, case, dataFrame2)
    dataFrameMFreq=rbind(dataFrameMFreqControl, dataFrameMFreqCase)
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreqVerboseTRUE(dataFrameMFreq, correctSigLevel)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreqVerboseFALSE(dataFrameMFreq, correctSigLevel)
    }
    
  }
  # WAVELET
  if(useWavelet == TRUE){
    dataFrameMWaveletControl = wavelet_analysis(filesControl, classControl, control, dataFrameMWavelet)
    dataFrameMWaveletCase = wavelet_analysis(filesCase, classCase, case, dataFrameMWavelet)
    dataFrameMWavelet=rbind(dataFrameMWaveletControl, dataFrameMWaveletCase)
    if(verbose == TRUE){
      listFreqStatysticalAnalysis = statistical_analysisFreqVerboseTRUE(dataFrameMFreq, correctSigLevel)
    }else{
      listFreqStatysticalAnalysis = statistical_analysisFreqVerboseFALSE(dataFrameMFreq, correctSigLevel)
    }    
    dataFrameMFreq = dataFrameMWavelet
  }
  results = list(dataFrameMTime, listTimeStatysticalAnalysis, dataFrameMFreq, listFreqStatysticalAnalysis)
  class(results) = "RHRVEasyResult"
  results
}
