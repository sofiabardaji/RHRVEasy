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
    # donde he puesto files ponia allfilesnormal
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
# dataFrameMFreq = data.frame()
# dataFrameMFreq = freq_analysis(files=list.files("rrs/normal/"), class = "normal", rrs2 = "rrs/normal/", dataFrameMFreq)
# dataFrameMFreq2 = freq_analysis(files=list.files("rrs/chf/"), class = "chf", rrs2 = "rrs/chf/", dataFrameMFreq)
# dataFrameMFreq=rbind(dataFrameMFreq, dataFrameMFreq2)



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
    dataFrame2 <- dataFrame2[!duplicated(dataFrame2), ]
  }
  dataFrame2
}
# dataFrameMWavelet = data.frame()
# dataFrameMWavelet = wavelet_analysis(files=list.files("rrs/new/"), class = "chf", rrs2 = "rrs/new/", dataFrameMWavelet)
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
    # dataFrameMWavelet <- dataFrameMWavelet[!duplicated(dataFrameMWavelet), ]
    
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
  list (ULF = posthoc.kruskal.dunn.test(ULF ~ clase, data=dfM, p.adjust="bonf"),
        VLF = posthoc.kruskal.dunn.test(VLF ~ clase, data=dfM, p.adjust="bonf"),
        LF = posthoc.kruskal.dunn.test(LF ~ clase, data=dfM, p.adjust="bonf"),
        HF = posthoc.kruskal.dunn.test(HF ~ clase, data=dfM, p.adjust="bonf") )
}
dunntime<-function(dfM){
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



statistical_analysisFreq<-function(dfst, dfM, correctSigLevel){
  shapiroFreqULF = shapiroFreqULF(dfst)
  shapiroFreqVLF = shapiroFreqVLF(dfst)
  shapiroFreqLF = shapiroFreqLF(dfst)
  shapiroFreqHF = shapiroFreqHF(dfst)
  
  anova = list(ULF = 0, VLF = 0, LF = 0, HF = 0)
  kruskal = list(ULF = 0, VLF = 0, LF = 0, HF = 0)
  dunn = 0
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)

  if (shapiroFreqULF$p.value > 0.05) {
    print("ULF Normal: Anova")
    list$anova$ULF = aov(ULF ~ clase, data = dfM)
  }
  if (shapiroFreqULF$p.value < 0.05) {
    print("ULF NOT normal: Kruskal")
    lista$kruskal$ULF = kruskalFreqULF(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqULF$p.value)
  }
  if (shapiroFreqVLF$p.value > 0.05) {
    print("VLF Normal: Anova")
    aov(VLF ~ clase, data = dfM)
    list$anova$VLF = aov(VLF ~ clase, data = dfM)
  }
  if (shapiroFreqVLF$p.value < 0.05) {
    print("VLF NOT normal: Kruskal")
    lista$kruskal$VLF = kruskalFreqVLF(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqVLF$p.value)
  }
  if (shapiroFreqLF$p.value > 0.05) {
    print("LF Normal: Anova")
    list$anova$LF = aov(LF ~ clase, data = dfM)  }
  
  if (shapiroFreqLF$p.value < 0.05) {
    print("LF NOT normal: Kruskal")
    lista$kruskal$LF = kruskalFreqLF(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqLF$p.value)
  }
  if (shapiroFreqHF$p.value > 0.05) {
    print("HF Normal: Anova")
    list$anova$HF = aov(HF ~ clase, data = dfM)  }
  if (shapiroFreqHF$p.value < 0.05) {
    print("HF NOT normal: Kruskal")
    lista$kruskal$HF = kruskalFreqHF(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqHF$p.value)
  }
  lista$dunn = dunnfreq(dfM)
  lista
  
}
statistical_analysisTime<-function(dfst, dfM, correctSigLevel){
  shapiroTimeSDNN = shapiroTimeSDNN(dfst)
  shapiroTimeSDANN = shapiroTimeSDANN(dfst)
  shapiroTimeSDNNIDX = shapiroTimeSDNNIDX(dfst)
  shapiroTimepNN50 = shapiroTimepNN50(dfst)
  shapiroTimeSDSD = shapiroTimeSDSD(dfst)
  shapiroTimerMSSD = shapiroTimerMSSD(dfst)
  shapiroTimeIRRR = shapiroTimeIRRR(dfst)
  shapiroTimeMADRR = shapiroTimeMADRR(dfst)
  shapiroTimeTINN = shapiroTimeTINN(dfst)
  shapiroTimeHRVi = shapiroTimeHRVi(dfst)
  
  anova = list(SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, MSSD = 0, IRRR = 0,
               MADRR = 0, TINN = 0, HRVi = 0)
  kruskal = list(SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, MSSD = 0, IRRR = 0,
                 MADRR = 0, TINN = 0, HRVi = 0)
  dunn = 0
  lista = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  if (shapiroTimeSDNN$p.value > 0.05) { 
    print("SDNN Normal: Anova")
    lista$anova$SDNN = aov(SDNN ~ clase, data = dfM)
  }
  if (shapiroTimeSDNN$p.value < 0.05) {
    print("SDNN NOT normal: Kruskal")
    lista$kruskal$SDNN = kruskalTimeSDNN(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDNN$p.value)
  }
  if (shapiroTimeSDANN$p.value > 0.05) { 
    print("SDANN Normal: Anova")
    lista$anova$SDANN = aov(SDANN ~ clase, data = dfM)
  }
  if (shapiroTimeSDANN$p.value < 0.05) {
    print("SDANN NOT normal: Kruskal")
    lista$kruskal$SDANN = kruskalTimeSDANN(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDANN$p.value)
  }
  
  if (shapiroTimeSDNNIDX$p.value > 0.05) { 
    print("SDNNIDX Normal: Anova")
    lista$anova$SDNNIDX = aov(SDNNIDX ~ clase, data = dfM)
  }
  if (shapiroTimeSDNNIDX$p.value < 0.05) {
    print("SDNNIDX NOT normal: Kruskal")
    lista$kruskal$SDNNIDX = kruskalTimeSDNNIDX(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDNNIDX$p.value)
  }
  
  if (shapiroTimepNN50$p.value > 0.05) { 
    print("pNN50 Normal: Anova")
    lista$anova$pNN50 = aov(pNN50 ~ clase, data = dfM)
  }
  if (shapiroTimepNN50$p.value < 0.05) {
    print("pNN50 NOT normal: Kruskal")
    lista$kruskal$pNN50 = kruskalTimepNN50(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimepNN50$p.value)
  }
  
  
  if (shapiroTimeSDSD$p.value > 0.05) { 
    print("SDSD Normal: Anova")
    lista$anova$SDSD = aov(SDSD ~ clase, data = dfM)
  }
  if (shapiroTimeSDSD$p.value < 0.05) {
    print("SDSD NOT normal: Kruskal")
    lista$kruskal$SDSD = kruskalTimeSDSD(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDSD$p.value)
  }
  
  
  if (shapiroTimerMSSD$p.value > 0.05) { 
    print("MSSD Normal: Anova")
    lista$anova$MSSD = aov(MSSD ~ clase, data = dfM)
  }
  if (shapiroTimerMSSD$p.value < 0.05) {
    print("MSSD NOT normal: Kruskal")
    lista$kruskal$MSSD = kruskalTimerMSSD(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimerMSSD$p.value)
  }
  
  
  if (shapiroTimeIRRR$p.value > 0.05) { 
    print("IRRR Normal: Anova")
    lista$anova$IRRR = aov(IRRR ~ clase, data = dfM)
  }
  if (shapiroTimeIRRR$p.value < 0.05) {
    print("IRRR NOT normal: Kruskal")
    lista$kruskal$IRRR = kruskalTimeIRRR(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeIRRR$p.value)
  }
  
  
  if (shapiroTimeMADRR$p.value > 0.05) { 
    print("MADRR Normal: Anova")
    lista$anova$MADRR = aov(MADRR ~ clase, data = dfM)
  }
  if (shapiroTimeMADRR$p.value < 0.05) {
    print("MADRR NOT normal: Kruskal")
    lista$kruskal$MADRR = kruskalTimeMADRR(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeMADRR$p.value)
  }
  
  if (shapiroTimeTINN$p.value > 0.05) { 
    print("TINN Normal: Anova")
    lista$anova$TINN = aov(TINN ~ clase, data = dfM)
  }
  if (shapiroTimeTINN$p.value < 0.05) {
    print("TINN NOT normal: Kruskal")
    lista$kruskal$TINN = kruskalTimeTINN(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeTINN$p.value)
  }
  
  if (shapiroTimeHRVi$p.value > 0.05) { 
    print("HRVi Normal: Anova")
    lista$anova$HRVi = aov(HRVi ~ clase, data = dfM)
  }
  if (shapiroTimeHRVi$p.value < 0.05) {
    print("HRVi NOT normal: Kruskal")
    lista$kruskal$HRVi = kruskalTimeHRVi(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeHRVi$p.value)
  }
  lista$dunn = dunntime(dfM)
  lista
  
}


RHRVEasy<-function(control, case, useWavelet = FALSE, correctSigLevel = TRUE){
  dataFrame3 = data.frame()
  dataFrame2 = data.frame()
  dataFrameMWavelet = data.frame()
  files1 = list.files(control)
  class1 <- strsplit(control, .Platform$file.sep)[[1]][2]
  files2 = list.files(case)
  listaFreqControl = list()
  listaFreqCase = list()
  listaTimeControl = list()
  listaTimeCase = list()
  class2 <- strsplit(case, .Platform$file.sep)[[1]][2]
  dataFrameMTime1 = time_analysis(files1, class1, control, dataFrame3)
  dataFrameMTime2 = time_analysis(files2, class2, case, dataFrame3)
  # Creating a DF with both in Time
  dataFrameMTime=rbind(dataFrameMTime1, dataFrameMTime2)
  # Statistical analysis of both
  listaTimeControl = statistical_analysisTime(dataFrameMTime1, dataFrameMTime, correctSigLevel)
  listaTimeCase = statistical_analysisTime(dataFrameMTime2, dataFrameMTime, correctSigLevel)
  # FREQUENCY:
  if(useWavelet == FALSE){
    dataFrameMFreq1 = freq_analysis(files1, class1, control, dataFrame2)
    dataFrameMFreq2 = freq_analysis(files2, class2, case, dataFrame2)
    dataFrameMFreq=rbind(dataFrameMFreq1, dataFrameMFreq2)
    listaFreqControl = statistical_analysisFreq(dataFrameMFreq1, dataFrameMFreq, correctSigLevel)
    listaFreqCase = statistical_analysisFreq(dataFrameMFreq2, dataFrameMFreq, correctSigLevel)
  }
  # WAVELET
  if(useWavelet == TRUE){
    dataFrameMWavelet1 = wavelet_analysis(files1, class1, control, dataFrameMWavelet)
    print("Done 1")
    dataFrameMWavelet2 = wavelet_analysis(files2, class2, case, dataFrameMWavelet)
    print("Done 2")
    dataFrameMWavelet=rbind(dataFrameMWavelet1, dataFrameMWavelet2)
    print("Done 3")
    listaFreqControl = statistical_analysisFreq(dataFrameMWavelet1, dataFrameMWavelet, correctSigLevel)
    listaFreqCase = statistical_analysisFreq(dataFrameMWavelet2, dataFrameMWavelet, correctSigLevel)
    dataFrameMFreq = dataFrameMWavelet
  }
  list(dataFrameMTime, listaTimeControl, listaTimeCase, dataFrameMFreq, listaFreqControl, listaFreqCase)
}

