
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
time_analysis<-function(files, class, rrs2){
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
freq_analysis<-function(files, class, rrs2){
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
    df = rbind(df, as.data.frame(row_list))
    dataFrame2=rbind(dataFrame2, df)
    dataFrame2 <- dataFrame2[!duplicated(dataFrame2), ]
  }
  dataFrame2
}

#  WAVELET ANALYSIS
wavelet_analysis<-function(files, class, rrs2){
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
    # resultsWavelet$HRV = NULL
    resultsWavelet$ULF = sum(hrv.data$FreqAnalysis[[index]]$ULF)
    resultsWavelet$VLF = sum(hrv.data$FreqAnalysis[[index]]$VLF)
    resultsWavelet$LF = sum(hrv.data$FreqAnalysis[[index]]$LF)
    resultsWavelet$HF = sum(hrv.data$FreqAnalysis[[index]]$HF)
    # resultsWavelet$LFHF = NULL
    # resultsWavelet$Time = NULL
    
    name_file = list ("filename" = file)
    x1 = as.list(resultsWavelet)
    clase = list ("clase" = class)
    row_list = c (name_file, x1, clase)
    df = rbind(df, as.data.frame(row_list))
    dataFrameMWavelet=rbind(dataFrameMWavelet, df)
    dataFrameMWavelet <- dataFrameMWavelet[!duplicated(dataFrameMWavelet), ]
    
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
shapiroTimeMSSD<-function(directorio){
  shapiroTimeMSSD = shapiro.test(directorio$rMSSD)
  shapiroTimeMSSD
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
kruskalfreq<-function(df){
  list (kruskal.test(ULF ~ clase, data = df), kruskal.test(VLF ~ clase, data = df),
        kruskal.test(LF ~ clase, data = df), kruskal.test(HF ~ clase, data = df))
}
kruskaltime<-function(df){
  list (kruskal.test(SDNN ~ clase, data = df), kruskal.test(SDANN ~ clase, data = df), 
        kruskal.test(SDNNIDX ~ clase, data = df), kruskal.test(pNN50 ~ clase, data = df),
        kruskal.test(SDSD ~ clase, data = df), kruskal.test(rMSSD ~ clase, data = df),
        kruskal.test(IRRR ~ clase, data = df),kruskal.test(MADRR ~ clase, data = df),
        kruskal.test(TINN ~ clase, data = df), kruskal.test(HRVi ~ clase, data = df))
}


# POST HOC DUNN TEST
library(dunn.test)
library(FSA)
library(PMCMR)

# POSTHOC DUNN
dunnfreq<-function(df){
  list (posthoc.kruskal.dunn.test(ULF ~ clase, data=df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(VLF ~ clase, data=df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(LF ~ clase, data=df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(HF ~ clase, data=df, p.adjust="bonf") )
}
dunntime<-function(df){
  list (posthoc.kruskal.dunn.test(SDNN ~ clase, data = df, p.adjust="bonf"), 
        posthoc.kruskal.dunn.test(SDANN ~ clase, data = df, p.adjust="bonf"), 
        posthoc.kruskal.dunn.test(SDNNIDX ~ clase, data = df, p.adjust="bonf"), 
        posthoc.kruskal.dunn.test(pNN50 ~ clase, data = df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(SDSD ~ clase, data = df, p.adjust="bonf"), 
        posthoc.kruskal.dunn.test(rMSSD ~ clase, data = df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(IRRR ~ clase, data = df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(MADRR ~ clase, data = df, p.adjust="bonf"),
        posthoc.kruskal.dunn.test(TINN ~ clase, data = df, p.adjust="bonf"), 
        posthoc.kruskal.dunn.test(HRVi ~ clase, data = df, p.adjust="bonf"))
}

statistical_analysisFreq<-function(df, dfM, correctSigLevel){
  shapiroFreqULF = shapiroFreqULF(df)
  shapiroFreqVLF = shapiroFreqVLF(df)
  shapiroFreqLF = shapiroFreqLF(df)
  shapiroFreqHF = shapiroFreqHF(df)
  
  if (shapiroFreqULF$p.value > 0.05) {
    print("ULF Normal: Anova")
    aov(ULF ~ clase, data = dfM)
  }
  if (shapiroFreqULF$p.value < 0.05) {
    print("ULF NOT normal: Kruskal")
    kruskalfreq(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqULF$p.value)
  }
  if (shapiroFreqVLF$p.value > 0.05) {
    print("VLF Normal: Anova")
    aov(VLF ~ clase, data = dfM)
  }
  if (shapiroFreqVLF$p.value < 0.05) {
    print("VLF NOT normal: Kruskal")
    kruskalfreq(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqVLF$p.value)
  }
  if (shapiroFreqLF$p.value > 0.05) {
    print("LF Normal: Anova")
    aov(LF ~ clase, data = dfM)
  }
  
  if (shapiroFreqLF$p.value < 0.05) {
    print("LF NOT normal: Kruskal")
    kruskalfreq(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqLF$p.value)
  }
  if (shapiroFreqHF$p.value > 0.05) {
    print("HF Normal: Anova")
    aov(HF ~ clase, data = dfM)
  }
  if (shapiroFreqHF$p.value < 0.05) {
    print("HF NOT normal: Kruskal")
    kruskalfreq(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroFreqHF$p.value)
  }
  dunnfreq(dfM)
  
}
statistical_analysisTime<-function(df, dfM, correctSigLevel){
  shapiroTimeSDNN = shapiroTimeSDNN(df)
  shapiroTimeSDANN = shapiroTimeSDANN(df)
  shapiroTimeSDNNIDX = shapiroTimeSDNNIDX(df)
  shapiroTimepNN50 = shapiroTimepNN50(df)
  shapiroTimeSDSD = shapiroTimeSDSD(df)
  shapiroTimeMSSD = shapiroTimeMSSD(df)
  shapiroTimeIRRR = shapiroTimeIRRR(df)
  shapiroTimeMADRR = shapiroTimeMADRR(df)
  shapiroTimeTINN = shapiroTimeTINN(df)
  shapiroTimeHRVi = shapiroTimeHRVi(df)
  
  if (shapiroTimeSDNN$p.value > 0.05) { 
    print("SDNN Normal: Anova")
    aov(SDNN ~ clase, data = dfM)
  }
  if (shapiroTimeSDNN$p.value < 0.05) {
    print("SDNN NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDNN$p.value)
  }
  if (shapiroTimeSDANN$p.value > 0.05) { 
    print("SDANN Normal: Anova")
    aov(SDANN ~ clase, data = dfM)
  }
  if (shapiroTimeSDANN$p.value < 0.05) {
    print("SDANN NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDANN$p.value)
  }
  
  if (shapiroTimeSDNNIDX$p.value > 0.05) { 
    print("SDNNIDX Normal: Anova")
    aov(SDNNIDX ~ clase, data = dfM)
  }
  if (shapiroTimeSDNNIDX$p.value < 0.05) {
    print("SDNNIDX NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDNNIDX$p.value)
  }
  
  if (shapiroTimepNN50$p.value > 0.05) { 
    print("pNN50 Normal: Anova")
    aov(pNN50 ~ clase, data = dfM)
  }
  if (shapiroTimepNN50$p.value < 0.05) {
    print("pNN50 NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimepNN50$p.value)
  }
  
  
  if (shapiroTimeSDSD$p.value > 0.05) { 
    print("SDSD Normal: Anova")
    aov(SDSD ~ clase, data = dfM)
  }
  if (shapiroTimeSDSD$p.value < 0.05) {
    print("SDSD NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeSDSD$p.value)
  }
  
  
  if (shapiroTimeMSSD$p.value > 0.05) { 
    print("MSSD Normal: Anova")
    aov(MSSD ~ clase, data = dfM)
  }
  if (shapiroTimeMSSD$p.value < 0.05) {
    print("MSSD NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeMSSD$p.value)
  }
  
  
  if (shapiroTimeIRRR$p.value > 0.05) { 
    print("IRRR Normal: Anova")
    aov(IRRR ~ clase, data = dfM)
  }
  if (shapiroTimeIRRR$p.value < 0.05) {
    print("IRRR NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeIRRR$p.value)
  }
  
  
  if (shapiroTimeMADRR$p.value > 0.05) { 
    print("MADRR Normal: Anova")
    aov(MADRR ~ clase, data = dfM)
  }
  if (shapiroTimeMADRR$p.value < 0.05) {
    print("MADRR NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeMADRR$p.value)
  }
  
  if (shapiroTimeTINN$p.value > 0.05) { 
    print("TINN Normal: Anova")
    aov(TINN ~ clase, data = dfM)
  }
  if (shapiroTimeTINN$p.value < 0.05) {
    print("TINN NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeTINN$p.value)
  }
  
  if (shapiroTimeHRVi$p.value > 0.05) { 
    print("HRVi Normal: Anova")
    aov(HRVi ~ clase, data = dfM)
  }
  if (shapiroTimeHRVi$p.value < 0.05) {
    print("HRVi NOT normal: Kruskal")
    kruskaltime(dfM)
  }
  if (correctSigLevel == TRUE){
    p.adjust(shapiroTimeHRVi$p.value)
  }
  
  dunntime(dfM)
  
}


RHRVEasy<-function(rrs2, rrs, useWavelet = FALSE, correctSigLevel = TRUE){

  files1 = list.files(rrs2)
  class1 <- strsplit(rrs2, "/")[[1]][2]
  files2 = list.files(rrs)
  class2 <- strsplit(rrs, "/")[[1]][2]
  dataFrameMTime1 = time_analysis(files1, class1, rrs2)
  dataFrameMTime2 = time_analysis(files2, class2, rrs)
  # Creating a DF with both in Time
  dataFrameMTime=rbind(dataFrameMTime1, dataFrameMTime2)
  # Statistical analysis of both
  statistical_analysisTime(dataFrameMTime1, dataFrameMTime, correctSigLevel)
  statistical_analysisTime(dataFrameMTime2, dataFrameMTime, correctSigLevel)
  # FREQUENCY:
  if(useWavelet == FALSE){
    dataFrameMFreq1 = freq_analysis(files1, class1, rrs2)
    dataFrameMFreq2 = freq_analysis(files2, class2, rrs)
    dataFrameMFreq=rbind(dataFrameMFreq1, dataFrameMFreq2)
    statistical_analysisFreq(dataFrameMFreq1, dataFrameMFreq, correctSigLevel)
    statistical_analysisFreq(dataFrameMFreq2, dataFrameMFreq, correctSigLevel)
  } 
  # WAVELET
  if(useWavelet == TRUE){
    dataFrameMWavelet1 = wavelet_analysis(files1, class1, rrs2)
    dataFrameMWavelet2 = wavelet_analysis(files2, class2, rrs)
    dataFrameMWavelet=rbind(dataFrameMWavelet1, dataFrameMWavelet2)
    statistical_analysisFreq(dataFrameMWavelet1, dataFrameMWavelet, correctSigLevel)
    statistical_analysisFreq(dataFrameMWavelet2, dataFrameMWavelet, correctSigLevel)
  } 
}

RHRVEasy ("rrs/normal/", "rrs/chf/", useWavelet=FALSE, correctSigLevel=TRUE)   
