#install.packages("RHRV")
library(RHRV)



dataFrameMNormal=data.frame()
dataFrameMFrec=data.frame()
dataFrameMCHF=data.frame()
dataFrame2=data.frame()
dataFrame3=data.frame()
df=data.frame()
df2=data.frame()
df_grande=data.frame()
dataFrameMNormalF=data.frame()
dataFrameMCHF_F=data.frame()
all_filesNormal = list.files("rrs/normal/")
all_filesCHF= list.files("rrs/chf/")

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
    hrv.data = preparing_analysis(directory_files = all_filesNormal, file = file, rrs = rrs2)
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
    hrv.data = preparing_analysis(directory_files = all_filesNormal, file = file, rrs = rrs2)
    # PlotNIHR(hrv.data)
    hrv.data=InterpolateNIHR(hrv.data)
    # Find the zeros in HR
    zero_indexes = which(hrv.data$HR == 0)
    # Compute the median of HR after removing 0s
    hr_median = median(hrv.data$HR[-zero_indexes])
    # Finally, substitute the 0s in HR by the median
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

# TIME ANALYSIS: 
dataFrameMNormal = time_analysis(files = all_filesNormal, class = "NORMAL", rrs2 = "rrs/normal")
dataFrameMCHF = time_analysis(files = all_filesCHF, class = "CHF", rrs2 = "rrs/chf")
# Creating a DF with both CHF and Normal in Time
dataFrameMTime=rbind(dataFrameMCHF, dataFrameMNormal)


# FREQUENCY ANALYSIS
dataFrameMNormalF = freq_analysis(files = all_filesNormal, class = "NORMAL", rrs2 = "rrs/normal")
dataFrameMCHF_F = freq_analysis(files = all_filesCHF, class = "CHF", rrs2 = "rrs/chf")
# Creating a DF with both CHF and Normal in Freq
dataFrameMFrec=rbind(dataFrameMCHF_F, dataFrameMNormalF)


# ANALYZING NORMALITY
shapirofreq<-function(directorio){
  list(shapiro.test(directorio$ULF), shapiro.test(directorio$VLF), shapiro.test(directorio$LF))
}
shapirotime<-function(directorio){
  list(shapiro.test(directorio$SDNN), shapiro.test(directorio$SDANN), shapiro.test(directorio$SDNNIDX),
       shapiro.test(directorio$pNN50),shapiro.test(directorio$SDSD),shapiro.test(directorio$rMSSD),
       shapiro.test(directorio$IRRR),shapiro.test(directorio$MADRR),shapiro.test(directorio$TINN),
       shapiro.test(directorio$HRVi))
}

shapirofreq(directorio = dataFrameMNormalF)
shapirofreq(directorio = dataFrameMCHF)
shapirotime(dataFrameMNormal)
shapirotime(dataFrameMCHF)


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

kruskalfreq(dataFrameMFrec)
kruskaltime(dataFrameMTime)

# POST HOC DUNN TEST
#install.packages("PMCMR)
#install.packages('FSA')
#install.packages("dunn.test")
library(dunn.test)
#install.packages("FSA")
library(FSA)
#install.packages("PMCMR")
library(PMCMR)

# POSTHOC DUNN
dunnfreq<-function(df, adj){
  list (posthoc.kruskal.dunn.test(ULF ~ clase, data=df, p.adjust=adj),
  posthoc.kruskal.dunn.test(VLF ~ clase, data=df, p.adjust=adj),
  posthoc.kruskal.dunn.test(LF ~ clase, data=df, p.adjust=adj),
  posthoc.kruskal.dunn.test(HF ~ clase, data=df, p.adjust=adj) )
}
dunntime<-function(df, adj){
  list (posthoc.kruskal.dunn.test(SDNN ~ clase, data = df, p.adjust=adj), 
        posthoc.kruskal.dunn.test(SDANN ~ clase, data = df, p.adjust=adj), 
        posthoc.kruskal.dunn.test(SDNNIDX ~ clase, data = df, p.adjust=adj), 
        posthoc.kruskal.dunn.test(pNN50 ~ clase, data = df, p.adjust=adj),
        posthoc.kruskal.dunn.test(SDSD ~ clase, data = df, p.adjust=adj), 
        posthoc.kruskal.dunn.test(rMSSD ~ clase, data = df, p.adjust=adj),
        posthoc.kruskal.dunn.test(IRRR ~ clase, data = df, p.adjust=adj),
        posthoc.kruskal.dunn.test(MADRR ~ clase, data = df, p.adjust=adj),
        posthoc.kruskal.dunn.test(TINN ~ clase, data = df, p.adjust=adj), 
        posthoc.kruskal.dunn.test(HRVi ~ clase, data = df, p.adjust=adj))
}
dunnfreq(dataFrameMFrec, "bonf")
dunntime(dataFrameMTime, "bonf")



