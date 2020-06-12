source("utils.R")

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
dataFrameMNormalWavelet=data.frame()
dataFrameMCHFWavelet=data.frame()
dataFrameMWavelet=data.frame()

all_filesNormal = list.files("rrs/normal/")
all_filesCHF= list.files("rrs/chf/")

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

# WAVELET ANALYSIS
dataFrameMNormalWavelet = wavelet_analysis(all_filesNormal, class = "NORMAL", rrs2 =  "rrs/normal")
dataFrameMCHFWavelet = wavelet_analysis(all_filesCHF, class = "CHF", rrs2 = "rrs/chf")
# Creating a DF with both CHF and Normal in Wavelet
dataFrameMWavelet=rbind(dataFrameMNormalWavelet, dataFrameMCHFWavelet)

# ANALYZING NORMALITY
shapiroNormalFreq = shapirofreq(dataFrameMNormalF)
shapiroCHFFreq = shapirofreq(dataFrameMCHF_F)
shapiroNormalTime = shapirotime(dataFrameMNormal)
shapiroCHFTime = shapirotime(dataFrameMCHF)
shapiroNormalWavelet = shapirofreq(dataFrameMNormalWavelet)
shapiroCHFWavelet = shapirofreq(dataFrameMCHFWavelet)

# KRUSKAL WALLIS
kruskalfreq(dataFrameMFrec)
kruskalfreq(dataFrameMWavelet)
kruskaltime(dataFrameMTime)

# POST HOC
dunnfreq(dataFrameMFrec)
dunnfreq(dataFrameMWavelet)
dunntime(dataFrameMTime)




