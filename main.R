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



# ANALYZING NORMALITY
shapirofreq(dataFrameMNormalF)
shapirofreq(dataFrameMCHF_F)
shapirotime(dataFrameMNormal)
shapirotime(dataFrameMCHF)

# KRUSKAL WALLIS
kruskalfreq(dataFrameMFrec)
kruskaltime(dataFrameMTime)

# POST HOC
dunnfreq(dataFrameMFrec, "bonf")
dunntime(dataFrameMTime, "bonf")




