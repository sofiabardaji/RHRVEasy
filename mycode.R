
c1 = list(ULF = 0, VLF = 0, LF = 0, HF = 0,
             SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, 
             rMSSD = 0, IRRR = 0, MADRR = 0, TINN = 0, HRVi = 0)


c2 = c1
c3 = c1
c4 = c1
c5 = c1
c6 = c1
l = list("none"=c1, "bonferroni"=c2, "holm"=c3, "hochberg"=c4, "fdr"=c5, "BY"=c6)


columNames = c('ULF', 'VLF', 'LF', 'HF', 'SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
               'MADRR', 'TINN', 'HRVi')
records =c("C:\\rrsc\\normal","C:\\rrsc\\chf")

corrections = c("none", "bonferroni", "holm", "hochberg", "fdr", "BY")

for(correct in corrections){
  for(i in 1:100){
    a=RHRVEasy(folders = records,correction = TRUE, method = correct )
    for (column in columNames){
      if(a$pValues[[column]] < 0.05){
        l[[correct]][[column]] = l[[correct]][[column]] +1;
      }
    }
    
  }
}

res =cbind(l$none, l$bonferroni,l$holm,l$hochberg,l$fdr,l$BY)
res




c2 = list(ULF = 0, VLF = 0, LF = 0, HF = 0,
         SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, 
         rMSSD = 0, IRRR = 0, MADRR = 0, TINN = 0, HRVi = 0)


for(i in 1:100){
  a=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE)
  
  for (column in c('ULF', 'VLF', 'LF', 'HF', 'SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    
    if(a$pValues[[column]] < 0.05){
      c2[[column]] = c2[[column]]+1;
    }
  }
  
}


c3 = list(ULF = 0, VLF = 0, LF = 0, HF = 0,
         SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, 
         rMSSD = 0, IRRR = 0, MADRR = 0, TINN = 0, HRVi = 0)


for(i in 1:100){
  a=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "fdr")
  
  for (column in c('ULF', 'VLF', 'LF', 'HF', 'SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    
    if(a$pValues[[column]] < 0.05){
      c3[[column]] = c3[[column]]+1;
    }
  }
  
}

con12_c=cbind(c1, c2,c3)

#------------------



RHRVEasy<-function(folders, correction = FALSE, method = "bonferroni", verbose=FALSE, 
                   format = "RR", typeAnalysis = 'fourier', significance_level = 0.05, ...) {
  dataFrameMWavelet = data.frame()
  dataFrameMTime = data.frame()
  dataFrameMFreq = data.frame()
  #We create a global variable signif_level with the level significance
  signif_level <<- significance_level
  
  files = list()
  
  f  <- vector(mode = "list", length = length(folders))
  i=1
  
  for (folder in folders){
    file_validation(folder)
    f[[i]]=  sample(list.files(folder), size =10, replace = FALSE)
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format, 
                                                                          f[[i]], split_path(folder)[1], folder, ...))
    
    i = i+1
  }
  
  numberOfExperimentalGroups = length(folders)
  # Statistical analysis of both
  
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime, 
                                                         verbose, numberOfExperimentalGroups)
  
  # FREQUENCY:
  if(typeAnalysis == "fourier"){
    i=1
    for (folder in folders){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format, 
                                                                            f[[i]], split_path(folder)[1], folder, ...))
      i = i+1
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
    i=1
    for (folder in folders){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format,
                                                                                        f[[i]], split_path(folder)[1], folder, 
                                                                                        type = typeAnalysis), ...)
      i = i+1
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
