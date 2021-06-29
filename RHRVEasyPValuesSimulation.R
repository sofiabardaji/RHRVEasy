

RHRVEasyPValues<-function(folders, correction = FALSE, method = "bonferroni", verbose=FALSE, 
                   format = "RR", typeAnalysis = 'fourier', significance_level = 0.05, nonLinear=FALSE, ...) {
  
  dataFrameMWavelet = data.frame()
  dataFrameMTime = data.frame()
  dataFrameMFreq = data.frame()
  dataFrameMNonLinear = data.frame()
  listNonLinearStatisticalAnalysis = list()
  listTimeStatysticalAnalysis = list()
  listFreqStatysticalAnalysis = list()
  
  #We create a global variable signif_level with the level significance
  signif_level <<- significance_level
  #We create a global variable verb with the verbose mode  
  verb <<- verbose
  
  files = list()
  
  f  <- vector(mode = "list", length = length(folders))
  i=1
  
  for (folder in folders){
    file_validation(folder)
    if(i==1){
      f[[i]] =  sample(list.files(folder), size =18, replace = FALSE)
    }
    else{
      f[[i]] =  sample(list.files(folder), size =10, replace = FALSE)
    }
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format,
                                                                          f[[i]], split_path(folder)[1], folder, ...))
    if(nonLinear == TRUE){
      cat("Performing non linear analysis...")
      dataFrameMNonLinear = rbind(dataFrameMNonLinear,non_linear_analysis(format, 
                                                                          f[[i]], split_path(folder)[1], folder, ...))
    }
    i = i+1
  }
  
  numberOfExperimentalGroups = length(folders)
  # Statistical analysis of both
  
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime,
                                                         numberOfExperimentalGroups, method, signif_level)
  
  # FREQUENCY:
  if(typeAnalysis == "fourier"){
    i=1
    for (folder in folders){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format,
                                                                            f[[i]], split_path(folder)[1], folder, ...))
      i=i+1
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq,
                                                           numberOfExperimentalGroups, method, signif_level)
  }
  
  # WAVELET
  if(typeAnalysis == "wavelet"){
    i=1
    for (folder in folders){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format,
                                                                                        f[[i]], split_path(folder)[1], folder,
                                                                                        type = typeAnalysis, ...))
      i=i+1
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet,
                                                           numberOfExperimentalGroups, method, signif_level)
    
    dataFrameMFreq = dataFrameMWavelet
  }
  
  
  if(!all(is.na(dataFrameMNonLinear))){
    #cat("DF non linear NOT empty: correct p values with those values")
    listNonLinearStatisticalAnalysis = statistical_analysisNonLinear(dataFrameMNonLinear,
                                                                     numberOfExperimentalGroups, method, signif_level)
    listpValues = correctpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis,
                                 listNonLinearStatisticalAnalysis,
                                 correction, method)
  }else{
    #cat("DF non linear empty")
    listNonLinearStatisticalAnalysis = NA
    listpValues = correctpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis,
                                 listNonLinearStatisticalAnalysis,
                                 correction, method)
  }
  
  
  
  
  results = list("TimeAnalysis" = dataFrameMTime, "StatysticalAnalysisTime" = listTimeStatysticalAnalysis,
                 "FrequencyAnalysis" = dataFrameMFreq, "NonLinearAnalysis" = dataFrameMNonLinear,
                 "StatysticalAnalysisNonLinear" = listNonLinearStatisticalAnalysis,
                 "StatysticalAnalysisFrequency" = listFreqStatysticalAnalysis, "pValues" = listpValues)
  
  class(results) = "RHRVEasyResult"
  results
}
