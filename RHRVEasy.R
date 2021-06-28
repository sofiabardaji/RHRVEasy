#install.packages("RHRV")
library(RHRV)

# Post hoc Dunn test
library(dunn.test)
library(FSA)
library(PMCMR)

source('scaling_region_estimation.R')

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
  
  hrv.data = LoadBeat(fileType = format, HRVData = hrv.data,  Recordname = file, RecordPath = rrs)
  
  hrv.data=BuildNIHR(hrv.data)
  hrv.data=FilterNIHR(hrv.data)
  hrv.data$Beat = hrv.data$Beat[2: nrow(hrv.data$Beat),]
  hrv.data
}

#Calls an RHRV function with hrv.data after cleaning the parameters
easy_call <- function(hrv.data, mf, ...) {
  args.list = plotrix::clean.args(list(...), mf)
  args.list$HRVData = hrv.data
  do.call(mf, args.list)
}

# Creating time analysis data frames
time_analysis<-function(format, files, class, rrs2, ...){
  dataFrame = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, CreateTimeAnalysis, ...)
    results=hrv.data$TimeAnalysis[]
    name_file = list ("filename" = file)
    group = list ("group" = class)
    # group_name = list("group" = group)
    row_list = c (name_file, results, group)
    df=as.data.frame(row_list)
    dataFrame=rbind(dataFrame, df)
  }
  dataFrame
}

# Frequency analysis
freq_analysis<-function(format, files, class, rrs2, ...){
  dataFrame = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    hrv.data = easy_call(hrv.data, CreateFreqAnalysis, ...)
    hrv.data = easy_call(hrv.data, CalculatePSD, doPlot = F, ...)
    name_file = list ("filename" = file)
    x1 = easy_call(hrv.data, CalculateEnergyInPSDBands, ...)
    names(x1) = c("ULF", "VLF", "LF", "HF")
    group = list ("group" = class)
    row_list = c (name_file, x1, group)
    df = data.frame()
    df = rbind(df, as.data.frame(row_list))
    dataFrame=rbind(dataFrame, df)
  }
  dataFrame
}

#  Wavelet analysis
wavelet_analysis<-function(format, files, class, rrs2, ...){
  dataFrameMWavelet = data.frame()
  for (file in files) {
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
    zero_indexes = which(hrv.data$HR == 0)
    hr_median = median(hrv.data$HR[-zero_indexes])
    hrv.data$HR[zero_indexes] = hr_median
    
    hrv.data = easy_call(hrv.data, CreateFreqAnalysis, ...)
    hrv.data=SetVerbose(hrv.data, FALSE)
    hrv.data = easy_call(hrv.data, CalculatePowerBand, ...)
    
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
    name_file = list ()
    x1 = as.list(resultsWavelet)
    group = list ("group" = class)
    row_list = c (name_file, x1, group)
    dataFrameMWavelet = rbind(dataFrameMWavelet, as.data.frame(row_list))
    
  }
  dataFrameMWavelet
}

attempToCalculateTimeLag <- function(hrv.data) {
  kTimeLag = NA
  tryCatch(
    {
      kTimeLag=CalculateTimeLag(hrv.data,technique = "acf", method = "first.minimum",
                                lagMax = 20, doPlot=FALSE)
    },
    error=function(cond) {
      tryCatch(
        {
          kTimeLag=CalculateTimeLag(hrv.data,technique = "acf", method = "first.e.decay",
                                    lagMax = 20, doPlot=FALSE)
        },
        error=function(cond) {
          
          tryCatch(
            {
              kTimeLag=CalculateTimeLag(hrv.data,technique = "ami", method = "first.minimum",
                                        lagMax = 20, doPlot=FALSE)
            },
            error=function(cond) {
              tryCatch(
                {
                  kTimeLag=CalculateTimeLag(hrv.data,technique = "ami", method = "first.e.decay",
                                            lagMax = 20, doPlot=FALSE)
                },
                error=function(cond) {
                  kTimeLag=NA
                }
              )
            }
          )
        }
      )
    }
  )
  kTimeLag
}

# Non Linear analysis
non_linear_analysis <- function(format, files, class, rrs2, ...){
  dataFrame = data.frame()
  for (file in files){
    hrv.data = preparing_analysis(format, file = file, rrs = rrs2)
    hrv.data = CreateNonLinearAnalysis(hrv.data)
    kTimeLag=attempToCalculateTimeLag(hrv.data)
    
    if(is.na(kTimeLag)){
      hrv.data$NonLinearAnalysis[[1]]$correlation$statistic = NA
      hrv.data$NonLinearAnalysis[[1]]$sampleEntropy$statistic = NA
      hrv.data$NonLinearAnalysis[[1]]$lyapunov$statistic = NA   
    }
    else{
      
      tryCatch(
        {
          kEmbeddingDim = CalculateEmbeddingDim(hrv.data, numberPoints = 10000,
                                                timeLag = kTimeLag, maxEmbeddingDim = 15, doPlot=FALSE)
          if(kEmbeddingDim == 0){
            hrv.data$NonLinearAnalysis[[1]]$correlation$statistic = NA
            hrv.data$NonLinearAnalysis[[1]]$sampleEntropy$statistic = NA
            hrv.data$NonLinearAnalysis[[1]]$lyapunov$statistic = NA   
          }
          else{
            
            hrv.data = CalculateCorrDim(hrv.data, indexNonLinearAnalysis = 1, minEmbeddingDim=kEmbeddingDim,
                                        maxEmbeddingDim = kEmbeddingDim + 2, timeLag = 1, minRadius = 1, maxRadius = 50,
                                        pointsRadius = 100, theilerWindow =10, corrOrder = 2, doPlot = FALSE)
            
            cd = hrv.data$NonLinearAnalysis[[1]]$correlation$computations
            
            filteredCd = nltsFilter(cd, threshold = 0.99)
            
            cdScalingRegion = 
              estimate_scaling_region(filteredCd, numberOfLinearRegions = 3, doPlot = FALSE)
            
            
            hrv.data = EstimateCorrDim(hrv.data, indexNonLinearAnalysis=1, regressionRange=cdScalingRegion, # ahi va variable
                                       useEmbeddings=(kEmbeddingDim):(kEmbeddingDim+2), 
                                       doPlot=FALSE)
            
            hrv.data = CalculateSampleEntropy(hrv.data, indexNonLinearAnalysis= 1, doPlot = FALSE)
            
            hrv.data = EstimateSampleEntropy(hrv.data, indexNonLinearAnalysis=1, doPlot = FALSE)
            
            hrv.data = CalculateMaxLyapunov(hrv.data, indexNonLinearAnalysis = 1,
                                            minEmbeddingDim= kEmbeddingDim, 
                                            maxEmbeddingDim= kEmbeddingDim+2,
                                            timeLag = kTimeLag,radius = 3, theilerWindow = 20,
                                            doPlot = TRUE)
            
            hrv.data = EstimateMaxLyapunov(hrv.data, indexNonLinearAnalysis = 1,
                                           regressionRange = c(1,6),
                                           useEmbeddings = (kEmbeddingDim):(kEmbeddingDim+2),
                                           doPlot = TRUE)   
          }
        },
        error=function(cond) {
          message("There has been a problem calculating some non lineal statystic  !!!!!!!!")
          message(cond)
          
        })
      
    }
    resultsCS = list("CorrelationStatistic" = mean(hrv.data$NonLinearAnalysis[[1]]$correlation$statistic, na.rm = TRUE))
    resultsSE = list("SampleEntropy" = mean(hrv.data$NonLinearAnalysis[[1]]$sampleEntropy$statistic, na.rm = TRUE))
    resultsML = list("MaxLyapunov" = mean(hrv.data$NonLinearAnalysis[[1]]$lyapunov$statistic, na.rm = TRUE))
    
    #as.data.frame considers that if the value of a list is NULL it does not exist. It must contain NA
    message(c("\nresultsCS ",resultsCS["CorrelationStatistic"]))
    if(is.null(resultsCS["CorrelationStatistic"])){
      resultsCS["CorrelationStatistic"] = NA
    }
    message(c("resultsSE ",resultsSE["SampleEntropy"]))
    if(is.null(resultsSE["SampleEntropy"])){
      resultsSE["SampleEntropy"] = NA
    }
    
    message(c("resultsML ",resultsML["MaxLyapunov"]))
    if(is.null(resultsML["MaxLyapunov"])){
      resultsML["MaxLyapunov"] = NA
    }
    
    name_file = list ("filename" = file)
    group = list ("group" = class)
    row_list = c (name_file, resultsCS, resultsSE, resultsML, group)
    df=as.data.frame(row_list)
    dataFrame=rbind(dataFrame, df)
  }
  dataFrame
}


# Statistical tests application
dunnNonLinar<-function(dfM, method){
  dfM$group = factor(dfM$group)
  CorrelationStatistic = NA
  SampleEntropy = NA
  MaxLyapunov  = NA
  
  tryCatch(
    {
      CorrelationStatistic = posthoc.kruskal.dunn.test(CorrelationStatistic ~ group, data=dfM, p.adjust=method)
    },
    error=function(cond) {
      #TODO : habrá que tener en cuenta si es NULL; no aquí sino donde corresponda
      
    })
  
  tryCatch(
    {
      SampleEntropy = posthoc.kruskal.dunn.test(SampleEntropy ~ group, data=dfM, p.adjust=method)
    },
    error=function(cond) {
      #TODO : habrá que tener en cuenta si es NULL; no aquí sino donde corresponda
    })
  
  tryCatch(
    {
      MaxLyapunov = posthoc.kruskal.dunn.test(MaxLyapunov ~ group, data=dfM, p.adjust=method)
    },
    error=function(cond) {
      #TODO : habrá que tener en cuenta si es NULL; no aquí sino donde corresponda
    })
  
  list (CorrelationStatistic, SampleEntropy, MaxLyapunov)
}

dunnfreq<-function(dfM, method){
  dfM$group = factor(dfM$group)
  list (ULF = posthoc.kruskal.dunn.test(ULF ~ group, data=dfM, p.adjust=method),
        VLF = posthoc.kruskal.dunn.test(VLF ~ group, data=dfM, p.adjust=method),
        LF = posthoc.kruskal.dunn.test(LF ~ group, data=dfM, p.adjust=method),
        HF = posthoc.kruskal.dunn.test(HF ~ group, data=dfM, p.adjust=method) )
}

dunntime<-function(dfM, method){
  dfM$group = factor(dfM$group)
  list (SDNN= posthoc.kruskal.dunn.test(SDNN ~ group, data = dfM, p.adjust=method),
        SDANN = posthoc.kruskal.dunn.test(SDANN ~ group, data = dfM, p.adjust=method),
        SDNNIDX = posthoc.kruskal.dunn.test(SDNNIDX ~ group, data = dfM, p.adjust=method),
        pNN50 = posthoc.kruskal.dunn.test(pNN50 ~ group, data = dfM, p.adjust=method),
        SDSD = posthoc.kruskal.dunn.test(SDSD ~ group, data = dfM, p.adjust=method),
        rMSSD = posthoc.kruskal.dunn.test(rMSSD ~ group, data = dfM, p.adjust=method),
        IRRR = posthoc.kruskal.dunn.test(IRRR ~ group, data = dfM, p.adjust=method),
        MADRR = posthoc.kruskal.dunn.test(MADRR ~ group, data = dfM, p.adjust=method),
        TINN = posthoc.kruskal.dunn.test(TINN ~ group, data = dfM, p.adjust=method),
        HRVi = posthoc.kruskal.dunn.test(HRVi ~ group, data = dfM, p.adjust=method))
}

shapiro.test.CheckAllVakuesEqual<-function(x){
  pval = 1 # If we cannot do the test, I see because all the numerical values are the same, we will return 1 since there are no differences between the populations.
  tryCatch(
    {
      pval = shapiro.test(x)$p.value
    },
    error=function(cond) {
      message("All values identical in shapiro.test; pvalue set to 1")
    }
  )
  pval
}



statistical_analysisFreq<-function(dfM, verbose, numberOfExperimentalGroups, method, signif_level){
  anova = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  kruskal = list(ULF = NA, VLF = NA, LF = NA, HF = NA)
  dunn = NA
  list = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  vec = c("group" = NA, "p-value ULF" = NA, "p-value VLF" = NA, "p-value LF" = NA, 
          "p-value HF" = NA)
  
  
  for(objeto in names(listDF)){
    vec$group = objeto
    
    for (column in c('ULF', 'VLF', 'LF', 'HF')){
      destino = paste0('p-value ', column)
      vec[[destino]] = shapiro.test.CheckAllVakuesEqual(listDF[[objeto]][[column]])
    }
    
    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
  }
  
  for (column in c('ULF', 'VLF', 'LF', 'HF')){
    p_values = formula_str = paste0("p.value.", column)
    formula_str = paste0(column, "~ group")
    formula = as.formula(formula_str)
    
    if (numberOfExperimentalGroups > 2 || all(dataFramePvalues[[p_values]] > signif_level)) {
      if (verbose == TRUE){
        cat(column, " Normal: Anova. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      list$anova[[column]] = aov(formula, data = dfM)
    }else {
      if (verbose == TRUE){
        cat(column, " NOT normal: Kruskal. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      list$kruskal[[column]] = kruskal.test(formula, data = dfM)
    }
    
  }
  
  list$dunn = dunnfreq(dfM, method)
  list
  
}

statistical_analysisTime<-function(dfM, verbose, numberOfExperimentalGroups, method, signif_level){
  
  anova = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
               MADRR = NA, TINN = NA, HRVi = NA)
  kruskal = list(SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                 MADRR = NA, TINN = NA, HRVi = NA)
  dunn = NA
  list = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  
  vec = list("group" = NA, "p-value SDNN" = NA, "p-value SDANN" = NA,
             "p-value SDNNIDX" = NA, "p-value pNN50" = NA, "p-value SDSD" = NA, 
             "p-value rMSSD" = NA, "p-value IRRR" = NA, "p-value MADRR" = NA, 
             "p-value TINN" = NA, "p-value HRVi" = NA)
  
  for(objeto in names(listDF)){
    vec$group = objeto
    
    for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                     'MADRR', 'TINN', 'HRVi')){
      destino = paste0('p-value ', column)
      vec[[destino]] =  shapiro.test.CheckAllVakuesEqual(listDF[[objeto]][[column]])
    }
    
    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
  }
  
  for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    p_values = formula_str = paste0("p.value.", column)
    formula_str = paste0(column, "~ group")
    formula = as.formula(formula_str)
    
    if (numberOfExperimentalGroups > 2 || all(dataFramePvalues[[p_values]] > signif_level)) {
      if (verbose == TRUE){
        cat(column, " Normal: Anova. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      list$anova[[column]] = aov(formula, data = dfM)
    }else {
      if (verbose == TRUE){
        cat(column, " NOT normal: Kruskal. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      list$kruskal[[column]] = kruskal.test(formula, data = dfM)
    }
  }
  
  list$dunn = dunntime(dfM, method)
  list
}

statistical_analysisNonLinear<-function(dfM, verbose, numberOfExperimentalGroups, method, signif_level){
  anova = list(CorrelationStatistic = NA, SampleEntropy = NA, MaxLyapunov = NA)
  kruskal =list(CorrelationStatistic = NA, SampleEntropy = NA, MaxLyapunov = NA)
  dunn = NA
  list = list(anova = anova, kruskal = kruskal, dunn = dunn)
  
  listDF = split(dfM, dfM$group)
  
  dataFramePvalues = data.frame()
  vec = c("group" = NA, "p-value CorrelationStatistic" = NA, "p-value SampleEntropy" = NA,  
          "p-value MaxLyapunov" = NA)
  
  for(objeto in names(listDF)){
    vec$group = objeto
    
    for (column in c('CorrelationStatistic', 'SampleEntropy', 'MaxLyapunov')){
      destino = paste0('p-value ', column)
      vec[[destino]] =  shapiro.test.CheckAllVakuesEqual(listDF[[objeto]][[column]])
    }
    
    df = data.frame(vec)
    
    dataFramePvalues = rbind(dataFramePvalues, df)
  }
  
  for (column in c('CorrelationStatistic', 'SampleEntropy', 'MaxLyapunov')){
    p_values = formula_str = paste0("p.value.", column)
    formula_str = paste0(column, "~ group")
    formula = as.formula(formula_str)
    if (numberOfExperimentalGroups > 2 || all(dataFramePvalues[[p_values]] > signif_level)) {
      if (verbose == TRUE){
        cat(column, " Normal: Anova. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      
      #ANOVA will fail if the statistic could not be calculated for all recordings in a group
      tryCatch(
        {
          list$anova[[column]] = aov(formula, data = dfM)
        },
        error=function(cond) {
          #TODO : habrá que tener en cuenta si es NULL
          list$anova[[column]] = NULL
        })
      
      
    }else {
      if (verbose == TRUE){
        cat(column, " NOT normal: Kruskal. P-values = ", dataFramePvalues[[p_values]], "\n")
      }
      #Krustal will fail if the statistic could not be calculated for all recordings in a group
      tryCatch(
        {      
          list$kruskal[[column]] = kruskal.test(formula, data = dfM)
        },
        error=function(cond) {
          #TODO : habrá que tener en cuenta si es NULL
          list$kruskal[[column]] = NULL
        })
    }
  }
  list$dunn = dunnNonLinar(dfM, method)
  list
  
}


correctpValues <- function(listTime, listFreq, listNonLinear, correction, method){
  
  listpValues = list(ULF = NA, VLF = NA, LF = NA, HF = NA,
                     SDNN = NA, SDANN = NA, SDNNIDX = NA, pNN50 = NA, SDSD = NA, 
                     rMSSD = NA, IRRR = NA, MADRR = NA, TINN = NA, HRVi = NA, 
                     CorrelationStatistic = NA, SampleEntropy = NA, MaxLyapunov = NA)
  
  listpValuesCorrected = list(ULF = NA, VLF = NA, LF = NA, HF = NA, SDNN = NA, SDANN = NA, 
                              SDNNIDX = NA, pNN50 = NA, SDSD = NA, rMSSD = NA, IRRR = NA,
                              MADRR = NA, TINN = NA, HRVi = NA,
                              CorrelationStatistic = NA, SampleEntropy = NA, MaxLyapunov = NA)
  
  for (column in c('ULF', 'VLF', 'LF', 'HF')){
    if(is.na(listFreq$anova[[column]])){
      listpValues[[column]] = listFreq$kruskal[[column]]$p.value
    }else{
      listpValues[[column]] = extract_ANOVA_pvalue(listFreq$anova[[column]])
    }
  }
  
  for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    if(is.na(listTime$anova[[column]])){
      listpValues[[column]] = listTime$kruskal[[column]]$p.value
    }else{
      listpValues[[column]] = extract_ANOVA_pvalue(listTime$anova[[column]])
    }
  }
  
  # In order for it to only be performed when there is non linear results: 
  if(!is.na(listNonLinear)){
    for (column in c('CorrelationStatistic', 'SampleEntropy', 'MaxLyapunov')){
      if(is.na(listNonLinear[["anova"]][[column]])){
        listpValues[[column]] = listNonLinear[["kruskal"]][[column]][["p.value"]]
      }else{
        listpValues[[column]] = extract_ANOVA_pvalue(listNonLinear[["anova"]][[column]])
      }
    }
  }
  
  if (correction == TRUE){
    listpValuesCorrected = p.adjust(listpValues, method)
    listpValuesCorrected <- as.list(listpValuesCorrected)
    
  }else if (correction == FALSE){
    listpValuesCorrected = listpValues
  }
  listpValuesCorrected
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
  
  listDF = split(results$TimeAnalysis, results$TimeAnalysis$group)
  
  differencesFound = FALSE
  
  cat("\n\nResult of the analysis of the variability of the heart rate of the group",
      levels(results$TimeAnalysis$group)[1], 
      "versus the group", levels(results$TimeAnalysis$group)[2], ":\n\n")
  
  for (column in c('SDNN', 'SDANN', 'SDNNIDX', 'pNN50', 'SDSD', 'rMSSD', 'IRRR',
                   'MADRR', 'TINN', 'HRVi')){
    if(all(is.na(results$StatysticalAnalysisTime$anova[[column]]))){
      #report kruskal
      if(results$pValues[[column]]<signif_level){#error pvalue 1
        differencesFound = TRUE
        cat("\nThere is a statistically significant difference in", column,  "; pvalue: ", 
            results$StatysticalAnalysisTime$kruskal[[column]]$p.value, "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
      }
    }
    #report anova
    else{
      if(results$pValues[[column]]<signif_level){
        differencesFound = TRUE
        cat("\nThere is a statistically significant difference in", column, "; pvalue: ", 
            extract_ANOVA_pvalue(results$StatysticalAnalysisTime$anova[[column]]), "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group ", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
        
      }
    }
    
    
    # Two Conditions to report Dunn:
    # 1. We have more than 2 groups, we check that by looking at the length of listDF
    if(length(listDF)>2){
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
      var = which(results$StatysticalAnalysisTime$dunn[[column]][["p.value"]]<signif_level, arr.ind = TRUE)
      if(length(var)>0){
        cat("The groups with stastically significant differences are:")
        for (i in 1:nrow(var)){
          cat('\n',results$StatysticalAnalysisTime$dunn[[column]][["p.value"]][var][i], 'is the p value of',
              row.names(results$StatysticalAnalysisTime$dunn[[column]][["p.value"]])[var[i,1]], 'vs',
              colnames(results$StatysticalAnalysisTime$dunn[[column]][["p.value"]])[var[i,2]], '')
        }
        cat('for', column, '\n')
        
      }
    }
    
  }
  
  listDF = split(results$FrequencyAnalysis, results$FrequencyAnalysis$group)
  
  cat("\n\nFrequency analysis:\n")
  
  for (column in c('ULF', 'VLF', 'LF', 'HF')){
    if(all(is.na(results$StatysticalAnalysisFrequency$anova[[column]]))){
      #report kruskal
      if(results$pValues[[column]]<signif_level){#error pvalue 1
        differencesFound = TRUE
        cat("\nThere is a statistically significant difference in", column,  "; pvalue: ", 
            results$StatysticalAnalysisFrequency$kruskal[[column]]$p.value, "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
      }
    }
    #report anova
    else{
      if(results$pValues[[column]]<signif_level){
        differencesFound = TRUE
        cat("\nThere is a statistically significant difference in", column, "; pvalue: ", 
            extract_ANOVA_pvalue(results$StatysticalAnalysisFrequency$anova[[column]]), "\n")
        
        for (i in 1:length(listDF)){
          group = levels(results$TimeAnalysis$group)[i]
          cat(column, " for the group ", levels(results$TimeAnalysis$group)[i], "is",
              mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
        }
        
      }
    }
    
    # Two Conditions to report Dunn:
    # 1. We have more than 2 groups, we check that by looking at the length of listDF
    
    if(length(listDF)>2){
      
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
      
      variable = which(results$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]]<signif_level, 
                       arr.ind = TRUE)
      if(length(variable)>0){
        cat("The groups with stastically significant differences are:")
        for (i in 1:nrow(variable)){
          cat('\n',results$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]][variable][i], 
              'is the p value of',
              row.names(results$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]])[variable[i,1]], 
              'vs',
              colnames(results$StatysticalAnalysisFrequency$dunn[[column]]$p.value)[variable[i,2]], '')
        }
        cat('for', column, '\n')
      }
      
    }
    
  }
  
  if(!all(is.na(results$NonLinearAnalysis))){
    listDF = split(results$NonLinearAnalysis, results$NonLinearAnalysis$group)
    
    cat("\n\nNon Linear analysis:\n")
    
    for (column in c('CorrelationStatistic', 'SampleEntropy', 'MaxLyapunov')){
      if(all(is.na(results$StatysticalAnalysisNonLinear$anova[[column]]))){
        #report kruskal
        if(results$pValues[[column]]<signif_level){#error pvalue 1
          differencesFound = TRUE
          cat("\nThere is a statistically significant difference in", column,  "; pvalue: ", 
              results$StatysticalAnalysisNonLinear$kruskal[[column]]$p.value, "\n")
          
          for (i in 1:length(listDF)){
            group = levels(results$TimeAnalysis$group)[i]
            cat(column, " for the group", levels(results$TimeAnalysis$group)[i], "is",
                mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
          }
        }
      }
      #report anova
      else{
        if(results$pValues[[column]]<signif_level){
          differencesFound = TRUE
          cat("\nThere is a statistically significant difference in", column, "; pvalue: ", 
              extract_ANOVA_pvalue(results$StatysticalAnalysisNonLinear$anova[[column]]), "\n")
          
          for (i in 1:length(listDF)){
            group = levels(results$TimeAnalysis$group)[i]
            cat(column, " for the group ", levels(results$TimeAnalysis$group)[i], "is",
                mean(listDF[[group]][[column]]), "+-", sd(listDF[[group]][[column]]), "\n")
          }
          
        }
      }
      
      # Two Conditions to report Dunn:
      # 1. We have more than 2 groups, we check that by looking at the length of listDF
      
      if(length(listDF)>2){
        
        # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
        
        variable = which(results$StatysticalAnalysisNonLinear$dunn[[column]][["p.value"]]<signif_level, 
                         arr.ind = TRUE)
        if(length(variable)>0){
          cat("The groups with stastically significant differences are:")
          for (i in 1:nrow(variable)){
            cat('\n',results$StatysticalAnalysisNonLinear$dunn[[column]][["p.value"]][variable][i], 
                'is the p value of',
                row.names(results$StatysticalAnalysisNonLinear$dunn[[column]][["p.value"]])[variable[i,1]], 
                'vs',
                colnames(results$StatysticalAnalysisNonLinear$dunn[[column]]$p.value)[variable[i,2]], '')
          }
          cat('for', column, '\n')
        }
        
      }
      
    }
  }

  if(!differencesFound){
    cat("No statistically significant difference were found\n")
  }
}


RHRVEasy<-function(folders, correction = FALSE, method = "bonferroni", verbose=FALSE, 
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
  
  files = list()
  
  for (folder in folders){
    file_validation(folder)
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format,
                                                                          list.files(folder), split_path(folder)[1], folder, ...))
    if(nonLinear == TRUE){
      dataFrameMNonLinear = rbind(dataFrameMNonLinear,non_linear_analysis(format, 
                                                                          list.files(folder), split_path(folder)[1], folder, ...))
      
    }
  }
  
  numberOfExperimentalGroups = length(folders)
  # Statistical analysis of both
  
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime,
                                                         verbose, numberOfExperimentalGroups, method, signif_level)
  
  # FREQUENCY:
  if(typeAnalysis == "fourier"){
    for (folder in folders){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format,
                                                                            list.files(folder), split_path(folder)[1], folder, ...))
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq,
                                                           verbose, numberOfExperimentalGroups, method, signif_level)
  }
  
  # WAVELET
  if(typeAnalysis == "wavelet"){
    for (folder in folders){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format,
                                                                                        list.files(folder), split_path(folder)[1], folder,
                                                                                        type = typeAnalysis, ...))
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet,
                                                           verbose, numberOfExperimentalGroups, method, signif_level)
    
    dataFrameMFreq = dataFrameMWavelet
  }
  
 
  if(!all(is.na(dataFrameMNonLinear))){
    listNonLinearStatisticalAnalysis = statistical_analysisNonLinear(dataFrameMNonLinear,
                                                                     verbose, numberOfExperimentalGroups, method, signif_level)
    listpValues = correctpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis,
                                 listNonLinearStatisticalAnalysis,
                                 correction, method)
  }else{
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
