createIntegrSiteObject = function(pathWorkingDirectory) {
  if(missing(pathWorkingDirectory)) pathWorkingDirectory = getwd()
  else stopifnot(is.character(pathWorkingDirectory) & (length(pathWorkingDirectory) == 1L))
  integrSiteObject = createEmptyIntegrSiteObject()
  integrSiteObject = setWorkingDirectory(integrSiteObject, pathWorkingDirectory)
  return(integrSiteObject)
}

setOrganism = function(integrSiteObject, organismType = "human", fileInfoGenes = "") {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  if(!is.character(organismType)) stop("error")
  if(length(organismType) != 1L) stop("wrong length")
  if(!any(organismType %in% c("simulation", "human", "mouse"))) stop("error")
  if(!is.character(fileInfoGenes)) stop("error")
  if(length(fileInfoGenes) != 1L) stop("error")
  
  if(organismType == "human") organism = includeOrganismHuman()
  if(organismType == "mouse") organism = includeOrganismMouse()
  if(organismType == "simulation") organism = includeOrganismSimulation()
  integrSiteObject$organism$type = organismType
  integrSiteObject$organism$infoChromosomes = organism$infoChromosomes
  integrSiteObject$organism$genome = organism$genome
  integrSiteObject$organism$fileInfoGenes = fileInfoGenes
  return(integrSiteObject)
}

setDesign = function(integrSiteObject, filesIntegrations = character(), covariatesIntegrations = data.frame(), id) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  if(!is.character(filesIntegrations)) stop("filesIntegrations must be a 'character'")
  if(!is.data.frame(covariatesIntegrations)) stop("error must be a 'data.frame'")
  if(length(covariatesIntegrations) > 0L) if(length(filesIntegrations) != nrow(covariatesIntegrations)) stop("if 'covariatesIntegrations' is not empty, its number of rows must be equal to the length of 'filesIntegrations'")
  if(missing(id)) id = seq(1, length(filesIntegrations), length.out = length(filesIntegrations))
  else stop("not yet implemented")
  
  if(length(filesIntegrations) > 0L) {
    if(length(filesIntegrations) != length(unique(filesIntegrations))) stop("error")
    integrSiteObject = setFilesIntegration(integrSiteObject, filesIntegrations, id)
    if(ncol(covariatesIntegrations) > 0L) integrSiteObject = setCovariatesIntegrations(integrSiteObject, covariatesIntegrations, id)
  }
  return(integrSiteObject)
}

StoreResults = function(integrSiteObject, iSOresult, nameFileResult = "iSiteResult.csv", nameFileTable = "iSiteTable.csv", storeTable = FALSE) {
  if(class(integrSiteObject) != "integrSiteObject") stop("argument 'integrSiteObject' must be of class 'integrSiteObject'")
  #if(class(iSOresult$result)[1] != "procTestsISObj") stop("argument 'processedTests' must be of class 'procTestsISObj'")
  if(!is.character(nameFileResult) | !is.character(nameFileTable)) stop("arguments 'nameFileResult' and 'nameFileTable' must be of class 'character'")
  if((length(nameFileResult) != 1L) | (length(nameFileTable) != 1L)) stop("arguments 'nameFileResult' and 'nameFileTable' must be of length 1")
  
  iSOresult$result = unique(iSOresult$result, by = "geneName")
  
  pwd = integrSiteObject$infoStorage$pathWorkingDirectory
  data.table :: fwrite(iSOresult$result, file = paste0(pwd, "/", nameFileResult), sep = "\t", row.names = FALSE, col.names = TRUE)
  if(storeTable) data.table :: fwrite(iSOresult$table,  file = paste0(pwd, "/", nameFileTable ), sep = "\t", row.names = FALSE, col.names = TRUE)
  return(invisible(NULL))
}


iSiteCloneFitness = function(integrSiteObject, 
                             idFiles, 
                             testedCovariate,
                             addPromoter = 0L,
                             selectedCovariates = character(), 
                             selCovAreFactors = rep(F, length(selectedCovariates)),
                             dummyForDatasets = FALSE,
                             testType = c("likRt", "score")[1], 
                             pvalueCorrectionMethod = "fdr",
                             nminSplit = 2L, 
                             numberCores = 1L, 
                             sizeChunks = 100L,
                             recursive = FALSE,
                             robust = FALSE,
                             excludeIntegrations = list(),
                             excludeInBothGroups = FALSE) {
  ags = ls()
  tA = "clfi"
  aA = setArgumentsAnalysis(tA, argumentsAnalysis = mget(ags[ags != "integrSiteObject"]))
  an = iSiteAnalysis(typeAnalysis = tA, integrSiteObject, argumentsAnalysis = aA)
  if(recursive) an$result = iSiteRecursAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  if(robust)    an$result = iSiteRobustAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  return(an)
}

iSiteCloneDiffFitness = function(integrSiteObject, 
                                 idFilesGroup1, 
                                 idFilesGroup2, 
                                 testedCovariate,
                                 addPromoter = 0L,
                                 selectedCovariates = character(), 
                                 selCovAreFactors = rep(F, length(selectedCovariates)),
                                 dummyForDatasets = FALSE,
                                 pvalueCorrectionMethod = "fdr",
                                 nminSplit = 2L, 
                                 numberCores = 1L, 
                                 sizeChunks = 100L,
                                 recursive = FALSE,
                                 robust = FALSE,
                                 excludeIntegrations = list(),
                                 excludeInBothGroups = FALSE) {
  ags = ls()
  tA = "cdif"
  aA = setArgumentsAnalysis(tA, argumentsAnalysis = mget(ags[ags != "integrSiteObject"]))
  an = iSiteAnalysis(typeAnalysis = tA, integrSiteObject, argumentsAnalysis = aA)
  if(recursive) an$result = iSiteRecursAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  if(robust)    an$result = iSiteRobustAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  return(an)
}

iSiteGeneOverTarget = function(integrSiteObject, 
                               idFiles, 
                               #outcomeVariable,
                               addPromoter = 0L,
                               selectedCovariates = character(), 
                               selCovAreFactors = rep(F, length(selectedCovariates)),
                               dummyForDatasets = FALSE,
                               typeRegressionModel = c("logistic", "poisson")[1], 
                               pvalueCorrectionMethod = "fdr",
                               nminSplit = 2L, 
                               numberCores = 1L,
                               recursive = FALSE,
                               robust = FALSE,
                               excludeIntegrations = list(),
                               excludeInBothGroups = FALSE) {
  ags = ls()
  tA = "govt"
  aA = setArgumentsAnalysis(tA, argumentsAnalysis = mget(ags[ags != "integrSiteObject"]))
  an = iSiteAnalysis(typeAnalysis = tA, integrSiteObject, argumentsAnalysis = aA)
  if(recursive) an$result = iSiteRecursAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  if(robust)    an$result = iSiteRobustAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  return(an)
}

iSiteGeneDiffTarget = function(integrSiteObject, 
                               idFilesGroup1, 
                               idFilesGroup2, 
                               #outcomeVariable,
                               addPromoter = 0L,
                               selectedCovariates = character(), 
                               selCovAreFactors = rep(F, length(selectedCovariates)),
                               dummyForDatasets = FALSE,
                               typeRegressionModel = c("logistic", "poisson")[1], 
                               pvalueCorrectionMethod = "fdr",
                               nminSplit = 2L, 
                               numberCores = 1L,
                               recursive = FALSE,
                               robust = FALSE,
                               excludeIntegrations = list(),
                               excludeInBothGroups = FALSE) {
  ags = ls()
  tA = "gdit"
  aA = setArgumentsAnalysis(tA, argumentsAnalysis = mget(ags[ags != "integrSiteObject"]))
  an = iSiteAnalysis(typeAnalysis = tA, integrSiteObject, argumentsAnalysis = aA)
  if(recursive) an$result = iSiteRecursAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  if(robust)    an$result = iSiteRobustAnalysis(typeAnalysis = tA, an$integrSiteObject, an$result, argumentsAnalysis = aA)
  return(an)
}

setArgumentsAnalysis = function(typeAnalysis, argumentsAnalysis) {
  for(n in names(argumentsAnalysis)) assign(n, argumentsAnalysis[[n]], envir = environment())
  aA = list()
  cf = typeAnalysis %in% c("clfi", "cdif")
  di = typeAnalysis %in% c("cdif", "gdit")

  if(di) {
    aA$idFilesGroup1 = idFilesGroup1
    aA$idFilesGroup2 = idFilesGroup2
  } else {
    aA$idFilesGroup1 = idFiles
    aA$idFilesGroup2 = integer()
  }
  aA$addPromoter = addPromoter
  aA$selectedCovariates = selectedCovariates
  aA$selCovAreFactors = selCovAreFactors
  aA$dummyForDatasets = dummyForDatasets
  if(cf) {
    aA$testedCovariate = testedCovariate
    if(typeAnalysis == "clfi") aA$testType = testType
    if(typeAnalysis == "cdif") aA$testType = "likRt"
  } else {
    #aA$outcomeVariable = outcomeVariable
    aA$typeRegressionModel = typeRegressionModel
  } 
  aA$pvalueCorrectionMethod = pvalueCorrectionMethod
  aA$nminSplit = nminSplit
  aA$numberCores = numberCores
  if(cf) aA$sizeChunks = sizeChunks
  aA$excludeIntegrations = excludeIntegrations
  aA$robust = robust
  aA$excludeInBothGroups = excludeInBothGroups
  return(aA)
}

iSiteAnalysis = function(typeAnalysis, integrSiteObject, argumentsAnalysis) {
  if(class(integrSiteObject) != "integrSiteObject") stop("argument 'integrSiteObject' must be of class 'integrSiteObject'")
  stopifnot(length(typeAnalysis) == 1L & typeAnalysis %in% c("clfi", "cdif", "govt", "gdit"))
  CheckValidArguments(integrSiteObject, typeAnalysis, argumentsAnalysis)
  
  iSO = integrSiteObject
  tA = typeAnalysis
  aA = argumentsAnalysis
  
  cf = tA %in% c("clfi", "cdif")
  df = tA %in% c("gdit", "cdif")
  
  iSO = setInfoSelectedCovariates(iSO, selectedCovariates = aA$selectedCovariates, selCovAreFactors = aA$selCovAreFactors)
  iSO = setGroups(iSO, group1 = as.integer(aA$idFilesGroup1), group2 = as.integer(aA$idFilesGroup2))
  iSO = loadIntegrations(iSO, aA$excludeIntegrations)
  iSO = loadLocationsGenes(iSO, aA$addPromoter)
  #if(aA$robust) iSO = excludeLargestClone3(iSO, aA$selectedCovariates)
  iSO = mergeData(iSO)
  #if(aA$robust) iSO = excludeLargestClone2(iSO, aA$excludeInBothGroups)
  
  #iSO = setFactorsCombinedData(iSO, variables = c(c("chr", "strand", "group"), aA$selectedCovariates[aA$selCovAreFactors]))
  #iSO = setFactorsCombinedData(iSO, variables = c(c("chr", "strand", "group"), aA$selectedCovariates[which(aA$selCovAreFactors)]))
  #iSO = setFactorsCombinedData(iSO, variables = c(c("chr", "strand", "group"), aA$selectedCovariates[which(aA$selectedCovariates %in% aA$selCovAreFactors)]))
  iSO = setFactorsCombinedData(iSO)
  iSO = computeSplitIntegrationsGenes(iSO, aA$nminSplit, removeDuplicates = tA %in% c("clfi", "cdif"))
  
  if(cf) {
    iSO = computeResponse(iSO)
    #iSO = computeDesignLogistic(iSO, covariates = c(aA$testedCovariate, "group"[tA == "cdif"]))
    iSO = computeDesignLogistic(iSO, testedCovariate = aA$testedCovariate, otherCovariates = aA$selectedCovariates, differential = df)
    if(tA == "clfi") iSO = fitLogisticRegressionUnderNull(iSO)
    
    ts = computeTestStatistics(iSO, aA$testedCovariate, aA$testType, aA$numberCores, aA$sizeChunks)
    if(aA$testType == "score") ts = computePvalueScore(ts, aA$pvalueCorrectionMethod)
    if(aA$testType == "likRt") ts = computePvalueLikRt(ts, aA$pvalueCorrectionMethod, namesParameters = iSO$model$namesParameters)
  } else {
    iSO = computeLengthInvestigatedGenome(iSO)
    iSO = computeSumLengthOutcomeVariable(iSO)
    iSO = computeDesignLogLinear(iSO, dummyForDatasets = aA$dummyForDatasets, otherCovariates = aA$selectedCovariates, typeRegressionModel = aA$typeRegressionModel)
    
    ts = computeTestStatisticLogLinear(iSO, aA$numberCores)
    ts = computePvalueLgLin(ts, aA$pvalueCorrectionMethod)
  }
  
  ts = processTests(iSO, ts, isTarget = !cf, overOnly = !cf & !df)
  return(list(result = ts, table = iSO$data$combined, integrSiteObject = iSO))
}

iSiteRobustAnalysis = function(typeAnalysis, integrSiteObject, resultAnalysis, argumentsAnalysis, thresholdAdjPvalue = 0.05) {
  iSO = integrSiteObject
  re = resultAnalysis
  aA = argumentsAnalysis
  
  #iSO = an$integrSiteObject
  #re = an$result
  #si = resultAnalysis$adjpvalue <= thresholdAdjPvalue
  if(typeAnalysis %in% c("clfi", "cdif")) {
    ts = computeTestStatistics(iSO, aA$testedCovariate, testType = "robst", aA$numberCores, aA$sizeChunks)
    ts = computePvalueRobst(ts, aA$pvalueCorrectionMethod)
  } 
  data.table :: setkey(re, testIndex)
  data.table :: setkey(ts, testIndex)
  re = merge(re, ts, by = "testIndex", all.x = TRUE)
  re = re[order(re$robustTestStat,decreasing = T),]
  return(re)
}

iSiteRecursAnalysis = function(typeAnalysis, integrSiteObject, resultAnalysis, argumentsAnalysis, thresholdAdjPvalue = 0.05, maxIterations = 10L) {
  if(class(integrSiteObject) != "integrSiteObject") stop("argument 'integrSiteObject' must be of class 'integrSiteObject'")
  stopifnot(length(typeAnalysis) == 1L & typeAnalysis %in% c("clfi", "cdif", "govt", "gdit"))
  CheckValidArguments(integrSiteObject, typeAnalysis, argumentsAnalysis)
  
  apv = resultAnalysis[1, "adjpvalue"]
  if(apv > thresholdAdjPvalue) {
    warning("first adjusted p-value greater than threshold, argument 'resultAnalysis' is returned")
    return(resultAnalysis)
  } 
  else {
    partialResult = resultAnalysis[1,]
    argumentsAnalysis = excludeSignificantIntegrations(argumentsAnalysis, resultAnalysis)
    i = 1L
    end = FALSE
    while(!end & i <= maxIterations) {
      print(paste("analysis", i, "/", maxIterations))
      newResult = iSiteAnalysis(typeAnalysis, integrSiteObject, argumentsAnalysis)$result
      apv = newResult[1, "adjpvalue"]
      partialResult = rbind(partialResult, newResult[1, ])
      if(apv > thresholdAdjPvalue) end = TRUE
      else argumentsAnalysis = excludeSignificantIntegrations(argumentsAnalysis, newResult)
      i = i+1L
    }
    if(!end) warning("maximum number of iterations reached")
    partialResult = rbind(partialResult, newResult[-1,])
    #partialResult = processTests(integrSiteObject, partialResult, logLinear = typeAnalysis %in% c("govt", "gdit"))
    return(partialResult)
  }
}

excludeSignificantIntegrations = function(argumentsAnalysis, result) {
  l = length(argumentsAnalysis$excludeIntegrations)
  argumentsAnalysis$excludeIntegrations[[l + 1L]] = c(result$chr[1], result$start[1], result$end[1])
  return(argumentsAnalysis)
}

CheckValidArguments = function(iSO, tA, aA) {
  idF1 = aA$idFilesGroup1
  if(!is.numeric(idF1)) stop("argument 'idFilesGroup1' must be of class 'numeric'")
  if(length(idF1) > 0L) {
    if(length(idF1) != length(unique(idF1))) stop("argument 'idFilesGroup1' must not contain duplicates")
    if(min(idF1) < 1L | max(idF1) > length(iSO$design$id)) stop(paste("argument 'idFilesGroup1' must contain values from 1 to", length(integrSiteObject$design$id), "(number of included tables)"))
  }
  if(tA %in% c("cdif", "gdit")) {
    idF2 = aA$idFilesGroup2
    if(!is.numeric(idF2)) stop("argument and 'idFilesGroup2' must be of class 'numeric'")
    if(length(idF2) != length(unique(idF2))) stop("argument 'idFilesGroup2' must not contain duplicates")
    if(min(idF2) < 1L | max(idF2) > length(iSO$design$id)) stop(paste("argument 'idFilesGroup2' must contain values from 1 to", length(integrSiteObject$design$id), "(number of included tables)"))
    if(length(intersect(idF1, idF2)) != 0L) stop("arguments 'idFilesGroup1' and 'idFilesGroup2' must be disjoint'")
  }
  addP = aA$addPromoter
  if(!is.numeric(addP)) stop("argument 'addPromoter' must be of class 'numeric' or 'integer'")
  if(length(addP) != 1L) stop("argument 'addPromoter' must be of length 1")
  if(as.integer(addP) != addP) stop("argument 'addPromoter' must be an integer")
  if(addP < 0L) stop("argument 'addPromoter' must be greater or equal than zero")
  selC = aA$selectedCovariates
  if(!is.character(selC)) stop("argument 'selectedCovariates' must be of class 'character'")
  if(length(selC) > 0L) {
    if(! all(selC %in% colnames(iSO$design$covariatesIntegrations)) ) stop("some elements of argument 'selectedCovariates' are not included in 'covariatesIntegration'")
    if(any(duplicated(selC))) stop("argument 'selectedCovariates' must not contain duplicates")
  } 
  sCAF = aA$selCovAreFactors
  if(!is.logical(sCAF)) stop("argument 'selCovAreFactors' must be of class 'logical'")
  if(length(selC) != length(sCAF)) stop("arguments 'selectedCovariates' and 'selCovAreFactors' must have the same length")
  if(tA %in% c("clfi", "cdif")) {
    tesC = aA$testedCovariate
    if(length(tesC) != 1L) stop("argument 'testedCovariate' must be of length 1")
    if(!(tesC %in% selC)) stop("argument 'testedCovariate' is not included in 'selectedCovariates'")
    tesT = aA$testType
    if(tA %in% "clfi")  if(!(tesT %in% c("likRt", "score"))) stop("argument 'testType' can be 'likRt' or 'score'")
    if(tA %in% "cdif") if(tesT != "likRt") stop("argument 'testType' must be 'likRt' for differential clone fitness analysis")
  }
  if(tA %in% c("govt", "gdit")) {
    #outC = aA$outcomeVariable
    #if(!(outC %in% c("count", selC))) stop("error")
    tyRM = aA$typeRegressionModel
    if(!(tyRM %in% c("logistic", "poisson"))) stop("argument 'typeRegressionModel' can be 'logistic' or 'poisson'")
  }
  return(invisible(NULL))
} 
