# clfi

argA = setArgumentsAnalysis("clfi", list(idFiles = 1:length(includedSamples), 
                                         testedCovariate = "timePointInMonths",
                                         addPromoter = addPromoter,
                                         selectedCovariates = "timePointInMonths", 
                                         selCovAreFactors = rep(F, 1), 
                                         dummyForDatasets = dummyForDatasets,
                                         testType = "likRt", 
                                         pvalueCorrectionMethod = pvalueCorrectionMethod, 
                                         nminSplit = nminSplit, 
                                         numberCores = nCores, 
                                         sizeChunks = sizeChunks,
                                         robust = F,
                                         excludeInBothGroups = F,
                                         excludeIntegrations = list()))

integrSiteObject = setWorkingDirectory(integrSiteObject, folSamples)
recres = iSiteRecursAnalysis("clfi", integrSiteObject, result$result, argA)

integrSiteObject = setWorkingDirectory(integrSiteObject, folRecRess)
nfres = paste0("r", nameFile, ".csv")
class(recres) = class(result$result)
StoreResults(integrSiteObject, list(result = recres), nameFileResult = nfres, storeTable = FALSE)



