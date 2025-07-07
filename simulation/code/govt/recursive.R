# govt

argA = setArgumentsAnalysis("govt", list(idFiles = 1:length(includedSamples), 
                                         #outcomeVariable = "count",
                                         addPromoter = addPromoter,
                                         selectedCovariates = character(), 
                                         selCovAreFactors = rep(F, 0), 
                                         dummyForDatasets = dummyForDatasets,
                                         typeRegressionModel = "logistic", 
                                         pvalueCorrectionMethod = pvalueCorrectionMethod, 
                                         nminSplit = nminSplit, 
                                         numberCores = nCores, 
                                         robust = F,
                                         excludeInBothGroups = F,
                                         excludeIntegrations = list()))

integrSiteObject = setWorkingDirectory(integrSiteObject, folSamples)
recres = iSiteRecursAnalysis("govt", integrSiteObject, result$result, argA)

integrSiteObject = setWorkingDirectory(integrSiteObject, folRecRess)
nfres = paste0("r", nameFile, ".csv")
class(recres) = class(result$result)
StoreResults(integrSiteObject, list(result = recres), nameFileResult = nfres, storeTable = FALSE)



