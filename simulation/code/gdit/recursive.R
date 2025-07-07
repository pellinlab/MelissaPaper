# gdit

argA = setArgumentsAnalysis("gdit", list(idFilesGroup1 = 1:length(includedSamples), 
                                         idFilesGroup2 = length(includedSamples) + 1:length(includedSamples), 
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
recres = iSiteRecursAnalysis("gdit", integrSiteObject, result$result, argA)

integrSiteObject = setWorkingDirectory(integrSiteObject, folRecRess)
nfres = paste0("r", nameFile, ".csv")
class(recres) = class(result$result)
StoreResults(integrSiteObject, list(result = recres), nameFileResult = nfres, storeTable = FALSE)



