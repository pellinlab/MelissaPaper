# cdif

argA = setArgumentsAnalysis("cdif", list(idFilesGroup1 = 1:length(includedSamples), 
                                         idFilesGroup2 = length(includedSamples) + 1:length(includedSamples), 
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
recres = iSiteRecursAnalysis("cdif", integrSiteObject, result$result, argA)

integrSiteObject = setWorkingDirectory(integrSiteObject, folRecRess)
nfres = paste0("r", nameFile, ".csv")
class(recres) = class(result$result)
StoreResults(integrSiteObject, list(result = recres), nameFileResult = nfres, storeTable = FALSE)



