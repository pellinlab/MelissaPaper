# clfi

integrSiteObject = createIntegrSiteObject(folSamples)
integrSiteObject = setOrganism(integrSiteObject, "simulation", "../../../../../data/gene.bed")

nameFile = paste0(c("", "0", "00", "000")[1 + numberDigits - nd], i)
nfs = paste0("s", nameFile, "T", includedSamples, ".bed")

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nfs), 
                             covariatesIntegrations = data.frame(timePointInMonths = c(includedSamples - 1L)))

result = iSiteCloneFitness(integrSiteObject, 
                           idFiles = includedSamples, 
                           testedCovariate = "timePointInMonths",
                           addPromoter = addPromoter,
                           selectedCovariates = "timePointInMonths",
                           testType = "likRt", 
                           pvalueCorrectionMethod = pvalueCorrectionMethod,
                           nminSplit = nminSplit, 
                           numberCores = nCores,
                           sizeChunks = sizeChunks,
                           recursive = FALSE)

integrSiteObject = setWorkingDirectory(integrSiteObject, folResults)
nfres = paste0("r", nameFile, ".csv")
StoreResults(integrSiteObject, result, nameFileResult = nfres, storeTable = FALSE)

