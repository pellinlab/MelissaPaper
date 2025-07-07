# cdif

integrSiteObject = createIntegrSiteObject(folSamples)
integrSiteObject = setOrganism(integrSiteObject, "simulation", "../../../../../data/gene.bed")

nameFile = paste0(c("", "0", "00", "000")[1 + numberDigits - nd], i)
nfs1 = paste0("s", nameFile, "T", includedSamples, "G1.bed")
nfs2 = paste0("s", nameFile, "T", includedSamples, "G2.bed")

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nfs1, nfs2), 
                             covariatesIntegrations = data.frame(timePointInMonths = rep(includedSamples - 1L, times = 2)))

result = iSiteCloneDiffFitness(integrSiteObject, 
                               idFilesGroup1 = includedSamples, 
                               idFilesGroup2 = length(includedSamples) + includedSamples, 
                               testedCovariate = "timePointInMonths",
                               addPromoter = addPromoter,
                               selectedCovariates = "timePointInMonths",
                               pvalueCorrectionMethod = pvalueCorrectionMethod,
                               nminSplit = nminSplit, 
                               numberCores = nCores,
                               sizeChunks = sizeChunks,
                               recursive = FALSE)

integrSiteObject = setWorkingDirectory(integrSiteObject, folResults)
nfres = paste0("r", nameFile, ".csv")
StoreResults(integrSiteObject, result, nameFileResult = nfres, storeTable = FALSE)

