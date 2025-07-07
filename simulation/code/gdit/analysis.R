# gdit

integrSiteObject = createIntegrSiteObject(folSamples)
integrSiteObject = setOrganism(integrSiteObject, "simulation", "../../../../../data/gene.bed")

nameFile = paste0(c("", "0", "00", "000")[1 + numberDigits - nd], i)
nfs1 = paste0("s", nameFile, "s", includedSamples, "g1.bed")
nfs2 = paste0("s", nameFile, "s", includedSamples, "g2.bed")

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nfs1, nfs2), 
                             covariatesIntegrations = data.frame())

result = iSiteGeneDiffTarget(integrSiteObject, 
                             idFilesGroup1 = includedSamples, 
                             idFilesGroup2 = length(includedSamples) + includedSamples, 
                             #outcomeVariable = "count",
                             addPromoter = addPromoter,
                             selectedCovariates = character(),
                             dummyForDatasets = dummyForDatasets,
                             typeRegressionModel = "logistic", 
                             pvalueCorrectionMethod = pvalueCorrectionMethod,
                             nminSplit = nminSplit, 
                             numberCores = nCores,
                             recursive = FALSE)

integrSiteObject = setWorkingDirectory(integrSiteObject, folResults)
nfres = paste0("r", nameFile, ".csv")
StoreResults(integrSiteObject, result, nameFileResult = nfres, storeTable = FALSE)
