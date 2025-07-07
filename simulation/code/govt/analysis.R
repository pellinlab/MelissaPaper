# govt

integrSiteObject = createIntegrSiteObject(folSamples)
integrSiteObject = setOrganism(integrSiteObject, "simulation", "../../../../../data/gene.bed")

nameFile = paste0(c("", "0", "00", "000")[1 + numberDigits - nd], i)
nfs = paste0("s", nameFile, "s", includedSamples, ".bed")

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nfs), 
                             covariatesIntegrations = data.frame())

result = iSiteGeneOverTarget(integrSiteObject, 
                             idFiles = includedSamples, 
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

