#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
cv = data.frame(time = rep(c(0,1,2,3), times = 1))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = nameBMSCdata, covariatesIntegrations = cv)

result = iSiteCloneFitness(integrSiteObject, idFiles = 1:4, 
                           addPromoter = 1000L,
                           testedCovariate = "time", 
                           selectedCovariates = "time", 
                           dummyForDatasets = T,
                           nminSplit = 2L,
                           numberCores = 16)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/cf_ov_bmsc.csv", storeTable = F)
