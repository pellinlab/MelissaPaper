#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
ch = data.frame(time = c(0L, 0L, 0L))

integrSiteObject = setDesign(integrSiteObject, filesIntegrations = dh, covariatesIntegrations = ch)

result = iSiteGeneOverTarget(integrSiteObject, 
                             idFiles = 1:length(dh), 
                             addPromoter = 1000L,
                             selectedCovariates = c("time"), 
                             selCovAreFactors = c(F), 
                             nminSplit = 1L,
                             numberCores = 18L,
                             robust = F)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_ov_hspc_all.csv", storeTable = F)
