#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = nameHSPCdata)

result = iSiteGeneOverTarget(integrSiteObject, idFiles = 1:3, 
                             addPromoter = 1000L,
                             dummyForDatasets = T,
                             nminSplit = 0L,
                             numberCores = 8L)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_ov_hspc.csv", nameFileTable = "/analyses/results/tab_hspc.csv", storeTable = T)

