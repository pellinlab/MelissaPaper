mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(nameAMSCdata, nameHSPCdata))

result = iSiteGeneDiffTarget(integrSiteObject, idFilesGroup1 = 1:4, idFilesGroup2 = 5:7, 
                             addPromoter = 1000L,
                             dummyForDatasets = T,
                             nminSplit = 1L, 
                             numberCores = 8L)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_di_amsc_vs_hspc.csv", storeTable = F)
