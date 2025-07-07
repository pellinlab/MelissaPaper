#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(nameBMSCdata, nameAMSCdata))

result = iSiteGeneDiffTarget(integrSiteObject, idFilesGroup1 = 1:4, idFilesGroup2 = 5:8, 
                             addPromoter = 1000L,
                             dummyForDatasets = T,
                             nminSplit = 0L, 
                             numberCores = 8L)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_di_bmsc_vs_amsc.csv", storeTable = F)
