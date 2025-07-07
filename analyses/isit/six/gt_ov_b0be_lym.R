#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

dataFolder = "analyses/data/pati/"
da = paste0(dataFolder, c("b0bE_11_B0BE_Sixetal_B_.bed" ,         
                          "b0bE_11_B0BE_Sixetal_T_.bed" ,
                          "b0bE_11_B0BE_Sixetal_NK_.bed", 
                          "b0bE_36_B0BE_Sixetal_B_.bed" ,
                          "b0bE_36_B0BE_Sixetal_T_.bed" ,
                          "b0bE_36_B0BE_Sixetal_NK_.bed",
                          "b0bE_48_B0BE_Sixetal_B_.bed" ,  
                          "b0bE_48_B0BE_Sixetal_T_.bed" ,
                          "b0bE_48_B0BE_Sixetal_NK_.bed"))
cv = data.frame(time = rep(c(11L, 36L, 48L) - 0L, each = 3),
                type = rep(c("B", "T", "NK"), times = 3))
#type = rep(c(1, 2, 3), times = 3))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = da, covariatesIntegrations = cv)

result = iSiteGeneOverTarget(integrSiteObject, 
                             idFiles = 1:length(da), 
                             addPromoter = 1000L,
                             selectedCovariates = c("time", "type"), 
                             selCovAreFactors = c(F, T), 
                             nminSplit = 1L,
                             numberCores = 18L,
                             robust = F)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_ov_b0be_lym.csv", storeTable = F)

