mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

dataFolder = "analyses/data/pati/"
da = paste0(dataFolder, c("bSbS_12_SCD_Sixetal_B_.bed"   ,       
                          "bSbS_12_SCD_Sixetal_T_.bed"   ,
                          "bSbS_12_SCD_Sixetal_NK_.bed"  ,
                          "bSbS_12_SCD_Sixetal_G_.bed"   ,         
                          "bSbS_12_SCD_Sixetal_Mono_.bed",
                          "bSbS_24_SCD_Sixetal_B_.bed"   ,  
                          "bSbS_24_SCD_Sixetal_T_.bed"   ,
                          "bSbS_24_SCD_Sixetal_NK_.bed"  ,
                          "bSbS_24_SCD_Sixetal_G_.bed"   ,
                          "bSbS_24_SCD_Sixetal_Mono_.bed"))
cv = data.frame(time = rep(c(12L, 24L) - 0L, each = 5),
                type = rep(c("B", "T", "NK", "G", "Mono"), times = 2))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = da, covariatesIntegrations = cv)

result = iSiteCloneFitness(integrSiteObject, 
                           testedCovariate = "time", 
                           idFiles = 1:length(da), 
                           addPromoter = 1000L,
                           selectedCovariates = c("time", "type"), 
                           selCovAreFactors = c(F, T), 
                           nminSplit = 1L,
                           numberCores = 16L,
                           robust = T)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/cf_ov_bsbs_all.csv", storeTable = F)
