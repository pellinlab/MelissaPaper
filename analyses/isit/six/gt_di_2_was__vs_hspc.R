#mainFolder = "~/Downloads/MELISSApaper"
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

dataFolder = "analyses/data/pati/"
d2 = c("WAS2_22_WAS_Sixetal_B_.bed" ,         
       "WAS2_22_WAS_Sixetal_T_.bed" ,
       "WAS2_22_WAS_Sixetal_NK_.bed", 
       "WAS2_22_WAS_Sixetal_G_.bed" ,         
       "WAS2_22_WAS_Sixetal_Mono_.bed" ,
       "WAS2_48_WAS_Sixetal_B_.bed" ,         
       "WAS2_48_WAS_Sixetal_T_.bed" ,
       "WAS2_48_WAS_Sixetal_NK_.bed", 
       "WAS2_48_WAS_Sixetal_G_.bed" ,         
       "WAS2_48_WAS_Sixetal_Mono_.bed" ,
       "WAS2_78_WAS_Sixetal_B_.bed" ,         
       "WAS2_78_WAS_Sixetal_T_.bed" ,
       "WAS2_78_WAS_Sixetal_NK_.bed",
       "WAS2_78_WAS_Sixetal_G_.bed" ,         
       "WAS2_78_WAS_Sixetal_Mono_.bed")
d4 = c("WAS4_12_WAS_Sixetal_B_.bed" ,         
       "WAS4_12_WAS_Sixetal_T_.bed" ,
       "WAS4_12_WAS_Sixetal_NK_.bed", 
       "WAS4_12_WAS_Sixetal_G_.bed" ,         
       "WAS4_12_WAS_Sixetal_Mono_.bed" ,
       "WAS4_36_WAS_Sixetal_B_.bed" ,         
       "WAS4_36_WAS_Sixetal_T_.bed" ,
       "WAS4_36_WAS_Sixetal_NK_.bed", 
       "WAS4_36_WAS_Sixetal_G_.bed" ,         
       "WAS4_36_WAS_Sixetal_Mono_.bed" ,
       "WAS4_48_WAS_Sixetal_B_.bed" ,         
       "WAS4_48_WAS_Sixetal_T_.bed" ,
       "WAS4_48_WAS_Sixetal_NK_.bed", 
       "WAS4_48_WAS_Sixetal_G_.bed" ,         
       "WAS4_48_WAS_Sixetal_Mono_.bed" ,
       "WAS4_60_WAS_Sixetal_B_.bed" ,         
       "WAS4_60_WAS_Sixetal_T_.bed" ,
       "WAS4_60_WAS_Sixetal_NK_.bed",
       "WAS4_60_WAS_Sixetal_G_.bed" ,         
       "WAS4_60_WAS_Sixetal_Mono_.bed")
d5 = c("WAS5_13_WAS_Sixetal_B_.bed" ,         
       "WAS5_13_WAS_Sixetal_T_.bed" ,
       "WAS5_13_WAS_Sixetal_NK_.bed", 
       "WAS5_13_WAS_Sixetal_G_.bed" ,         
       "WAS5_13_WAS_Sixetal_Mono_.bed" ,
       "WAS5_36_WAS_Sixetal_B_.bed" ,         
       "WAS5_36_WAS_Sixetal_T_.bed" ,
       "WAS5_36_WAS_Sixetal_NK_.bed", 
       "WAS5_36_WAS_Sixetal_G_.bed" ,         
       "WAS5_36_WAS_Sixetal_Mono_.bed" ,
       "WAS5_43_WAS_Sixetal_B_.bed" ,         
       "WAS5_43_WAS_Sixetal_T_.bed" ,
       "WAS5_43_WAS_Sixetal_NK_.bed", 
       "WAS5_43_WAS_Sixetal_G_.bed" ,         
       "WAS5_43_WAS_Sixetal_Mono_.bed" ,
       "WAS5_55_WAS_Sixetal_B_.bed" ,         
       "WAS5_55_WAS_Sixetal_T_.bed" ,
       "WAS5_55_WAS_Sixetal_NK_.bed",
       "WAS5_55_WAS_Sixetal_G_.bed" ,         
       "WAS5_55_WAS_Sixetal_Mono_.bed")
d7 = c("WAS7_12_WAS_Sixetal_B_.bed" ,         
       "WAS7_12_WAS_Sixetal_T_.bed" ,
       # missing NK time 12
       "WAS7_12_WAS_Sixetal_G_.bed" ,         
       "WAS7_12_WAS_Sixetal_Mono_.bed" ,
       "WAS7_30_WAS_Sixetal_B_.bed" ,         
       "WAS7_30_WAS_Sixetal_T_.bed" ,
       "WAS7_30_WAS_Sixetal_NK_.bed", 
       "WAS7_30_WAS_Sixetal_G_.bed" ,         
       "WAS7_30_WAS_Sixetal_Mono_.bed" ,
       "WAS7_48_WAS_Sixetal_B_.bed" ,         
       "WAS7_48_WAS_Sixetal_T_.bed" ,
       "WAS7_48_WAS_Sixetal_NK_.bed",
       "WAS7_48_WAS_Sixetal_G_.bed" ,         
       "WAS7_48_WAS_Sixetal_Mono_.bed")

da = paste0(dataFolder, c(d2, d4, d5, d7))
cv = data.frame(time = c(rep(c(22L, 48L, 78L) - 0L, each = 5),
                         rep(c(12L, 36L, 48L, 60L) - 0L, each = 5),
                         rep(c(13L, 36L, 43L, 55L) - 0L, each = 5),
                         rep(c(12L, 30L, 48L) - 0L, each = 5)),
                type = rep(c("B", "T", "NK", "G", "Mono"), times = sum(c(3, 4, 4, 3))),
                pati = rep(c("W2", "W4", "W5", "W7"), times = 5 * c(3, 4, 4, 3)))
cv = cv[-58,]; row.names(cv) = 1:nrow(cv) # missing NK time 12 patient 7

dataFolder = "analyses/data/expr/"
dh = paste0(dataFolder, c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
ch = data.frame(time = c(0L, 0L, 1L), type = rep("HSPC", 3), pati = "HS")
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(da, dh), covariatesIntegrations = rbind(cv, ch))

result = iSiteGeneDiffTarget(integrSiteObject, 
                             idFilesGroup1 = 1:length(da), 
                             idFilesGroup2 = 1:length(dh) + length(da),
                             addPromoter = 1000L,
                             selectedCovariates = c("time", "type", "pati"), 
                             selCovAreFactors = c(F, T, T), 
                             nminSplit = 1L,
                             numberCores = 18,                           
                             robust = F,
                             excludeInBothGroups = FALSE)

StoreResults(integrSiteObject, result, nameFileResult = "/analyses/results/gt_di_2_was__vs_hspc.csv", storeTable = F)


