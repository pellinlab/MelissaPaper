mainFolder = "~/Downloads/MELISSApaper"

source(paste0(mainFolder,"/package/descriptive.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
cv = data.frame(time = c(rep(c(0,1,2,3), times = 2), rep(0,3)), 
                type = c(rep("BM MSC", 4), rep("Ad MSC", 4), rep("HSPC", 3)),
                idty = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("r", c(1:3))))

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nameBMSCdata, nameAMSCdata, nameHSPCdata),
                             covariatesIntegrations = cv)

# A
CloneSizesPlot(integrSiteObject, filePlot = paste0(mainFolder, "/figures/f3/cloneSize.png"))

# B
DiversityIndexPlot(integrSiteObject, filePlot = paste0(mainFolder, "/figures/f3/diversityIndex.png"))

# C
AnnotationDistributionPlot(integrSiteObject, filePlot = paste0(mainFolder, "/figures/f3/annotationDistribution.png"))
