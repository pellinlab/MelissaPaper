mainFolder = "~/Downloads/MELISSApaper"

source(paste0(mainFolder,"/package/heatmapGenes.R"))
source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))
source(paste0(mainFolder,"/package/epigenetic.R"))

# A
HeatmapGenes(imputedFiles = paste0(mainFolder, "/analyses/results/", c("cf_ov_bmsc.csv", "cf_ov_amsc.csv")), 
             filePlot = paste0(mainFolder, "/figures/f6/heat_cfov.png"),
             cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
             includeGenes = c("LMO2", "CCND2", "BCL2", "MECOM", "PRDM16", "HMGA2", "HOXB7", "HOXB4", "NRAS", "KRAS", "ABL1"),
             typeCells = c("BM MSC", "Ad MSC"),
             plotOnlyInList = FALSE,
             numberTopGenes = 10L,
             cloneFitness = TRUE)

# B
MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "cf_di_bmsc_vs_amsc.csv"),
                      filePlot   = paste0(mainFolder, "/figures/f6/miam_cf_di_bmsc_vs_amsc"),
                      cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                      cloneFitness = TRUE, 
                      thresholdPvalue = 0.10,
                      formatPlot = "pdf",                      
                      colors1 = c("#0054ff", "#6c9cff"),
                      colors2 = c("#008b0b", "#72ba78"))

# C
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
cv = data.frame(time = c(rep(c(0,1,2,3), times = 2), rep(0,3)), 
                type = c(rep("BM MSC", 4), rep("Ad MSC", 4), rep("HSPC", 3)),
                idty = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("r", c(1:3))))

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nameBMSCdata, nameAMSCdata, nameHSPCdata),
                             covariatesIntegrations = cv)

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig6.bed")

FocusGenesTrackingPlot(fileGenes = fileGenes,
                       folderStoredPlots = paste0(mainFolder, "/figures/f6/"),
                       integrSiteObject = integrSiteObject,
                       covariatesInLabels = c("type", "idty"),
                       formatPlot = ".svg",
                       sizeInchesPlots = c(6, 4),
                       colors = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909"))

# D
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
cv = data.frame(time = c(rep(c(0,1,2,3), times = 2), rep(0,3)), 
                type = c(rep("BM MSC", 4), rep("Ad MSC", 4), rep("HSPC", 3)),
                idty = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("r", c(1:3))))

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nameBMSCdata, nameAMSCdata, nameHSPCdata),
                             covariatesIntegrations = cv)
epigData = paste0(mainFolder, "/analyses/data/epig/")
#epigData = paste0(mainFolder, "/../../../data/old/epigen/")
chainFile = paste0(epigData, "hg19ToHg38.over.chain")

epigeneticPlot(integrSiteObject, chainFile = chainFile, 
               inputFiles = paste0(epigData, c("GSM670028_BI.Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K9me3.60.bed.gz",
                                               "GSM621398_BI.Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K9me3.1.bed.gz",
                                               "GSM773049_BI.Mobilized_CD34_Primary_Cells.H3K9me3.RO_01562.bed.gz")),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K9me3", ".png"))
epigeneticPlot(integrSiteObject, chainFile = chainFile,
               inputFiles = paste0(epigData, c("GSM669957_BI.Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K27me3.60.bed.gz",
                                               "GSM621420_BI.Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K27me3.1.bed.gz",
                                               "GSM773047_BI.Mobilized_CD34_Primary_Cells.H3K27me3.RO_01562.bed.gz")),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K27me3", ".png"))
epigeneticPlot(integrSiteObject, chainFile = chainFile,
               inputFiles = paste0(epigData, c("GSM669920_BI.Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K4me1.59.bed.gz",
                                               "GSM772748_BI.Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K4me1.92.bed.gz",
                                               "GSM773043_BI.Mobilized_CD34_Primary_Cells.H3K4me1.RO_01562.bed.gz")),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K4me1", ".png"))
epigeneticPlot(integrSiteObject, chainFile = chainFile,
               inputFiles = paste0(epigData, c("GSM670019_BI.Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K4me3.60.bed.gz",
                                               "GSM772747_BI.Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K4me3.92.bed.gz",
                                               "GSM773041_BI.Mobilized_CD34_Primary_Cells.H3K4me3.RO_01562.bed.gz")),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K4me3", ".png"))
epigeneticPlot(integrSiteObject, chainFile = chainFile,
               inputFiles = paste0(epigData, c("GSM670037_BI.Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K36me3.58.bed.gz",
                                               "GSM772820_BI.Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells.H3K36me3.92.bed.gz",
                                               "GSM773042_BI.Mobilized_CD34_Primary_Cells.H3K36me3.RO_01562.bed.gz")),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K36me3", ".png"))
