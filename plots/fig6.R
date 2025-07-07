#mainFolder = paste0(baseFolder, "/MELISSApaper")

source(paste0(mainFolder,"/package/heatmapGenes.R"))
source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))
source(paste0(mainFolder,"/package/epigenetic.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))

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
genTable = data.table :: as.data.table(read.table(paste0(mainFolder, "/analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed"), header = T, sep = "\t"))
source(paste0(mainFolder,"/package/heatChromosomes.R"))

genTable$positive = "fal"
genTable$positive[1] = "tru"
genTable$positive = factor(genTable$positive)

resultFolder = paste0(mainFolder,"/analyses/results")
whichResult = c("cf_di_bmsc_vs_amsc")[1]
resultFile = paste0(resultFolder, "/", whichResult, ".csv")
p3 = MiamiDifferentialPlot(fileResult = resultFile, filePlot = "", storePlot = F,
                           #cancerGeneNameList = cancerGeneNameList,
                           formatPlot = "svg", cloneFitness = T, robust = F, bpsDomain = T, emptyRobust = F,
                           thresholdPvalue = 0.10, scaleByPvalue = F,
                           #colors1 = c("#0054ff", "#6c9cff"), colors2 = c("#008b0b", "#72ba78"))
                           colors1 = c("#0054ff", "#0054ff"), colors2 = c("#72ba78", "#72ba78"))
p3 = p3 + 
  ggplot2 :: geom_vline(data = genTable[-1], ggplot2 :: aes(xintercept = center), color = "#72ba78", alpha = .3) + 
  ggplot2 :: geom_vline(data = genTable[ 1], ggplot2 :: aes(xintercept = center), color = "#0054ff", alpha = .3) + 
  #scale_color_manual(values = c("tru" = "blue", "fal" = "green")) +
  ggplot2 :: xlab(NULL)

breaks_1 = ggplot2 :: ggplot_build(p2)$layout$panel_scales_x[[1]]$breaks
labels_1 = ggplot2 :: ggplot_build(p2)$layout$panel_scales_x[[1]]$labels
x_limits_1 <- ggplot2 :: ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range
x_limits_2 <- ggplot2 :: ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range
breaks_2 = ggplot2 :: ggplot_build(p3)$layout$panel_scales_x[[1]]$breaks
x_range_max <- range(c(x_limits_1, x_limits_2))
p2 = p2 + ggplot2 :: scale_x_continuous(breaks = breaks_1, labels = labels_1, limits = x_range_max)
p3 = p3 + ggplot2 :: scale_x_continuous(breaks = breaks_2, labels = NULL, limits = x_range_max)

prop_plots = 0.65
p = p3 + p2 + patchwork :: plot_layout(ncol = 1, height = c(prop_plots, 1 - prop_plots), axes = "collect")
p

ggplot2 :: ggsave(paste0(mainFolder, "/figures/f6/miam_cf_di_bmsc_vs_amsc_comb.pdf"), plot = p, dpi = 300, height = 6, width = 12)

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
                       formatPlot = ".pdf",
                       timeTransparency = T,
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
chainFile = paste0(epigData, "hg19ToHg38.over.chain")
bldt = c("Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
         "Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
         "Mobilized_CD34_Primary_Cells")

epigeneticPlot(integrSiteObject, chainFile = chainFile, 
               inputFiles = paste0(epigData, c(paste0("GSM670028_BI.", bldt[1], ".H3K9me3.60.bed.gz"),
                                               paste0("GSM621398_BI.", bldt[2], ".H3K9me3.1.bed.gz"),
                                               paste0("GSM773049_BI.", bldt[3], ".H3K9me3.RO_01562.bed.gz"))),
               lengthWindow = 100000L, filePlot = paste0(mainFolder, "/figures/f6/epig_", "H3K9me3", ".png"))
epigeneticPlot(integrSiteObject, chainFile = chainFile,
               inputFiles = paste0(epigData, c(paste0("GSM669957_BI.", bldt[1], ".H3K27me3.60.bed.gz"),
                                               paste0("GSM621420_BI.", bldt[2], ".H3K27me3.1.bed.gz"),
                                               paste0("GSM773047_BI.", bldt[3], ".H3K27me3.RO_01562.bed.gz"))),
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
