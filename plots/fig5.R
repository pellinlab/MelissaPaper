mainFolder = "~/Downloads/MELISSApaper"

source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/GoKeggReactome.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))

# A
MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
                      filePlot   = paste0(mainFolder, "/figures/f5/miam_gt_di_bmsc_vs_hspc"),
                      cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                      formatPlot = "pdf",
                      colors1 = c("#0054ff", "#6c9cff"),
                      colors2 = c("#090909", "#7c7c7c"))

# B
BushmanDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
                        filePlot   = paste0(mainFolder, "/figures/f5/mibu_gt_di_bmsc_vs_hspc.pdf"),
                        cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                        maxOverlapPlot = 10)

# C # randomness in how Normalized Enrichment Score is computed
GoKeggPlots(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
            folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
            annotateByFdrPvalue = T,
            namePlots = paste0("goke_bm_vs_hs_", c("keggScore","","","",""), ".png"),
            cloneFitness = F,
            sizeInchesPlots = list(ks = c(10, 4), kh = c(), re = c(), gs = c(), gh = c()))

# D # randomness in how Normalized Enrichment Score is computed
GoKeggPlots(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_amsc_vs_hspc.csv"),
            folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
            annotateByFdrPvalue = T,
            namePlots = paste0("goke_ad_vs_hs_", c("","","","goScore",""), ".png"),
            cloneFitness = F,
            sizeInchesPlots = list(ks = c(), kh = c(), re = c(), gs = c(15, 4), gh = c()))

# E
MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_amsc.csv"),
                      filePlot   = paste0(mainFolder, "/figures/f5/miam_gt_di_bmsc_vs_amsc"),
                      cancerGeneNameList = "/home/giacomo/Work/Boston/Melissa/Old/oncogenesBushman",
                      cloneFitness = FALSE,
                      thresholdPvalue = 0.1,
                      formatPlot = "pdf",
                      colors1 = c("#0054ff", "#6c9cff"),
                      colors2 = c("#008b0b", "#72ba78"))

# F
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

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig5.bed")

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
                          integrSiteObject = integrSiteObject,
                          aes_y_char = "time", 
                          aes_shape_char = "strand", 
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          #selectedFactorCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "Passage",
                          covariatesInLabels = "idty",
                          formatPlot = ".pdf",
                          transparency = 1,
                          jitter = 0,
                          sizeInchesPlots = c(9, 2),
                          showLegend = F,
                          colors = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909"))
