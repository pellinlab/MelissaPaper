#mainFolder = paste0(baseFolder, "/MELISSApaper")

source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/GoKeggReactome.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))

# A
MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
                      storePlot = T, formatPlot = "pdf",
                      filePlot   = paste0(mainFolder, "/figures/f5/miam_gt_di_bmsc_vs_hspc"),
                      cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                      colors1 = c("#0054ff", "#6c9cff"),
                      colors2 = c("#090909", "#7c7c7c"),
                      bpsDomain = F, sizePlot = c(9,5))

# B
transcript_table = data.table::data.table(read.table(paste0(mainFolder, "/analyses/data/soft/transcriptTable.bed"), sep = "\t", header = T))$geneName
source(paste0(mainFolder,"/package/plotDExprDTar.R"))
ggplot2 :: ggsave(paste0(mainFolder, "/figures/f5/bm_hs_diff.pdf"), plot = plt, height = 5, width = 3.5)

# C
BushmanDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
                        filePlot   = paste0(mainFolder, "/figures/f5/mibu_gt_di_bmsc_vs_hspc.pdf"),
                        cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                        maxOverlapPlot = 20,
                        sizePlot = c(3.8, 5))

# D
GoKeggPlots(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_hspc.csv"),
            folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
            annotateByFdrPvalue = T,
            namePlots = paste0("goke_bm_vs_hs_", c("keggScore","","","",""), ".png"),
            cloneFitness = F,
            sizeInchesPlots = list(ks = c(10, 4), kh = c(), re = c(), gs = c(), gh = c()))

# E
GoKeggPlots(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_amsc_vs_hspc.csv"),
            folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
            annotateByFdrPvalue = T,
            namePlots = paste0("goke_ad_vs_hs_", c("","","","goScore",""), ".png"),
            cloneFitness = F,
            sizeInchesPlots = list(ks = c(), kh = c(), re = c(), gs = c(15, 4), gh = c()))

# F
pm = MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/", "gt_di_bmsc_vs_amsc.csv"),
                           storePlot = F, formatPlot = "svg", bpsDomain = T,
                           cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                           cloneFitness = FALSE,
                           thresholdPvalue = 0.1,
                           colors1 = c("#0054ff", "#6c9cff"),
                           colors2 = c("#008b0b", "#72ba78"))
source(paste0(mainFolder,"/package/plot_gt_di_chrom.R"))
pc = p

breaks_1 = ggplot2 :: ggplot_build(pm)$layout$panel_scales_x[[1]]$breaks
breaks_2 = ggplot2 :: ggplot_build(pc)$layout$panel_scales_x[[1]]$minor_breaks
x_limits_1 <- ggplot2 :: ggplot_build(pm)$layout$panel_scales_x[[1]]$range$range
x_limits_2 <- ggplot2 :: ggplot_build(pc)$layout$panel_scales_x[[1]]$range$range
x_range_max <- range(c(x_limits_1, x_limits_2))
pm = pm + ggplot2 :: scale_x_continuous(breaks = breaks_1, labels = NULL, limits = x_range_max)
pc = pc + ggplot2 :: scale_x_continuous(breaks = breaks_1, labels = NULL, limits = x_range_max, minor_breaks = breaks_2)

pm = pm + ggplot2 :: xlab(NULL)
prop_plots = 0.5
p = pm + pc + patchwork::plot_layout(ncol = 1, height = c(prop_plots, 1 - prop_plots), axes = "collect")
#p
ggplot2 :: ggsave(paste0(mainFolder, "/figures/f5/miam_gt_di_bmsc_vs_amsc_comb.pdf"), plot = p, dpi = 300, height = 6, width = 9)

# G
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
cv = data.frame(time = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("P", c(1,1,1))), 
                type = c(rep("BM MSC", 4), rep("Ad MSC", 4), rep("HSPC", 3)),
                idty = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("r", c(1:3))))

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nameBMSCdata, nameAMSCdata, nameHSPCdata),
                             covariatesIntegrations = cv)

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig5.bed")
fileAnnotations = paste0(mainFolder, "/analyses/data/othr/gencode.v45.annotation.gtf")
cat("Import gencode annotation\n")
annotations = rtracklayer::import(fileAnnotations)

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f5/"),
                          integrSiteObject = integrSiteObject,
                          transcriptsDataBase = annotations,
                          one_transcript_plotted = TRUE,
                          proportion_transcripts = .3,
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
                          sizeInchesPlots = c(8, 2.5),
                          colors = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909"))
