#mainFolder = paste0(baseFolder, "/MELISSApaper")

source(paste0(mainFolder,"/package/geneStatisticPlot.R"))
source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))

# library(data.table)
# library(ggplot2)
# library(purrr)

# A
fileGeneList = paste0(mainFolder, "/analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")
fileResults = paste0(mainFolder, "/analyses/results/", c("gt_ov_b0be_lym.csv","gt_ov_b0be_mye.csv","gt_ov_bsbs_lym.csv","gt_ov_bsbs_mye.csv","gt_ov_was__lym.csv","gt_ov_was__mye.csv","gt_ov_hspc_all.csv"))[7:1]
idSamples = c("\u03b2-thal\nLym", "\u03b2-thal\nMye", "SCD\nLym", "SCD\nMye", "WAS\nLym", "WAS\nMye", "HSPC")[7:1]
colourChromosomes = list(c("#530e53", "#8A458A"), c("#311557", "#6B4E90"), c("#803315", "#D4886A"), c("#721330", "#BE5F7C"), c("#806515", "#D4BA6A"), c("#805215", "#D4A76A") ,c("#090909", "#7c7c7c")) [7:1]
p = genesStatisticPlot(fileGeneList, fileResults, idSamples, thresholdPvalue = 0.05, lineWidth = 0.4, transparency = 1, maxLineLength = Inf, inBasePairs = T)
p = p + ggplot2 :: xlab("Genes over-targeted (adjusted p-value < 0.05) in at least one analysis")
ggplot2 :: ggsave(filename =paste0(mainFolder, "/figures/f7/gene_stat.pdf"), plot = p, width = 15, height = 4.5, dpi = 300)

# B / C
col_b0be = list(all = c("#3D1255", "#764B8E"), lym = c("#530e53", "#8A458A"), mie = c("#311557", "#6B4E90"))
col_bsbs = list(all = c("#801515", "#D46A6A"), lym = c("#803315", "#D4886A"), mie = c("#721330", "#BE5F7C"))
col_was_ = list(all = c("#805c15", "#D4B16A"), lym = c("#806515", "#D4BA6A"), mie = c("#805215", "#D4A76A"))
col_hspc = list(all = c("#090909", "#7c7c7c"))

p = MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/gt_di_2_was__vs_hspc.csv"), storePlot = F,
                          #filePlot   = paste0(mainFolder, "/figures/f7/miam_gt_di_2_was__vs_hspc"),
                          cancerGeneNameList = "/home/giacomo/Work/Boston/Melissa/Old/oncogenesBushman",
                          formatPlot = "pdf", bpsDomain = T,
                          cloneFitness = F, colors1 = col_was_$all, colors2 = col_hspc$all)
p = p + ggplot2 :: xlab(NULL)
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f7/miam_gt_di_2_was__vs_hspc", ".pdf"), plot = p, width = 15, height = 4, dpi = 300)

p = MiamiOvertargetPlot(fileResult = paste0(mainFolder, "/analyses/results/cf_ov_bsbs_all.csv"), storePlot = F,
                    #filePlot   = paste0(mainFolder, "/figures/f7/miam_cf_ov_bsbs_all"),
                    cancerGeneNameList = "/home/giacomo/Work/Boston/Melissa/Old/oncogenesBushman",
                    formatPlot = "pdf", log10Scale = F, keepAboveThreshold = 20, bpsDomain = T,
                    cloneFitness = TRUE, colors = col_bsbs$all)
p = p + ggplot2 :: xlab(NULL)
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f7/miam_cf_ov_bsbs_all", ".pdf"), plot = p, width = 15, height = 4, dpi = 300)



# D / E
fileAnnotations = paste0(mainFolder, "/analyses/data/othr/gencode.v45.annotation.gtf")
cat("Import gencode annotation\n")
annotations = rtracklayer::import(fileAnnotations)

source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

dataFolder = "analyses/data/pati/"
d1 = paste0(dataFolder, c("b0bE_11_B0BE_Sixetal_B_.bed","b0bE_11_B0BE_Sixetal_T_.bed","b0bE_11_B0BE_Sixetal_NK_.bed",
                          "b0bE_36_B0BE_Sixetal_B_.bed","b0bE_36_B0BE_Sixetal_T_.bed","b0bE_36_B0BE_Sixetal_NK_.bed",
                          "b0bE_48_B0BE_Sixetal_B_.bed","b0bE_48_B0BE_Sixetal_T_.bed","b0bE_48_B0BE_Sixetal_NK_.bed"))
d2 = paste0(dataFolder, c("b0bE_11_B0BE_Sixetal_G_.bed","b0bE_11_B0BE_Sixetal_Mono_.bed",
                          "b0bE_36_B0BE_Sixetal_G_.bed","b0bE_36_B0BE_Sixetal_Mono_.bed",
                          "b0bE_48_B0BE_Sixetal_G_.bed","b0bE_48_B0BE_Sixetal_Mono_.bed"))
c1 = data.frame(time = rep(c(11L, 36L, 48L) - 0L, each = 3), type = rep(c("B", "T", "NK"), times = 3))
c2 = data.frame(time = rep(c(11L, 36L, 48L) - 0L, each = 2), type = rep(c("G", "Mono"), times = 3))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(d1, d2)[c(1:3, 10:11, 4:6, 12:13, 7:9, 14:15)], covariatesIntegrations = rbind(c1, c2)[c(1:3, 10:11, 4:6, 12:13, 7:9, 14:15),])

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig7D.bed")

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                          integrSiteObject = integrSiteObject,
                          transcriptsDataBase = annotations,
                          one_transcript_plotted = TRUE,
                          proportion_transcripts = .3,
                          aes_y_char = "time",
                          aes_shape_char = "strand",
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "time",
                          formatPlot = ".pdf",
                          transparency = 1,
                          jitter = 0.2,
                          sizeInchesPlots = c(8, 2.5))
FocusGenesTrackingPlot(fileGenes = fileGenes,
                       folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                       integrSiteObject = integrSiteObject,
                       covariatesInLabels = c("type", "time"),
                       timeTransparency = T,
                       formatPlot = ".pdf",
                       sizeInchesPlots = c(7,3))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

dataFolder = "analyses/data/pati/"
d1 = paste0(dataFolder, c("bSbS_12_SCD_Sixetal_B_.bed","bSbS_12_SCD_Sixetal_T_.bed","bSbS_12_SCD_Sixetal_NK_.bed",
                          "bSbS_24_SCD_Sixetal_B_.bed","bSbS_24_SCD_Sixetal_T_.bed","bSbS_24_SCD_Sixetal_NK_.bed"))
d2 = paste0(dataFolder, c("bSbS_12_SCD_Sixetal_G_.bed","bSbS_12_SCD_Sixetal_Mono_.bed",
                          "bSbS_24_SCD_Sixetal_G_.bed","bSbS_24_SCD_Sixetal_Mono_.bed"))
c1 = data.frame(time = rep(c(12L, 24L) - 0L, each = 3), type = rep(c("B", "T", "NK"), times = 2))
c2 = data.frame(time = rep(c(12L, 24L) - 0L, each = 2), type = rep(c("G", "Mono"), times = 2))
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(d1, d2)[c(1:3, 7:8, 4:6, 9:10)], covariatesIntegrations = rbind(c1, c2)[c(1:3, 7:8, 4:6, 9:10),])

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig7F.bed")

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                          integrSiteObject = integrSiteObject,
                          transcriptsDataBase = annotations,
                          one_transcript_plotted = TRUE,
                          proportion_transcripts = .3,
                          aes_y_char = "time",
                          aes_shape_char = "strand",
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "time",
                          formatPlot = ".pdf",
                          transparency = 1,
                          jitter = 0.2,
                          sizeInchesPlots = c(8, 2.5))
FocusGenesTrackingPlot(fileGenes = fileGenes,
                       folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                       integrSiteObject = integrSiteObject,
                       covariatesInLabels = c("type", "time"),
                       timeTransparency = T,
                       formatPlot = ".pdf",
                       sizeInchesPlots = c(7, 3))
