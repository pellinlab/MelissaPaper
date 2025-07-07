#mainFolder = paste0(baseFolder, "/MELISSApaper")

source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/venn.R"))
source(paste0(mainFolder,"/package/heatmapGenes.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))

# B
VennPlot(filesResult = paste0(mainFolder, "/analyses/results/", c("gt_ov_bmsc.csv", "gt_ov_amsc.csv", "gt_ov_hspc.csv")),
         nameCellTypes = c("BM MSC", "Adip MSC", "HSPC"),
         filePlot = paste0(mainFolder, "/figures/f4/venn.pdf"))

# A 
MiamiOvertargetPlot(fileResult = paste0(mainFolder, "/analyses/results/gt_ov_hspc.csv"), 
                    storePlot = T, formatPlot = "pdf", bpsDomain = T,
                    filePlot   = paste0(mainFolder, "/figures/f4/miam_gt_ov_hspc"),
                    cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"),
                    colors = c("#090909", "#7c7c7c"), sizePlot = c(16, 6))

# C
dev.off()
HeatmapGenes(imputedFiles = paste0(mainFolder, "/analyses/results/", c("gt_ov_bmsc.csv", "gt_ov_amsc.csv", "gt_ov_hspc.csv")), 
             filePlot = paste0(mainFolder, "/figures/f4/heat_govt.png"),
             cancerGeneNameList = paste0(mainFolder, "/analyses/data/othr/oncogenesBushman"), 
             includeGenes = c("LMO2", "CCND2", "BCL2", "MECOM", "PRDM16", "HMGA2", "HOXB7", "HOXB4", "NRAS", "KRAS", "ABL1"),
             typeCells = c("BM MSC", "Ad MSC", "HSPC"),
             plotOnlyInList = TRUE,
             cloneFitness = FALSE,
             numberTopGenes = 10L)

# D
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")

nameBMSCdata = paste0("analyses/data/expr/", c("BM_MSC_P1.bed", "BM_MSC_P4.bed", "BM_MSC_P6.bed", "BM_MSC_P8.bed"))
nameAMSCdata = paste0("analyses/data/expr/", c("Ad_MSC_P1.bed", "Ad_MSC_P4.bed", "Ad_MSC_P6.bed", "Ad_MSC_P8.bed"))
nameHSPCdata = paste0("analyses/data/expr/", c("CD34_HSPC_r1.bed", "CD34_HSPC_r2.bed", "CD34_HSPC_r3.bed"))
cv = data.frame(time = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("P", c(1,1,1))),
                #time = c(rep(c(0,1,2,3), times = 2), rep(0,3)), 
                type = c(rep("BM MSC", 4), rep("Ad MSC", 4), rep("HSPC", 3)),
                idty = c(rep(paste0("P", c(1,4,6,8)), 2), paste0("r", c(1:3))))

integrSiteObject = setDesign(integrSiteObject, 
                             filesIntegrations = c(nameBMSCdata, nameAMSCdata, nameHSPCdata),
                             covariatesIntegrations = cv)

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig4.bed")
fileAnnotations = paste0(mainFolder, "/analyses/data/othr/gencode.v45.annotation.gtf")
cat("Import gencode annotation\n")
annotations <- rtracklayer::import(fileAnnotations)

cat("Gene focus plot\n")
FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f4/"),
                          integrSiteObject = integrSiteObject,
                          transcriptsDataBase = annotations,
                          one_transcript_plotted = TRUE,
                          proportion_transcripts = .3,
                          aes_y_char = "time", 
                          aes_shape_char = "strand", 
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "Passage",
                          covariatesInLabels = "idty",
                          formatPlot = ".pdf",
                          transparency = 1,
                          jitter = 0,
                          sizeInchesPlots = c(8, 2.5),
                          colors = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909"))

# E
transcript_table = data.table::data.table(read.table(paste0(mainFolder, "/analyses/data/soft/transcriptTable.bed"), sep = "\t", header = T))$geneName
resuIsFile = paste0(mainFolder, "/analyses/results/gt_ov_hspc.csv")
resuGeFile = paste0(mainFolder, "/analyses/data/soft/hsoc_mean_expr.csv") 
colPoints = c("#090909", "#7c7c7c")[2:1]
labaxs = "HSPC"
xpos = 7.5e04
ypos = c(-70, -85)
lxpos = 3

source(paste0(mainFolder,"/package/plotExpressionTarget.R"))
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f4/expr_hspc.pdf"), plot = pll, width = 8, height = 4, dpi = 300)

transcript_table = data.table::data.table(read.table(paste0(mainFolder, "/analyses/data/soft/transcriptTable.bed"), sep = "\t", header = T))$geneName
resuIsFile = paste0(mainFolder, "/analyses/results/gt_ov_bmsc.csv")
resuGeFile = paste0(mainFolder, "/analyses/data/soft/bmsc_mean_expr.csv") 
colPoints = c("#0054ff", "#6c9cff")[2:1]
labaxs = "BM MSC"
xpos = 1.5e04
ypos = c(-55, -65)
lxpos = 2

source(paste0(mainFolder,"/package/plotExpressionTarget.R"))
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f4/expr_bmsc.pdf"), plot = pll, width = 8, height = 4, dpi = 300)
