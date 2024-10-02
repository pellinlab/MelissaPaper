mainFolder = "~/Downloads/MELISSApaper"

source(paste0(mainFolder,"/package/geneStatisticPlot.R"))
source(paste0(mainFolder,"/package/miamiBushman.R"))
source(paste0(mainFolder,"/package/focusGenes.R"))


library(data.table)
library(ggplot2)
library(purrr)

# A
dataFolder = paste0(mainFolder, "/analyses/data/pati/")
d_bt_ = paste0(dataFolder, c("b0bE_11_B0BE_Sixetal_B_.bed","b0bE_11_B0BE_Sixetal_T_.bed","b0bE_11_B0BE_Sixetal_NK_.bed","b0bE_11_B0BE_Sixetal_G_.bed","b0bE_11_B0BE_Sixetal_Mono_.bed",
                               "b0bE_36_B0BE_Sixetal_B_.bed","b0bE_36_B0BE_Sixetal_T_.bed","b0bE_36_B0BE_Sixetal_NK_.bed","b0bE_36_B0BE_Sixetal_G_.bed","b0bE_36_B0BE_Sixetal_Mono_.bed",
                               "b0bE_48_B0BE_Sixetal_B_.bed","b0bE_48_B0BE_Sixetal_T_.bed","b0bE_48_B0BE_Sixetal_NK_.bed","b0bE_48_B0BE_Sixetal_G_.bed","b0bE_48_B0BE_Sixetal_Mono_.bed"))
d_scd = paste0(dataFolder, c("bSbS_12_SCD_Sixetal_B_.bed","bSbS_12_SCD_Sixetal_T_.bed","bSbS_12_SCD_Sixetal_NK_.bed","bSbS_12_SCD_Sixetal_G_.bed","bSbS_12_SCD_Sixetal_Mono_.bed",
                               "bSbS_24_SCD_Sixetal_B_.bed","bSbS_24_SCD_Sixetal_T_.bed","bSbS_24_SCD_Sixetal_NK_.bed","bSbS_24_SCD_Sixetal_G_.bed","bSbS_24_SCD_Sixetal_Mono_.bed"))
dw1 = paste0(dataFolder, c("WAS2_22_WAS_Sixetal_B_.bed","WAS2_22_WAS_Sixetal_T_.bed","WAS2_22_WAS_Sixetal_NK_.bed","WAS2_22_WAS_Sixetal_G_.bed","WAS2_22_WAS_Sixetal_Mono_.bed" ,
                             "WAS2_48_WAS_Sixetal_B_.bed","WAS2_48_WAS_Sixetal_T_.bed","WAS2_48_WAS_Sixetal_NK_.bed","WAS2_48_WAS_Sixetal_G_.bed","WAS2_48_WAS_Sixetal_Mono_.bed",
                             "WAS2_78_WAS_Sixetal_B_.bed","WAS2_78_WAS_Sixetal_T_.bed","WAS2_78_WAS_Sixetal_NK_.bed","WAS2_78_WAS_Sixetal_G_.bed","WAS2_78_WAS_Sixetal_Mono_.bed"))
dw2 = paste0(dataFolder, c("WAS4_12_WAS_Sixetal_B_.bed","WAS4_12_WAS_Sixetal_T_.bed","WAS4_12_WAS_Sixetal_NK_.bed","WAS4_12_WAS_Sixetal_G_.bed","WAS4_12_WAS_Sixetal_Mono_.bed",
                             "WAS4_36_WAS_Sixetal_B_.bed","WAS4_36_WAS_Sixetal_T_.bed","WAS4_36_WAS_Sixetal_NK_.bed","WAS4_36_WAS_Sixetal_G_.bed","WAS4_36_WAS_Sixetal_Mono_.bed",
                             "WAS4_48_WAS_Sixetal_B_.bed","WAS4_48_WAS_Sixetal_T_.bed","WAS4_48_WAS_Sixetal_NK_.bed","WAS4_48_WAS_Sixetal_G_.bed","WAS4_48_WAS_Sixetal_Mono_.bed",
                             "WAS4_60_WAS_Sixetal_B_.bed","WAS4_60_WAS_Sixetal_T_.bed","WAS4_60_WAS_Sixetal_NK_.bed","WAS4_60_WAS_Sixetal_G_.bed","WAS4_60_WAS_Sixetal_Mono_.bed"))
dw3 = paste0(dataFolder, c("WAS5_13_WAS_Sixetal_B_.bed","WAS5_13_WAS_Sixetal_T_.bed","WAS5_13_WAS_Sixetal_NK_.bed","WAS5_13_WAS_Sixetal_G_.bed","WAS5_13_WAS_Sixetal_Mono_.bed",
                             "WAS5_36_WAS_Sixetal_B_.bed","WAS5_36_WAS_Sixetal_T_.bed","WAS5_36_WAS_Sixetal_NK_.bed","WAS5_36_WAS_Sixetal_G_.bed","WAS5_36_WAS_Sixetal_Mono_.bed",
                             "WAS5_43_WAS_Sixetal_B_.bed","WAS5_43_WAS_Sixetal_T_.bed","WAS5_43_WAS_Sixetal_NK_.bed","WAS5_43_WAS_Sixetal_G_.bed","WAS5_43_WAS_Sixetal_Mono_.bed",
                             "WAS5_55_WAS_Sixetal_B_.bed","WAS5_55_WAS_Sixetal_T_.bed","WAS5_55_WAS_Sixetal_NK_.bed","WAS5_55_WAS_Sixetal_G_.bed","WAS5_55_WAS_Sixetal_Mono_.bed"))
dw4 = paste0(dataFolder, c("WAS7_12_WAS_Sixetal_B_.bed","WAS7_12_WAS_Sixetal_T_.bed","","WAS7_12_WAS_Sixetal_G_.bed","WAS7_12_WAS_Sixetal_Mono_.bed",
                             "WAS7_30_WAS_Sixetal_B_.bed","WAS7_30_WAS_Sixetal_T_.bed","WAS7_30_WAS_Sixetal_NK_.bed","WAS7_30_WAS_Sixetal_G_.bed","WAS7_30_WAS_Sixetal_Mono_.bed",
                             "WAS7_48_WAS_Sixetal_B_.bed","WAS7_48_WAS_Sixetal_T_.bed","WAS7_48_WAS_Sixetal_NK_.bed","WAS7_48_WAS_Sixetal_G_.bed","WAS7_48_WAS_Sixetal_Mono_.bed"))

data_bt = data.table()
for(i in 1:length(d_bt_)) {
  da = as.data.table(read.table(paste0(d_bt_[i])))
  data_bt = rbind(data_bt, data.table(nint = nrow(da), scnt = sum(da$V5)))
}
data_bt$time = rep(c(11, 36, 48), each = 5)
data_bt$type = rep(c("B", "T", "NK", "G", "Mono"), times = 3)
data_bt$pati = "BT"
rm(da)
data_sc = data.table()
for(i in 1:length(d_scd)) {
  da = as.data.table(read.table(paste0(d_scd[i])))
  data_sc = rbind(data_sc, data.table(nint = nrow(da), scnt = sum(da$V5)))
}
data_sc$time = rep(c(12, 48), each = 5)
data_sc$type = rep(c("B", "T", "NK", "G", "Mono"), times = 2)
data_sc$pati = "SCD"
rm(da)
data_w1 = data.table()
for(i in 1:length(dw1)) {
  da = as.data.table(read.table(paste0(dw1[i])))
  data_w1 = rbind(data_w1, data.table(nint = nrow(da), scnt = sum(da$V5)))
}
data_w1$time = rep(c(22, 48, 78), each = 5)
data_w1$type = rep(c("B", "T", "NK", "G", "Mono"), times = 3)
data_w1$pati = "WAS1"
rm(da)
data_w2 = data.table()
for(i in 1:length(dw2)) {
  da = as.data.table(read.table(paste0(dw2[i])))
  data_w2 = rbind(data_w2, data.table(nint = nrow(da), scnt = sum(da$V5)))
}
data_w2$time = rep(c(12, 36, 48, 60), each = 5)
data_w2$type = rep(c("B", "T", "NK", "G", "Mono"), times = 4)
data_w2$pati = "WAS2"
rm(da)
data_w3 = data.table()
for(i in 1:length(dw3)) {
  da = as.data.table(read.table(paste0(dw3[i])))
  data_w3 = rbind(data_w3, data.table(nint = nrow(da), scnt = sum(da$V5)))
}
data_w3$time = rep(c(13, 36, 43, 55), each = 5)
data_w3$type = rep(c("B", "T", "NK", "G", "Mono"), times = 4)
data_w3$pati = "WAS3"
rm(da)
data_w4 = data.table()
for(i in (1:length(dw4))) {
  if(i != 3) {
    da = as.data.table(read.table(paste0(dw4[i])))
    data_w4 = rbind(data_w4, data.table(nint = nrow(da), scnt = sum(da$V5)))
  }
  else data_w4 = rbind(data_w4, data.table(nint = 0L, scnt = 0L))
}
data_w4$time = rep(c(12, 30, 48), each = 5)
data_w4$type = rep(c("B", "T", "NK", "G", "Mono"), times = 3)
data_w4$pati = "WAS4"
rm(da)
data_al = rbind(data_bt, data_sc, data_w1, data_w2, data_w3, data_w4)
#data_wa = rbind(data_w1, data_w2, data_w3, data_w4)

nobs = c(15, 10, 15, 20, 20, 15)
csumData = cumsum(c(0, head(nobs, -1)))
set_x = function(n) unlist(lapply(1:5, function(x) x + (0:(n-1))*5))
ordTable = rep(csumData, times = nobs) + c(set_x(3), set_x(2), set_x(3), set_x(4), set_x(4), set_x(3))
data_al = data_al[ordTable,]
shiftPatients = 0.5
shiftDisease  = 0.5
data_al$x = 1:95 + shiftPatients * rep(0:5, times = nobs) + shiftDisease * rep(0:2, times = c(15, 10, 70))
wicol = .9

p = ggplot(data = data_al)
p = p + geom_col(aes(x = x, y = scnt, fill = type), width = wicol)
p = p + geom_segment(aes(x = x - wicol/2, xend = x + wicol/2, y = nint, yend = nint))
p = p + scale_x_continuous(breaks = data_al$x, labels = map2_chr(c("","\n")[c(rep(1:2, times = 47),1)], data_al$time, ~ paste0(.x,.y)))
p = p + labs(x = "Time (m)", y = "nIS / Cum. clone size")
p = p + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = NA),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.minor.x = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
  panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
  axis.text.x = ggplot2::element_text(size = 11), 
  axis.ticks.x = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_text(size = 11))
p
ggsave(filename = paste0(mainFolder, "/figures/f7/patientData.pdf"), plot = p, width = 16.5, height = 2.5, dpi = 300)


# B
fileGeneList = paste0(mainFolder, "/analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed")
fileResults = paste0(mainFolder, "/analyses/results/", c("gt_ov_b0be_lym.csv","gt_ov_b0be_mye.csv","gt_ov_bsbs_lym.csv","gt_ov_bsbs_mye.csv","gt_ov_was__lym.csv","gt_ov_was__mye.csv","gt_ov_hspc_all.csv"))
idSamples = c("BT lym", "BT mye", "SCD lym", "SCD mye", "WAS lym", "WAS mye", "HSPC")
colourChromosomes = list(c("#530e53", "#8A458A"), c("#311557", "#6B4E90"), c("#803315", "#D4886A"), c("#721330", "#BE5F7C"), c("#806515", "#D4BA6A"), c("#805215", "#D4A76A"),c("#090909", "#7c7c7c")) 
p = genesStatisticPlot(fileGeneList, fileResults, idSamples, thresholdPvalue = 0.05, lineWidth = 0.4, transparency = 1, maxLineLength = Inf, inBasePairs = T)
ggsave(filename =paste0(mainFolder, "/figures/f7/gene_stat.pdf"), plot = p, width = 15, height = 4.5, dpi = 300)

# C / E
col_bsbs = list(all = c("#801515", "#D46A6A"), lym = c("#803315", "#D4886A"), mie = c("#721330", "#BE5F7C"))
col_was_ = list(all = c("#805c15", "#D4B16A"), lym = c("#806515", "#D4BA6A"), mie = c("#805215", "#D4A76A"))
col_hspc = list(all = c("#090909", "#7c7c7c"))
MiamiOvertargetPlot(fileResult = paste0(mainFolder, "/analyses/results/cf_ov_bsbs_all.csv"), 
                    filePlot   = paste0(mainFolder, "/figures/f7/miam_cf_ov_bsbs_all"),
                    cancerGeneNameList = "/home/giacomo/Work/Boston/Melissa/Old/oncogenesBushman",
                    formatPlot = "pdf",
                    cloneFitness = TRUE, typePlot = "2", colors = col_bsbs$all, robust = T, keepRobustStat = T, removeRobustDifferentSign = T)
MiamiDifferentialPlot(fileResult = paste0(mainFolder, "/analyses/results/gt_di_2_was__vs_hspc.csv"),
                      filePlot   = paste0(mainFolder, "/figures/f7/miam_gt_di_2_was__vs_hspc"),
                      cancerGeneNameList = "/home/giacomo/Work/Boston/Melissa/Old/oncogenesBushman",
                      formatPlot = "pdf",
                      cloneFitness = F, colors1 = col_was_$all, colors2 = col_hspc$all)


# D / F
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
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(d1, d2), covariatesIntegrations = rbind(c1, c2))

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig7D.bed")

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                          integrSiteObject = integrSiteObject,
                          aes_y_char = "time",
                          aes_shape_char = "strand",
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          #selectedFactorCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "time",
                          formatPlot = ".pdf",
                          transparency = .3,
                          jitter = 0.05,
                          showLegend = F)
FocusGenesTrackingPlot(fileGenes = fileGenes,
                       folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                       integrSiteObject = integrSiteObject,
                       covariatesInLabels = c("type", "time"),
                       formatPlot = ".pdf")

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
integrSiteObject = setDesign(integrSiteObject, filesIntegrations = c(d1, d2), covariatesIntegrations = rbind(c1, c2))

fileGenes = paste0(mainFolder, "/analyses/data/othr/focusFig7F.bed")

FocusGenesIntegrationPlot(fileGenes = fileGenes,
                          folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                          integrSiteObject = integrSiteObject,
                          aes_y_char = "time",
                          aes_shape_char = "strand",
                          aes_color_char = "type",
                          selectedCovariates = c("time", "type"),
                          selectedFactorCovariates = c(T, T),
                          ylab = "time",
                          formatPlot = ".pdf",
                          transparency = .3,
                          jitter = 0.05,
                          showLegend = F,
                          sizeInchesPlots = c(7, 2))
FocusGenesTrackingPlot(fileGenes = fileGenes,
                       folderStoredPlots = paste0(mainFolder, "/figures/f7/"),
                       integrSiteObject = integrSiteObject,
                       covariatesInLabels = c("type", "time"),
                       formatPlot = ".pdf")
