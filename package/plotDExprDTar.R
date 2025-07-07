resuIsFile = paste0(mainFolder, "/analyses/results/gt_di_bmsc_vs_hspc.csv") # + bmsc, - hspc
resuGeFile = paste0(mainFolder, "/analyses/data/soft/bmsc_vs_hspc.csv") # - bmsc, + hspc
reverseGeneStat = T

reIs = data.table :: as.data.table(read.table(resuIsFile, header = T))
reGe = data.table :: as.data.table(read.table(resuGeFile, header = T, sep = ","))
if(reverseGeneStat) reGe$test_statistic = - reGe$test_statistic

reGe = reGe[!is.na(test_statistic)]
reGe$adj_pvalue[is.na(reGe$adj_pvalue)] = 1
names(reGe) = c("transcriptId", "expr_signed_teststat", "expr_adj_pvalue")

transcript_table = data.table :: data.table(
  transcriptId = stringr::str_split_fixed(transcript_table, "_", 2)[, 1],   
  geneName = stringr::str_split_fixed(transcript_table, "_", 2)[, 2]  
)

reGe = merge(reGe, transcript_table, by = "transcriptId", all = F)

mgd = merge(reIs[,.(geneName, signedTestStat, adjpvalue)], reGe, by = "geneName", all = TRUE)
mgd = mgd[!is.na(mgd$signedTestStat)]
data.table :: setorder(mgd, -signedTestStat)
mgd = mgd[!duplicated(mgd$geneName)]
mgd$expr_signed_teststat[is.na(mgd$expr_signed_teststat)] = 0
mgd$expr_adj_pvalue[is.na(mgd$expr_adj_pvalue)] = 1
mgd$transcriptId = NULL
mgd = mgd[!grepl("\\.", mgd$geneName)]

colPlot = c("#0054ff", "#6c9cff", "#090909", "#7c7c7c")

merged = mgd#[,.(signedTestStat, expr_signed_teststat, transcriptId, geneId, adjpvalue, expr_adj_pvalue, testStat)]
merged$annotate = merged$geneName #paste0(merged$transcriptId[integer()], "", merged$geneId)
merged$sigmin = (merged$adjpvalue < 0*1e-03) | (merged$expr_adj_pvalue < 0*1e-16) | (merged$adjpvalue < 5e-02 & merged$expr_adj_pvalue < 5e-02)
merged$cols = "hspc2"
merged$cols[merged$expr_signed_teststat < 0] = "hspc1"
merged$cols[merged$signedTestStat > 0] = "bmsc2"
merged$cols[merged$signedTestStat > 0 & merged$expr_signed_teststat > 0] = "bmsc1"

thr05int = min(abs(merged$signedTestStat)[merged$adjpvalue < 0.05])
thr05exp = min(abs(merged$expr_signed_teststat[merged$expr_adj_pvalue < 0.05]))

textSize = c(14)
textsizeGeneNames = 3
maxOverlapGeneNames = 20


plt = ggplot2 :: ggplot(data = merged, mapping = ggplot2 :: aes(x = expr_signed_teststat, y = signedTestStat)) +
  ggplot2 :: geom_vline(xintercept = 0) + ggplot2 :: geom_hline(yintercept = 0) +
  ggplot2 :: geom_point(ggplot2 :: aes(colour = cols)) +
  ggplot2 :: scale_color_manual(values = colPlot) +
  ggplot2 :: geom_hline(yintercept = c(1,-1) * thr05int, col = rep("#f40000",2), linetype = "dashed") +
  ggplot2 :: geom_vline(xintercept = c(1,-1) * thr05exp, col = rep("#f40000",2), linetype = "dashed") +
  ggplot2 :: geom_smooth(method = "lm", color = "red", fill = "pink") +
  ggrepel :: geom_text_repel(data=subset(merged, sigmin == T), ggplot2 :: aes(label = annotate), 
                             size = textsizeGeneNames, color = "grey20", max.overlaps = getOption("ggrepel.max.overlaps", default = maxOverlapGeneNames)) + 
  ggplot2 :: ylab(paste0("HSPC        BM MSC")) + 
  ggplot2 :: xlab(paste0("Gene differential expression score,")) + 
  ggplot2 :: theme(panel.background = ggplot2 :: element_rect(fill = "white"), 
                   legend.position = "none", 
                   panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor.x = ggplot2::element_line(colour = "grey90"),
                   axis.title.y = ggplot2::element_text(size = textSize),
                   axis.title.x = ggplot2::element_text(size = textSize),
                   axis.text.x = ggplot2::element_text(size = textSize), 
                   axis.text.y = ggplot2::element_text(size = textSize)) +
  ggplot2 :: ylim(c(-1,1) * max(abs(merged$signedTestStat)))
#plt

