reIs = data.table::as.data.table(read.table(resuIsFile, header = T))
reGe = data.table::as.data.table(read.table(resuGeFile, header = T, sep = ","))

reGe$test_statistic[is.na(reGe$test_statistic)] = 0
names(reGe) = c("transcriptId", "expr_signed_teststat")#, "expr_adj_pvalue")

transcript_table = data.table::data.table(
  transcriptId = stringr::str_split_fixed(transcript_table, "_", 2)[, 1],   
  geneName = stringr::str_split_fixed(transcript_table, "_", 2)[, 2]  
)

reGe = merge(reGe, transcript_table, by = "transcriptId", all = F)

mgd = merge(reIs[,.(geneName, signedTestStat, adjpvalue)], reGe, by = "geneName", all = TRUE)
mgd = mgd[!is.na(mgd$signedTestStat)]
data.table :: setorder(mgd, -signedTestStat)
mgd = mgd[!duplicated(mgd$geneName)]
mgd$expr_signed_teststat[is.na(mgd$expr_signed_teststat)] = 0
mgd$transcriptId = NULL

textsizeGeneNames = 3
maxOverlapGeneNames = 20


merged = mgd #mgd[,.(signedTestStat, expr_signed_teststat, transcriptId, geneId, adjpvalue, testStat)]
merged$annotate = merged$geneName #paste0(merged$transcriptId[integer()], "", merged$geneId)
data.table::setorder(merged, -signedTestStat)
merged$sigmin = F
merged$sigmin[1:20] = T
data.table::setorder(merged, -expr_signed_teststat)
merged$sigmin[1:20] = T
data.table::setorder(merged, -signedTestStat)
merged$sig05 = merged$adjpvalue < 5e-02 & merged$signedTestStat > 0
merged = merged[!grepl("\\.", merged$geneName)]

#merged$cols = c("bmsc", "hspc")[1L + (merged$signedTestStat < 0)]

thr05int = min(abs(merged$signedTestStat)[merged$adjpvalue < 0.05])
#thr05exp = min(abs(merged$expr_signed_teststat[merged$expr_adj_pvalue < 0.05]))

plt = ggplot2 :: ggplot(data = merged, mapping = ggplot2 :: aes(x = expr_signed_teststat, y = signedTestStat)) +
  ggplot2 :: geom_vline(xintercept = 0) + ggplot2 :: geom_hline(yintercept = 0) +
  ggplot2 :: geom_point(ggplot2 :: aes(color = sig05)) +
  ggplot2 :: scale_color_manual(values = colPoints) +
  ggplot2 :: geom_smooth(method = "lm", color = "red", fill = "pink") +
  ggplot2 :: geom_hline(yintercept = c(1) * thr05int, col = rep("#f40000",1), linetype = "dashed") +
  ggrepel :: geom_text_repel(data=subset(merged, sigmin == T), ggplot2 :: aes(x = expr_signed_teststat, y = signedTestStat, label = annotate), 
                  size = textsizeGeneNames, color = "grey20", max.overlaps = getOption("ggrepel.max.overlaps", default = maxOverlapGeneNames)) + 
  #geom_text(x = xpos, y = ypos[1], label = labeq[[1]], parse = TRUE, color = "red") +
  #geom_text(x = xpos, y = ypos[2], label = labeq[[2]], parse = TRUE, color = "red") +
  ggplot2 :: ylab(paste0("Gene target score, ", labaxs)) + ggplot2 :: xlab(paste0("Gene expression score, ", labaxs)) + 
  ggplot2 :: theme(panel.background = ggplot2 :: element_rect(fill = "white"), legend.position = "none")
#plt

pll = ggplot2 :: ggplot(data = merged, mapping = ggplot2 :: aes(x = expr_signed_teststat, y = signedTestStat)) +
  ggplot2 :: geom_vline(xintercept = 0) + ggplot2 :: geom_hline(yintercept = 0) +
  #geom_vline(xintercept = c(1,-1) * thr05exp, col = rep("#f40000",2), linetype = "dashed") +
  ggplot2 :: geom_point(ggplot2 :: aes(color = sig05)) +
  ggplot2 :: scale_color_manual(values = colPoints) +
  #geom_point(aes(alpha = sig05)) +
  #scale_alpha_manual(values = c(.3, 1)) +
  ggplot2 :: scale_x_log10() +
  #scale_y_log10() +
  #scale_color_manual(values = colPlot) +
  #stat_summary(fun.data=mean_cl_normal) +
  ggplot2 :: geom_smooth(method = "lm", color = "red", fill = "pink") +
  ggplot2 :: geom_hline(yintercept = c(1) * thr05int, col = rep("#f40000",1), linetype = "dashed") +
  ggrepel :: geom_text_repel(data=subset(merged, sigmin == T), ggplot2 :: aes(x = expr_signed_teststat, y = signedTestStat, label = annotate), 
                  size = textsizeGeneNames, color = "grey20", max.overlaps = getOption("ggrepel.max.overlaps", default = maxOverlapGeneNames)) + 
  #geom_text(x = lxpos, y = ypos[1], label = labeq[[1]], parse = TRUE, color = "red") +
  #geom_text(x = lxpos, y = ypos[2], label = labeq[[2]], parse = TRUE, color = "red") +
  #geom_text(x = 1, y = -70, label = , parse = TRUE) +
  ggplot2 :: ylab(paste0("Gene target score, ", labaxs)) + ggplot2 :: xlab(paste0("Gene expression log10 score, ", labaxs)) + 
  ggplot2 :: theme(panel.background = ggplot2 :: element_rect(fill = "white"), legend.position = "none",
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor.y = ggplot2::element_line(colour = "grey90"))
#pll
