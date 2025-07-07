fileResult = paste0(mainFolder, "/analyses/results/gt_di_bmsc_vs_amsc_chrom.csv")
result = data.table :: as.data.table(read.table(fileResult, header = T, sep = "\t"))
chrTable = data.table :: as.data.table(read.table(paste0(mainFolder, "/analyses/data/othr/hg38_chromosomeAnnotations.bed"), header = T, sep = "\t"))
genTable = data.table :: as.data.table(read.table(paste0(mainFolder, "/analyses/data/othr/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed"), header = T, sep = "\t"))

convert_geneName_to_value = function(name) {
  if(grepl("^c[0-9XYM]+$", name)) return(1L) # grepl("^c[0-9XYM]+$", "cM")
  if(grepl("^c[0-9XYM]+[pq]$", name)) return(2L) # grepl("^c[0-9XYM]+[pq]$", "c12q")
  if(grepl("^c[0-9XYM]+[pq][1-9]$", name)) return(3L) # grepl("^c[0-9XYM]+[pq][1-9]$", "c12q1")
  if(grepl("^c[0-9XYM]+[pq][1-9][1-9]$", name)) return(4L) # grepl("^c[0-9XYM]+[pq][1-9][1-9]$", "cXq12")
  if(grepl("^c[0-9XYM]+[pq][1-9][1-9].[1-9]$", name)) return(5L) # grepl("^^c[0-9XYM]+[pq][1-9][1-9].[1-9]$", "c12q23.1")
  if(grepl("^c[0-9XYM]+[pq][1-9][1-9].[1-9][1-9]$", name)) return(6L) # grepl("^c[0-9XYM]+[pq][1-9][1-9].[1-9][1-9]$", "c12q12.32")
  return(7L)
}

chrTable$level = purrr :: map_int(chrTable$geneName, convert_geneName_to_value)
chrTable = chrTable[level < 7]

result$level = purrr :: map_int(result$geneName, convert_geneName_to_value)
result = result[level < 7]

res = result[, .(geneName, testStat, signedTestStat, pvalue, adjpvalue)]
mer = merge(chrTable, res, by = "geneName", all.x = TRUE)
mer[is.na(testStat), testStat := 0]
mer[is.na(signedTestStat), signedTestStat := 0]
mer[is.na(pvalue), pvalue := 1]
mer[is.na(adjpvalue), adjpvalue := 1]
mer = mer[chr_int != 25]

data.table :: setorder(mer, -signedTestStat, chr_int, start)
#threshold05 = min(mer$testStat[which(mer$adjpvalue <= 0.05)])
threshold05 = min(mer$testStat[which(mer$adjpvalue <= 0.1)])

mer$start_bp = as.numeric(mer$start)
mer$end_bp = as.numeric(mer$end)

cumchr = vector("numeric", 25)
cum = cumchr[2] = max(mer[chr_int == 1L]$end)
for(i in 2:24) {
  mer[chr_int == i]$start_bp = mer[chr_int == i]$start_bp + cum
  mer[chr_int == i]$end_bp = mer[chr_int == i]$end_bp + cum
  cum = max(mer[chr_int == i]$end_bp)
  cumchr[i+1] = cum
}
mer$mid_bp = trunc((mer$start_bp + mer$end_bp) / 2)

data.table :: setorder(mer, level, start_bp)

mer$sig05 = mer$testStat >= min(mer$testStat[which(mer$adjpvalue <= 0.05)])
mer$posst = mer$signedTestStat >= 0

genTable = genTable[geneName %in% c("SETD2", "DMD")]
genTable = genTable[c(1, 10)]
genTable$pos = genTable$start + cumchr[genTable$chr_int]

maxOverlapGeneNames = 20
p = ggplot2 :: ggplot(data = mer) +
  ggplot2 :: geom_hline(yintercept = c(-1,1) * threshold05, color = "#f40000",linetype = "dashed") +
  ggplot2 :: geom_vline(xintercept = genTable$pos, color = "#f40000", alpha = .5, linewidth = 0.4) +
  ggplot2 :: geom_rect(ggplot2 :: aes(xmin = start_bp, xmax = end_bp, ymin = 0, ymax = signedTestStat, alpha = level/6, fill = posst)) +
  ggplot2 :: scale_fill_manual(values = c("#008b0b", "#0054ff")) +
  ggrepel :: geom_text_repel(data=subset(mer, sig05 == T), 
                             ggplot2 :: aes(x = mid_bp, y = signedTestStat, label = geneName, color = posst, alpha = sqrt(level)/sqrt(6), size = I(abs(signedTestStat) / max(abs(signedTestStat)) * 2 + 3)), 
                             max.overlaps = getOption("ggrepel.max.overlaps", default = maxOverlapGeneNames)) + 
  ggplot2 :: scale_color_manual(values = c("#008b0b", "#0054ff")) +
  ggplot2 :: scale_x_continuous(breaks = sort(unique(c(mer[level == 1L]$start_bp, mer[level == 1L]$end_bp))), labels = NULL,
                                minor_breaks = sort(unique(c(mer[level == 2L]$start_bp, mer[level == 2L]$end_bp)))) +
  ggplot2 :: ylim(c(-1,1) * max(abs(mer$signedTestStat))) + 
  ggplot2 :: xlab("Chromosomes") + ggplot2 :: ylab("Band targeting score") +
  ggplot2 :: theme(panel.background = ggplot2::element_rect(fill = "white"),
                   legend.position = "none",
                   panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor.x = ggplot2::element_line(colour = "grey90"),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 14))
#p
