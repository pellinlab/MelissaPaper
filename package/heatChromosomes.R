pltGenes =  c("COBLL1", "LIMD1", "CNBP", "PLD1", "ICE1", "RAD18", "FKBP5", "NSD1", "PDLIM7", "KMT2C", "PTEN", "DOCK1", "SP1", "RORA", "STK11", "SBNO2", "SPAG9", "NFIC", "MECP2", "DAZAP1")
genTable = genTable[geneName %in% pltGenes]
genTable$length = genTable$end - genTable$start
data.table :: setorder(genTable, -length)
genTable = genTable[!duplicated(genTable$geneName)]
genTable$posstat = "#72ba78"
genTable$posstat[genTable$geneName == "RORA"] = "#0054ff"
#genTable$center = as.integer((genTable$end + genTable$start) / 2)
genTable$center = genTable$start

differential = T

name_analysis = "cf_di_bm_vs_am_largewin"
result = data.table :: as.data.table(read.table(paste0(mainFolder,"/analyses/results/", name_analysis,".csv"), header = T))
labels_groups = c("+BM MSC", "+Ad MSC")
result$winlength = result$end - result$start
data.table :: setorder(result, chrInt, winlength, start)

#resmax = result[end - start == 1e08]
#setorder(resmax, chrInt)
# qplot((end + start) / 2, chrs_len[resmax$chrInt] / 2, data = resmax) + geom_smooth(method = "lm")

res_chrs = split(result, by = "chr")
res_chrs = purrr::map(res_chrs, ~ split(.x, by = "winlength"))

chrs_len = includeOrganismHuman()$infoChromosomes$length
chrs_nam = includeOrganismHuman()$infoChromosomes$idCharacter

chrs_select = 1:25
chrs_whi = intersect(which(chrs_nam %in% unique(result$chr)), chrs_select)

chrs_len = chrs_len[chrs_whi]
chrs_nam = chrs_nam[chrs_whi]

separate_chrs = 0e05
chrs_str = cumsum(c(0, as.numeric(chrs_len[-length(chrs_len)])))
chrs_str = chrs_str + (seq_len(length(chrs_len)) - 1L) * separate_chrs
chrs_end = chrs_str + chrs_len

plt_chrs = purrr :: map(res_chrs, ~ purrr :: map(.x, ~ .x[, .(start, end, signedTestStat)]))
plt_chrs = purrr :: map2(plt_chrs, chrs_str, function(x, y) purrr :: map(x, ~ .x + matrix(c(rep(y, 2*nrow(.x)), rep(0, nrow(.x))), nrow = nrow(.x), ncol = 3)))
plt_chrs = do.call("rbind", purrr :: map(plt_chrs, ~ do.call("rbind", .x)))
plt_chrs$whlen = sapply(plt_chrs$end - plt_chrs$start, function(x) which(x == sort(unique(plt_chrs$end - plt_chrs$start))))

#if(!differential) plt_chrs$signedTestStat[plt_chrs$signedTestStat < 0] = 0
genTable$center = chrs_str[genTable$chr_int] + genTable$center

range_test = range(plt_chrs$signedTestStat)
data.table :: setorder(plt_chrs, whlen, start, end)

plt_chrs$whlm1 = plt_chrs$whlen - 1L

win_len = sort(unique(plt_chrs$end - plt_chrs$start))
max_range = max(abs(range_test))
threshold = min(result$testStat[result$adjpvalue < 0.05])
    
shift_labels = 3e07
size_legend = 2e07
shift_legend = 3e07
length_legend = 101

if(differential) {
  #color_scale <- scales::col_numeric(palette = c("blue", "whitesmoke", "red"), domain = c(-1,1))
  color_scale <- scales::col_numeric(palette = c("#72ba78", "grey99", "#0054ff"), domain = c(-1,1)) #008b0b
  #color_scale <- scales::col_numeric(palette = c("red2", "grey99", "black"), domain = c(-1,1))
  plt_chrs[, cols := color_scale(signedTestStat / max(abs(signedTestStat)))]
  grid_legend = seq(from = -1, to = 1, length.out = length_legend)
  
  data_legend = data.table :: data.table(xmin = rep(tail(chrs_end,1) + shift_legend, length_legend), 
                                         xmax = rep(tail(chrs_end,1) + shift_legend + size_legend, length_legend), 
                                         ymin = length(win_len) * (0:(length_legend-1)) / length_legend, 
                                         ymax = length(win_len) * (1:(length_legend-0)) / length_legend, 
                                         col  = color_scale(grid_legend))
  
  max_legend = as.integer(substr(as.character(max_range), 0, 1)) * 10 ^ (nchar(as.character(as.integer(max_range))) - 1)
  labels_legend = data.table :: data.table(x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                           y = c(0,length(win_len))[1+differential] / 2 + list(0:4/4, (-2:2)/4)[[1+differential]] * length(win_len) * max_legend / max_range,
                                           l = list(0:4/4, (-2:2)/2)[[1+differential]] * max_legend,
                                           xmin = tail(chrs_end,1) + shift_legend,
                                           xmax = tail(chrs_end,1) + shift_legend + size_legend)
  
  data_threshold = data.table :: data.table(xmin = tail(chrs_end,1) + shift_legend,
                                            xmax = tail(chrs_end,1) + shift_legend + size_legend,
                                            x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                            y = c(list(numeric(), 0)[[1+differential]], c(0,length(win_len))[1+differential] / 2 + list(1, c(-1,1)/2)[[1+differential]] * length(win_len) * threshold / max_range, length(win_len)),
                                            l = c(list(character(), labels_groups[2])[[1+differential]], rep("pv < 0.05", c(1,2)[1+differential]), labels_groups[1]))
}
if(!differential) {
  color_scale <- scales::col_numeric(palette = c("aquamarine4", "whitesmoke", "black"), domain = c(-1, 1))
  plt_chrs[, cols := color_scale(signedTestStat / abs(range_test[1 + (signedTestStat > 0)]))]
  # plot(plt_chrs$signedTestStat, plt_chrs$signedTestStat, col = plt_chrs$cols)
  # plt_chrs[which.min(signedTestStat)]
  grid_legend = c(-10:-1 / 10, 0, 1:90 / 90) #seq(from = range_test[1], to = range_test[2], length.out = length_legend)
  length_legend = length(grid_legend)
  
  data_legend = data.table :: data.table(xmin = rep(tail(chrs_end,1) + shift_legend, length_legend), 
                                         xmax = rep(tail(chrs_end,1) + shift_legend + size_legend, length_legend), 
                                         ymin = length(win_len) * (0:(length_legend-1)) / length_legend, 
                                         ymax = length(win_len) * (1:(length_legend-0)) / length_legend, 
                                         col  = color_scale(grid_legend))
  
  max_legend = as.integer(substr(as.character(max_range), 0, 1)) * 10 ^ (nchar(as.character(as.integer(max_range))) - 1)
  labels_legend = data.table :: data.table(x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                           y = c(0,length(win_len))[1+differential] / 2 + list(0:4/4, (-2:2)/4)[[1+differential]] * length(win_len) * max_legend / max_range,
                                           l = (0:4/4) * max_legend,
                                           xmin = tail(chrs_end,1) + shift_legend,
                                           xmax = tail(chrs_end,1) + shift_legend + size_legend)
  
  data_threshold = data.table :: data.table(xmin = tail(chrs_end,1) + shift_legend,
                                            xmax = tail(chrs_end,1) + shift_legend + size_legend,
                                            x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                            y = c(list(numeric(), 0)[[1+differential]], c(0,length(win_len))[1+differential] / 2 + list(1, c(-1,1)/2)[[1+differential]] * length(win_len) * threshold / max_range, length(win_len)),
                                            l = c(list(character(), labels_groups[2])[[1+differential]], rep("pv < 0.05", c(1,2)[1+differential]), labels_groups[1]))
  
}

plt_chrs[, significant := abs(signedTestStat) >= threshold]

max_legend = as.integer(substr(as.character(max_range), 0, 1)) * 10 ^ (nchar(as.character(as.integer(max_range))) - 1)
labels_legend = data.table :: data.table(x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                         y = c(0,length(win_len))[1+differential] / 2 + list(0:4/4, (-2:2)/4)[[1+differential]] * length(win_len) * max_legend / max_range,
                                         l = list(0:4/4, (-2:2)/2)[[1+differential]] * max_legend,
                                         xmin = tail(chrs_end,1) + shift_legend,
                                         xmax = tail(chrs_end,1) + shift_legend + size_legend)

data_threshold = data.table :: data.table(xmin = tail(chrs_end,1) + shift_legend,
                                          xmax = tail(chrs_end,1) + shift_legend + size_legend,
                                          x = tail(chrs_end,1) + shift_legend + size_legend + shift_labels,
                                          y = c(list(numeric(), 0)[[1+differential]], c(0,length(win_len))[1+differential] / 2 + list(1, c(-1,1)/2)[[1+differential]] * length(win_len) * threshold / max_range, length(win_len)),
                                          l = c(list(character(), labels_groups[2])[[1+differential]], rep("pv < 0.05", c(1,2)[1+differential]), labels_groups[1]))

if(differential)  sig_chrs = plt_chrs[significant == TRUE, .(start, end, y = (whlen + whlm1) / 2, col = color_scale(sign(signedTestStat)))]
if(!differential) sig_chrs = plt_chrs[significant == TRUE, .(start, end, y = (whlen + whlm1) / 2, col = color_scale(ifelse(signedTestStat>0, range_test[2], range_test[1])))]

p2 = ggplot2 :: ggplot()
p2 = p2 + ggplot2 :: geom_vline(xintercept = unique(c(chrs_str, chrs_end)), color = "grey30", linewidth = .2)
p2 = p2 + ggplot2 :: geom_vline(data = genTable, mapping = ggplot2 :: aes(xintercept = center, color = posstat)) + ggplot2 :: scale_color_identity()
p2 = p2 + ggplot2 :: geom_rect(mapping = ggplot2 :: aes(xmin = 0 - 1e07, xmax = max(chrs_end) + 1e07, ymin = 0, ymax = max(plt_chrs$whlen)), fill = "white") #whitesmoke, grey
p2 = p2 + ggplot2 :: geom_rect(data = plt_chrs, mapping = ggplot2 :: aes(xmin = start, xmax = end, ymin = whlm1, ymax = whlen, fill = cols)) + ggplot2 :: scale_fill_identity()
#p2 = p2 + geom_segment(data = sig_chrs, mapping = aes(x = start, xend = end, y = y, yend = y), size = 1.2, color = "white")
#p2 = p2 + geom_segment(data = sig_chrs, mapping = aes(x = start, xend = end, y = y, yend = y), size = 0.8, color = "black")
if(differential) {
  p2 = p2 + ggplot2 :: geom_point(data = sig_chrs, mapping = ggplot2 :: aes(x = (start + end) / 2, y = y), size = 1.2, color = "white")
  #p2 = p2 + geom_point(data = sig_chrs, mapping = aes(x = (start + end) / 2, y = y), size = 1.0, color = "#f40000")
  p2 = p2 + ggplot2 :: geom_point(data = sig_chrs, mapping = ggplot2 :: aes(x = (start + end) / 2, y = y, colour = col), size = 0.8) + ggplot2 :: scale_color_identity()
}
p2 = p2 + ggplot2 :: geom_rect(data = data_legend, mapping = ggplot2 :: aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col)) + ggplot2 :: scale_fill_identity()
p2 = p2 + ggplot2 :: geom_text(data = labels_legend, mapping = ggplot2 :: aes(x = x, y = y, label = l), size = 3)
p2 = p2 + ggplot2 :: geom_segment(data = labels_legend, mapping = ggplot2 :: aes(x = xmin, y = y, xend = xmax), color = "white")
p2 = p2 + ggplot2 :: geom_text(data = data_threshold, mapping = ggplot2 :: aes(x = x, y = y, label = l), size = 2.5, hjust = 0.2)
p2 = p2 + ggplot2 :: geom_segment(data = data_threshold, mapping = ggplot2 :: aes(x = xmin, y = y, xend = xmax), color = "black")

p2 = p2 + ggplot2 :: scale_x_continuous(breaks = head(c((chrs_str + chrs_end) / 2, tail(chrs_end,1) + shift_legend + size_legend / 2), -1), labels = c(1:22, "X", "Y"))
p2 = p2 + ggplot2 :: scale_y_continuous(breaks = 1:length(win_len) - .5, labels = format(win_len, scientific = T, digits = 3))
p2 = p2 + ggplot2 :: xlab("Chromosomes") + ggplot2 :: ylab("Length tested windows (bps)")
p2 = p2 + ggplot2 :: theme(legend.position = "none", 
                           panel.background = ggplot2 :: element_blank(), 
                           axis.ticks.x = ggplot2 :: element_blank(), 
                           axis.ticks.y = ggplot2 :: element_blank(), 
                           panel.grid.major.x = ggplot2 :: element_blank(), 
                           panel.grid.major.y = ggplot2 :: element_blank(), 
                           panel.grid.minor.y = ggplot2 :: element_blank(), 
                           panel.grid.minor.x = ggplot2 :: element_blank(),
                           axis.title.y = ggplot2::element_text(size = 14))
#p2
