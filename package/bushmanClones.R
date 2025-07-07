plotOneChromosome = F
whichChrInt = 12 # 9

whichPatients = c(1, 3, 4, 5, 6, 7, 8) # c(5, 6, 7) # 
colPatients = c(6, NA, 5, 8, 2, 4, 3, 7)[whichPatients]
patients = paste0("p", whichPatients)
filesResult = paste0(mainFolder, "/analyses/results/", "cf_ov_", patients, "_tall_winclones.csv")
resultMelissa = purrr::map(filesResult , data.table::fread, sep="\t", header = T)


if(plotOneChromosome) {
  resuChrs = map(resultMelissa, ~ split(.x, f = as.character(.x$chrInt)))
  resuChrs = map(resuChrs, ~ .x[[as.character(whichChrInt)]])
} else resuChrs = resultMelissa

log10scale = T
differential = F
apvThreshold = 0.1
numLabelClones = 5

resuChrs = purrr :: map(resuChrs, ~ .x[, whichPositive := signedTestStat >= 0])
resuChrs = purrr :: map(resuChrs, ~ .x[, whichBlwTreshold := adjpvalue <= apvThreshold])
if(!differential) resuChrs = purrr :: map(resuChrs, ~ .x[, whichBlwTreshold := whichBlwTreshold & whichPositive])

plt = purrr :: map(resuChrs, ~ subset(.x, whichBlwTreshold))
plt = purrr :: map(plt, ~ .x[, plottedStat := testStat * whichPositive])
plt = purrr :: map(plt, ~ data.table :: setorder(.x, start, -strand))
plt = purrr :: map(plt, ~ .x[, cloneIndex := (0:(nrow(.x)-1))])
plt = purrr :: map(plt, ~ data.table :: setorder(.x, -plottedStat, -strand))
plt = purrr :: map(plt, ~ .x[, statIndex  := (0:(nrow(.x)-1))])
plt = purrr :: map(plt, ~ .x[, plotLabelClone := rep(c(T,F), times = c(min(numLabelClones, nrow(.x)), max(nrow(.x)-numLabelClones, 0)))])
plt = purrr :: map2(plt, patients, ~ .x[, patient := .y])

bpsChrhg38 = c(0, 248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417, 3088269832, 3088286401)[-26]

shiftTss = function(tss, chr) {
  bps = bpsChrhg38
  return(tss + bps[chr])
}


pl_globSizeTextRepel = 3
pl_globOverlapTextRepel = 15
pl_axisTitleSize = 14


# pl1
pl_shiftXcoord = (0 : 6) #/ 2

plt = purrr :: map(plt, ~ data.table :: setorder(.x, statIndex))
plt = purrr :: map2(plt, pl_shiftXcoord, ~ .x[,xcoord := (statIndex / (nrow(.x)-1) + .y)])

pl = do.call("rbind", plt)
pl = pl[testStat >= 0]
pl$patient = factor(pl$patient)
pl$basepair = shiftTss(pl$start, pl$chrInt)
pl$hmga2 = (pl$start >= 65824130 & pl$start < 65966295 & pl$chrInt == 12L)

data.table :: setorder(pl, -plottedStat)

pl$adjpvalue = pl$adjpvalue + 1e-16

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

colsLabel = purrr :: map(colPatients, ~ colorRampPalette(c(.x, "black")))
colsLabel = purrr :: map_chr(colsLabel, ~ .x(10)[3])

pl$Patient = pl$patient
ggplot2::ggplot(data = pl, ggplot2::aes(x = xcoord, y = plottedStat)) +
  ggplot2::geom_point(ggplot2::aes(size = I(.1 - 2 * log10(adjpvalue) / 16), color = Patient)) +
  ggplot2::geom_point(data = subset(pl, hmga2), mapping = ggplot2::aes(size = I(.1 - 2 * log10(adjpvalue) / 16)), color = "black", shape = 21, fill = "white") +
  ggrepel::geom_text_repel(data = subset(pl, plotLabelClone), ggplot2::aes(label = geneName, color = patient),
                           size = pl_globSizeTextRepel, max.overlaps = Inf, direction = "both") +
  ggplot2::scale_color_manual(values = colsLabel) +
  ggplot2::scale_y_log10(breaks = 10 ^ (1:6), minor_breaks = log10_minor_break()) +
  ggplot2::scale_x_continuous(breaks = NULL) +
  ggplot2::coord_cartesian(ylim = range(pl$plottedStat)) +
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA),
                 axis.title.x = ggplot2::element_blank(),
                 #axis.text.x = ggplot2::element_blank(), 
                 axis.ticks.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                 panel.grid.minor.y = ggplot2::element_line(colour = "grey90")) +
  ggplot2::ylab("Individual clone growth statistic") -> p1
#p1

whichpv = c(1, 4, 7, 10, 13)
adjpv = (10 ^ -(1:15) + 1e-16)[whichpv]
ggplot2::ggplot(data.table::data.table(xax = 1, yax = (1:15)[whichpv], `Adjusted\n pvalue` = as.character(adjpv - 1e-16), `Clone in\n HMGA2` = c("HMGA2", rep("in", length(whichpv) - 1)))) +
  ggplot2::geom_point(ggplot2::aes(xax, yax, size = `Adjusted\n pvalue`, color = `Clone in\n HMGA2`)) +
  ggplot2::scale_color_manual(values = c(1,2)) +
  ggplot2::scale_size_manual(values = .1 - 2 * log10(adjpv) / 16) +
  #geom_text(x = 1.025, y = 5, label = "sdf", size = ) +
  ggplot2::theme_minimal() -> pleg
