six = six[db == "Sixetal"]
six = six[seqnames == "chr12"]
six = six[start >= 64700000] # 12q14.3
six = six[start <  67300000] # 12q14.3
six$inHmga2 = six$start >= 65824130 & six$start < 65966295
six = six[patient == "b0_bE"]
six = six[six$inHmga2]
table(six$timePointMonths)
six$timePointMonths[six$timePointMonths == 11L] = 12L

rav = rav[cell_type != "HSPC"]
rav = rav[seqnames == "chr12"]
rav = rav[start >= 64700000]
rav = rav[start <  67300000]
rav$inHmga2 = rav$start >= 65824130 & rav$start < 65966295
rav = rav[rav$inHmga2]
table(rav$patient, rav$time_mt)
rav$time_mt[rav$time_mt == 35L] = 36L
rav = rav[rav$time_mt %in% c(12L, 36L, 48L)]

s = six[, .(pos = start, str = strand, cnt = estAbund, pat = "B-thal", time = timePointMonths)]
s = s[, .(tot = sum(cnt)), by = .(pos, str, pat, time)]
#r = rav[, .(pos = start, str = strand, cnt = count, pat = purrr::map_chr(rav$patient, ~ paste0("SCID-X1_", .x)), time = time_mt)]
r = rav[, .(pos = start, str = strand, cnt = count, pat = purrr::map_chr(rav$patient, ~ paste0("", .x)), time = time_mt)]
r = r[, .(tot = sum(cnt)), by = .(pos, str, pat, time)]

p = rbind(s, r)
p$pos = p$pos * c(-1, 1)[1 + (p$str == "+")]
data.table::setorder(p, pos)
p = split(p, by = "pat")
p = purrr::map(p, split, by = "time")
p = purrr::map(p, ~ purrr::map(.x, split, by = "str"))

calc_csm = function(x) {
  strplus = x$str[1] == "+"
  x0 = x[c(1,1),]
  x0$pos = (c(65824130L, 65966295L) * c(-1L, 1L)[1 + strplus])[list(2:1, 1:2)[[1 + strplus]]]
  x0$tot = 0L
  
  x = rbind(x0[1], x, x0[2])
  x$pos = as.integer(abs(x$pos))
  x$csm = cumsum(x$tot)
  x$ecd = x$csm / tail(x$csm, 1)
  return(x)
}

p = purrr::map(p, ~ purrr::map(.x, ~ purrr::map(.x, calc_csm)))
p = purrr::map(p, ~ purrr::map(.x, ~ do.call("rbind", .x)))

rescale_ecd = function(x) {
  xm = x[str == "-"]
  xp = x[str == "+"]
  pp = max(c(xp$csm, 0)) / (max(c(xp$csm,0)) + max(c(xm$csm,0)))
  xp$ecd = xp$ecd * pp
  xm$ecd = max(xp$ecd) + xm$ecd * (1-max(xp$ecd))
  
  x = rbind(xp, xm)
  return(x)
}

p = purrr::map(p, ~ purrr::map(.x, rescale_ecd))

p = purrr::map(p, ~ do.call("rbind", .x))
p = do.call("rbind", p)

#p = p[str == "+"]

p$pat = factor(p$pat, levels = sort(unique(p$pat)))
p$time = factor(p$time, levels = c("48", "36", "12"), labels = c("48m", "36m", "12m"))
p$isbtl = p$pat == "B-thal"

colPatients = c(6, NA, 5, 8, 2, 4, 3, 7)[c(1,3,4,5,6,7,8)]
colsLabel = purrr::map(colPatients, ~ colorRampPalette(c(.x, "black")))
colsLabel = purrr::map_chr(colsLabel, ~ .x(10)[3])

annotations = rtracklayer::import(fileAnnotations)

transcripts = data.table::as.data.table(subset(annotations, gene_name == "HMGA2" & type == "exon"))
ggplot2 :: ggplot() + #data = p, aes(x = pos, y = ecd, color = pat, group = time)) + 
  ggplot2 :: geom_rect(mapping = ggplot2 :: aes(xmin = start, xmax = end, ymin = 0, ymax = 1), data = transcripts, color = "grey90", fill = "grey90") +
  #geom_vline(xintercept = c(65824130L, 65966295L), color = "grey90") + 
  ggplot2 :: geom_point(ggplot2 :: aes(x = pos, y = ecd, color = pat), data = subset(p, !(pos %in% (c(65824130L, 65966295L))) & (str == "-")), size = .6) +
  ggplot2 :: geom_step(ggplot2 :: aes(x = pos, y = ecd, color = pat, group = pat, linetype = isbtl), data = subset(p, str == "-"), direction = "vh") + 
  ggplot2 :: geom_point(ggplot2 :: aes(x = pos, y = ecd, color = pat), data = subset(p, !(pos %in% (c(65824130L, 65966295L))) & (str == "+")), size = .8) +
  ggplot2 :: geom_step(ggplot2 :: aes(x = pos, y = ecd, color = pat, group = pat, linetype = isbtl), data = subset(p, str == "+")) + 
  ggplot2 :: facet_wrap(~ time + str, ncol = 1) +
  #ggplot2 :: facet_wrap(~ time, ncol = 1) +
  ggplot2 :: scale_color_manual(values = c("#3D1255", colsLabel)) +
  ggplot2 :: scale_x_continuous(expand = c(0, 0)) + 
  ggplot2 :: coord_cartesian(xlim = c(65824130L, 65966295L)) +
  #geom_rect(mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), data = transcripts, color = NA, fill = "grey90") +
  ggplot2 :: theme(legend.position = "bottom",
        #axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        #axis.text.x = ggplot2::element_blank(), 
        #axis.ticks.x = ggplot2::element_blank(),
        #axis.text.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = NA),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor.y = ggplot2::element_blank()
  ) + ggplot2 :: guides(colour = ggplot2 :: guide_legend(nrow = 1)) +
  ggplot2 :: ylab("cumulative proportion strand + clones in HMGA2") -> pl
pl

one_transcript_plotted = T
if(one_transcript_plotted) {
  tps = split(transcripts, by = "transcript_id")
  rtp = purrr::map(tps, ~ range(unlist(.x[,.(start, end)])))
  rgn = c(65824130L, 65966295L)
  wtp = which.min(purrr::map_int(rtp, ~ sum(abs(.x - rgn))))
  transcripts = tps[[wtp]]
}

p2 = ggplot2 :: ggplot(data = transcripts) + ggplot2 :: geom_vline(xintercept = c(65824130L, 65966295L), color = "grey90")
#if(one_transcript_plotted) {
  p2 = p2 + 
    #ggtranscript :: geom_range(aes(xstart = start, xend = end, y = transcript_id, fill = NULL, colour = NULL)) +
    ggtranscript :: geom_intron(data = ggtranscript :: to_intron(transcripts, group_var = "transcript_id"), ggplot2 :: aes(xstart = start, xend = end, y = transcript_id, strand = strand)) +
    ggtranscript :: geom_range(ggplot2 :: aes(xstart = start, xend = end, y = transcript_id, fill = NULL, colour = NULL))
#} else {
#   p2 = p2 + 
#     ggtranscript :: geom_intron(data = ggtranscript :: to_intron(transcripts, group_var = "transcript_name"), aes(xstart = start, xend = end, y = transcript_name, strand = strand)) +
#     ggtranscript :: geom_range(aes(xstart = start, xend = end, y = transcript_name, fill = transcript_type, colour = transcript_type))
# }
p2 = p2 + ggplot2 :: coord_cartesian(xlim=c(65824130L, 65966295L)) +
  ggplot2 :: ylab(NULL) +
  ggplot2 :: scale_x_continuous(expand = c(0, 0)) + 
  ggplot2 :: theme(plot.background = ggplot2 :: element_rect(fill = 'white', color = "white"),
        panel.background = ggplot2 :: element_rect(fill = 'white', color = "white"),
        axis.text.x = ggplot2 :: element_blank())
#p2

