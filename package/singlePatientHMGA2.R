idpt = list.files(path = resuFolder, pattern = "^cf_ov_p.*\\.csv$", full.names = F)
pt = purrr::map_chr(idpt, substr, start = 7, stop = 8)
tm = as.integer(sub(".*_tmax(\\d+).*", "\\1", idpt))

results = purrr::map(paste0(resuFolder, idpt), read.table, header = T, sep = "\t")
results = purrr::map(results, data.table::as.data.table)
infodata = data.table::data.table(file = idpt, patient = pt, time = tm)
rm(idpt, pt, tm)

re = purrr::map(results, ~ .x[,.(geneName, testStat, signedTestStat, adjpvalue, robustTestStat, robustSignedTestStat, robustAdjPvalue)])
re = purrr::map(re, ~ .x[testStat >= 0])
invisible(purrr::map(re, ~ data.table::setorder(.x, -signedTestStat)))

rp = purrr::map(re, ~ .x[signedTestStat >= 1])
rp = purrr::map(rp, ~ .x[geneName != "HMGA2"])

dp = cbind(infodata[unlist(purrr::map2(1:length(rp), rp, ~ rep(.x, nrow(.y))))], 
           do.call("rbind", purrr::map(rp, ~ .x[,.(geneName, signedTestStat)])))

infodata$thr5m2 = purrr::map_dbl(re, ~ min(.x$testStat[.x$adjpvalue <= 0.05]))
#infodata$thr1m4 = purrr::map_dbl(re, ~ min(.x$testStat[.x$adjpvalue <= 0.0001]))
infodata$hmga2 =  purrr::map_dbl(re, ~ .x$signedTestStat[.x$geneName == "HMGA2"])

data.table::setorder(infodata, patient, time)

infodata$hmga2[infodata$hmga2 < 1] = 1

infodata$time[c(20, 37, 53, 61)] = c(54L, 54L, 42L, 38L)
dp$time[dp$patient == "p3" & dp$time == 60L] = 54L
dp$time[dp$patient == "p5" & dp$time == 60L] = 54L
dp$time[dp$patient == "p7" & dp$time == 48L] = 42L
dp$time[dp$patient == "p8" & dp$time == 48L] = 38L


setTheme = function(p, margins = ggplot2::theme_grey()$plot.margin) p + ggplot2::theme(
  legend.position = "none",
  #axis.title.y = ggplot2::element_blank(),
  #axis.title.x = ggplot2::element_blank(),
  #axis.text.x = ggplot2::element_blank(), 
  #axis.ticks.x = ggplot2::element_blank(),
  #axis.text.y = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_text(size = 14),
  axis.title.x = ggplot2::element_text(size = 14),
  panel.background = ggplot2::element_rect(fill = NA),
  panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
  panel.grid.minor.x = ggplot2::element_line(colour = "grey90"),
  panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
  panel.grid.minor.y = ggplot2::element_line(colour = "grey90")
  #plot.margin = unit(margins, "pt")
)

log10_minor_break = function (...) {
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

p = ggplot2::ggplot(data = infodata, mapping = ggplot2::aes(x = as.factor(time), group = patient)) +
  ggplot2::geom_jitter(ggplot2::aes(x = as.factor(time), y = signedTestStat, group = patient), width = .2, height = 0, data = dp, pch = 1, color = "grey70") +
  ggplot2::geom_path(ggplot2::aes(y = hmga2)) +
  ggplot2::geom_point(ggplot2::aes(y = hmga2)) +
  ggplot2::geom_tile(ggplot2::aes(y = thr5m2), height = .01, color = "#f40000", fill = "#f40000") +
  ggplot2::facet_grid(rows = . ~ patient, scales = "free") +
  #scale_x_continuous(breaks = c(6,12,18,24,36,48,60,72,84), minor_breaks = setdiff(3 * (1:28), c(6,12,18,24,36,48,60,72,84))) +
  ggplot2::scale_y_continuous(breaks = 10 ^ (0:6), minor_breaks = log10_minor_break(), trans = c("log10")) +
  ggplot2::xlab("max included time") + ggplot2::ylab("positive fitness score")
p = setTheme(p)
#p
