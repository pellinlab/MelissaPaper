mainFolder = "~/Downloads/MELISSApaper"

library(data.table)
library(ggplot2)
library(purrr)

idType = c("govt", "gdit", "clfi", "cdif")
idSize = c("s01", "s02", "s04", "s08")
idHeig = c("h01", "h02", "h04", "h08", "h16", "h32", "h64")
idGene = c("g20", "g10", "g80", "g40", "g05")

ppv = cov = fwt = vector("list", length(idType))
names(ppv) = names(cov) = names(fwt) = idType

h=i=j=1L
for(h in 1:4) {
  ppv[[idType[h]]] = vector("list", length(idSize))
  names(ppv[[idType[h]]]) = idSize
  for(i in 1:4) {
    idFile = paste0(mainFolder, "/simulation/resu/ppv/", idType[h], "_", idSize[i], ".csv")
    tab = read.csv(idFile, header = T)
    ppv[[idType[h]]][[idSize[i]]] = lapply(tab, identity)
  }  
}
for(h in 1:4) {
  cov[[idType[h]]] = fwt[[idType[h]]] = vector("list", length(idSize))
  names(cov[[idType[h]]]) = names(fwt[[idType[h]]]) = idSize
  for(i in 1:4) {
    cov[[idType[h]]][[idSize[i]]] = fwt[[idType[h]]][[idSize[i]]] = vector("list", length(idGene))
    names(cov[[idType[h]]][[idSize[i]]]) = names(fwt[[idType[h]]][[idSize[i]]]) = idGene
    for(j in 1:5) {
      idFile = paste0(mainFolder, "/simulation/resu/cov/", idType[h], "_", c("g1","g2","g3","g4","g5")[j], "_", idSize[i], ".csv")
      tab = read.csv(idFile, header = T)
      names(tab) = idHeig
      cov[[idType[h]]][[idSize[i]]][[idGene[j]]] = lapply(tab, identity)
      idFile = paste0(mainFolder, "/simulation/resu/cov/", idType[h], "_", c("f1","f2","f3","f4","f5")[j], "_", idSize[i], ".csv")
      tab = read.csv(idFile, header = T)
      names(tab) = idHeig
      fwt[[idType[h]]][[idSize[i]]][[idGene[j]]] = lapply(tab, identity)
    }  
  }  
}
rm(idFile, tab)

meanStat = \(x) unlist(x |> map(~ map_dbl(.x, ~ mean(.x, na.rm = T))))
sumStat = \(x) (x |> map(~.x |> map_dbl(sum)))
titlesPlot = c("A) gene over-targeting", "B) differential gene targeting", "C) clone fitness", "D) differential clone fitness")
breaksy = (0:5)/5
ordGen = c(3,2,5,4,1)

h=1
for(h in 1:length(idType)) {
  pp = ppv[[idType[h]]]
  co = cov[[idType[h]]]
  fw = fwt[[idType[h]]]
  
  isCloneFitness = h >= 3L
  isDifferential = h %in% c(2L, 4L)
  idLegend = paste0("sample size  [ ", list("", "2 x ")[[1 + isDifferential]], list("4 x", "6 x")[[1 + isCloneFitness]], " ]")
  
  db = data.frame(ppv = meanStat(pp))
  db$x = log2(as.numeric(substring(idHeig, 2)))
  db$s = rep(as.character(c(1, 2, 4, 8) * c(100, 1000)[1 + isCloneFitness]), each  = length(idHeig))
  db$h = rep(idHeig, times = length(idSize))
  p = ggplot(db)
  p = p + scale_x_continuous(breaks = log2(c(1, 2, 4, 8, 16, 32, 64)), labels = c(1, 2, 4, 8, 16, 32, 64), minor_breaks = NULL)
  p = p + scale_y_continuous(breaks = breaksy, minor_breaks = NULL)
  p = p + labs(x = "effect size h (log2 scale)", 
               y = "<span style = 'color:#0072B2;'>PPV</span>", 
               title = titlesPlot[h], 
               linetype = idLegend)
  p = p +  theme(axis.title.y = ggtext::element_textbox_simple(width = NULL, orientation = "left-rotated"),
                 panel.background = element_rect(fill = NA),
                 panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 legend.position = c("none","bottom")[2]
  )
  p = p + geom_line(aes(x = x, y = ppv, linetype = s), color = "#0072B2")
  p = p + geom_point(aes(x = x, y = ppv), color = "#0072B2")
  p = p + guides(linetype = guide_legend(override.aes = list(colour = "black")))
  ggsave(paste0(mainFolder, "/figures/f2/sim_ppv_", idType[h], ".pdf"), width = 5.0, height = 4.0, units = "in")
  
  
  db = data.frame(cov = unlist(lapply(co, sumStat)), fwt = unlist(lapply(fw, sumStat)), gen = rep(rep(idGene, each = 7), 4))
  db$x = rep(rep(ordGen, each = 7), 4)
  db$s = rep(as.character(c(1, 2, 4, 8) * c(100, 1000)[1 + isCloneFitness]), each  = 7 * 5) 
  db$h = rep(rep(idHeig, times = length(idSize)), 5)
  
  # choose size 400 / 4000
  db = db[which(db$s == 4*c(100, 1000)[1 + isCloneFitness]),]
  db$x = db$x + rep((-3:3)/10, 5)
  db$cl = rep(rainbow(20)[10 + 1:7], 5)
  
  p = ggplot(db)
  for(i in 1:nrow(db)) p = p + annotate("segment", x = db$x[i], xend = db$x[i], y = 0, yend = db$fwt[i], colour = adjustcolor(db$cl[i], alpha = .2), linewidth = 2.5)
  for(i in 1:nrow(db)) p = p + annotate("segment", x = db$x[i], xend = db$x[i], y = 0, yend = db$cov[i], colour = db$cl[i], linewidth = 2.5)
  p = p + scale_x_continuous(breaks = c(1,2,3,4,5), labels = c("5kb", "10kb", "20kb", "40kb", "80kb"), minor_breaks = NULL)
  p = p +  theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 #axis.title.y = ggtext::element_textbox_simple(width = NULL, orientation = "left-rotated"),
                 legend.position = "none"
                 #legend.position = "bottom"
  )
  p = p + labs(x = "length target gene", y = "coverage")
  
  ggsave(paste0(mainFolder, "/figures/f2/sim_cov_", idType[h], ".pdf"), width = 5.0, height = 4.0, units = "in")
}




