FocusGenesIntegrationPlot = function(fileGenes, folderStoredPlots, 
                                      integrSiteObject, includeOnly,
                                      aes_y_char, aes_shape_char, aes_color_char,# aes_size_char,
                                     selectedCovariates = character(),
                                     selectedFactorCovariates = character(),
                                     covariatesInLabels = "time",
                                     ylab = "Passage",
                                     formatPlot = c(".png", ".svg", ".pdf")[1],
                                     showLegend = FALSE,
                                     transparency = 1,
                                     jitter = 0,
                                     colors = list(),
                                     sizeInchesPlots = c(9, 2)) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  if(length(iSO$data$combined) != 0L) stop("error") # not yet implemented
  if(!missing(includeOnly)) {
    iSO$design$filesIntegrations = it = it[includeOnly]
    iSO$design$covariatesIntegrations = cv = cv[includeOnly,]
    iSO$design$id = 1:length(it)
  }
  
  #selectedCovariates = unique(c(aes_y_char, aes_shape_char, aes_color_char, selectedCovariates))#, aes_size_char)
  if(!all(selectedCovariates %in% c("chr", "strand", "count", names(cv)))) stop("error")
  #if(!all(selectedFactorCovariates %in% selectedCovariates)) stop("error")
  
  fg = data.table :: fread(fileGenes, sep = "\t", header = T)
  
  iSO = setGroups(iSO, group1 = 1:length(it), group2 = integer(0L))
  iSO = setInfoSelectedCovariates(iSO, selectedCovariates = selectedCovariates, selCovAreFactors = selectedFactorCovariates)
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  iSO = loadLocationsGenes(iSO)
  iSO = mergeData(iSO)
  #iSO$design$selectedCovariates = selectedCovariates
  #iSO$design$selCovAreFactors   = sapply(selectedCovariates, function(x) x %in% selectedFactorCovariates)
  #class(iSO)
  #iSO = setFactorsCombinedData(iSO)
  #iSO = mergeData(iSO, selectedCovariates[which(!(selectedCovariates %in% c("chr", "start", "end", "strand", "group")))])
  #iSO = setFactorsCombinedData(iSO, variables = c(c("chr", "strand", "group"), selectedCovariates[which(selectedCovariates %in% selectedFactorCovariates)]))
  
  selectedFactorCovariates = c("strand", selectedCovariates[selectedFactorCovariates])
  if(length(selectedFactorCovariates) > 0L) for(i in 1:length(selectedFactorCovariates)) iSO$data$combined[[selectedFactorCovariates[i]]] = as.factor(iSO$data$combined[[selectedFactorCovariates[i]]])
  
  iSO$data$combined$relativeAbundance = iSO$data$combined$count / iSO$data$combined$totalCount
  aes_size_char = "relativeAbundance"
  
  textSize = 16
  #colorCell = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909")
  
  i = 1L
  for(i in 1:nrow(fg)) {
    da = iSO$data$combined[chr == fg$chr[i],]
    #da = iSO$data$combined[start >= fg$start[i] & end < fg$end[i],]
    da = da[start >= fg$start[i] & end < fg$end[i],]
    p = ggplot2 :: ggplot(data = da)
    #p = p + geom_point(data = da, aes(x = pos, y = timePointMonths, shape = strand, color = cellType, size = relAbund ))
    p = p + ggplot2 :: geom_point(ggplot2 :: aes(x = start, 
                                                 y = .data[[aes_y_char]], 
                                                 shape = .data[[aes_shape_char]], 
                                                 color = .data[[aes_color_char]], 
                                                 size = .data[[aes_size_char]]), 
                                  position = ggplot2 :: position_jitter(width = 0, height = jitter),
                                  alpha = transparency)
    if(length(colors)) p = p + ggplot2 :: scale_color_manual(values = as.list(colors))
    p = p + ggplot2 :: scale_shape_manual(values = c("-"="<","+"=">"))
      #scale_size_continuous(breaks = c(0.005,0.01,0.05,0.1,0.2,0.3))+
    p = p + ggplot2 :: scale_size_continuous(range = c(4, 10))
      # geom_segment(aes(x=startS,y=0.5,xend=endS,yend=0.5),size=2)+
    p = p + ggplot2 :: scale_y_discrete(labels = unique(cv[, covariatesInLabels]))
    p = p + ggplot2 :: xlim(fg$start[i],fg$end[i])
    p = p + ggplot2 :: theme(text= ggplot2 :: element_text(family="sans"),
                  axis.title.x = ggplot2 :: element_text(size =textSize),
                  axis.title.y = ggplot2 :: element_text(size =textSize), 
                  panel.background = ggplot2 :: element_rect(fill = NA), 
                  panel.grid.major.y = ggplot2 :: element_line(colour = "grey90"),
                  panel.grid.minor.y = ggplot2 :: element_line(colour =NA),
                  axis.text.x = ggplot2 :: element_text(size=textSize), 
                  axis.text.y = ggplot2 :: element_text(size=textSize))#,
                  #legend.position = "bottom")
            #    legend.position = "none",
            #    panel.grid.major.x = element_line(colour = "grey90"),
            #    panel.grid.minor.x = element_line(colour = "grey90"),
            #axis.ticks.x=element_text(size=14),
    p = p + ggplot2 :: xlab(fg$chr[i])
    p = p + ggplot2 :: ylab(ylab)
    p = p + ggplot2 :: guides(shape = ggplot2 :: guide_legend(override.aes=list(size = 3)), 
                              size  = ggplot2 :: guide_legend(override.aes=list(shape = ">"))#,
                              #color = ggplot2 :: guide_legend(override.aes=list(alpha = transparency, size = 3))
                              )
    ggplot2 :: ggsave(filename = paste0(folderStoredPlots, 
                                        c("integrations_", "legend_")[1+showLegend], 
                                        paste(fg$geneName[i],fg$chr[i],fg$start[i],fg$end[i],"2", formatPlot ,sep="_",collapse = "")), 
                      #plot = p, width = 9,height = c(2,7)[1+showLegend],units = "in")
                      plot = p, width = sizeInchesPlots[1], height = sizeInchesPlots[2], units = "in")
                      
  }
  
  return(invisible(NULL))  
}


FocusGenesTrackingPlot = function(fileGenes, folderStoredPlots, integrSiteObject, includeOnly,
                                  formatPlot = c(".png", ".svg", ".pdf")[1],
                                  covariatesInLabels = c("time"),
                                  compressUniqueCovariates = TRUE,
                                  sizeInchesPlots = c(6, 4),
                                  colors = list()) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  if(length(iSO$data$combined) != 0L) stop("error") # not yet implemented
  if(!missing(includeOnly)) {
    iSO$design$filesIntegrations = it = it[includeOnly]
    iSO$design$covariatesIntegrations = cv = cv[includeOnly,]
    iSO$design$id = 1:length(it)
  }
  
  fg = data.table :: fread(fileGenes, sep = "\t", header = T)
  
  if(compressUniqueCovariates) {
    dp = duplicated(cv)
    cv = unique(cv)
    if(any(dp)) {
      ss = which(!dp)
      es = which(c(diff(dp), -1L) == -1L)
      cv$series = purrr::map2(ss, es, ~ .x : .y)
    }
    else {
      cv$series = 1:nrow(cv)
    }
  }
  else {
    cv$series = 1:nrow(cv)
  }
  
  iSO = setGroups(iSO, group1 = 1:length(it), group2 = integer(0L))
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  iSO = loadLocationsGenes(iSO)
  iSO = mergeData(iSO)
  
  iSO$data$combined$relativeAbundance = iSO$data$combined$count / iSO$data$combined$totalCount
  
  ic = iSO$data$combined
  
  textSize = 13
  pchSize = 2
  colorCell = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909")
  minAlpha = .2
  maxAlpha = .8
  minTime = min(cv$time)
  maxTime = max(cv$time)
  
  i = 1L
  for(i in 1:nrow(fg)) {
    totAb = rep(0, nrow(cv))
    relAb = replicate(n = nrow(cv), expr = numeric(), simplify = "list")
    for(j in 1:nrow(cv)) {
      #integrInGene_i = ic[which(ic$series == j & ic$chr == fg$chr[i] & ic$start >= fg$start[i] & ic$start < fg$end[i]),]
      integrInGene_i = ic[which((ic$series %in% cv$series[[j]]) & (ic$chr == fg$chr[i]) & (ic$start >= fg$start[i]) & (ic$start < fg$end[i])),]
      if(nrow(integrInGene_i) > 0L) {
        relAb[[j]] = integrInGene_i$relativeAbundance
        totAb[j] = sum(relAb[[j]])
      }
    }
    lenAb = sapply(relAb, length)
    selAb = rep(1:nrow(cv), times = lenAb)
    #alpAb = unlist(lapply(split(totAb, cv$type), function(x) (x - min(x)) / (max(x) - min(x))))[order(cv$type)]
    alpAb = sapply(cv$time, function(x) (x - minTime) / (maxTime - minTime))
    
    #d1 = data.frame(series = 1:nrow(cv), total = totAb, type = cv$type, labels = paste(cv$type, cv$time), alpha = alpAb)
    d1 = data.frame(series = 1:nrow(cv), total = totAb, type = cv$type, 
                    labels = apply(cv[,covariatesInLabels], 1, function(x) paste(x, collapse = " ")), 
                    alpha = alpAb)
    d2 = data.frame(series = selAb, relative = unlist(relAb), type = cv$type[selAb]) 
    
    p = ggplot2 :: ggplot(data = d1)
    p = p + ggplot2 :: geom_col(ggplot2 :: aes(x = series, y = total, fill = type, alpha = alpha))
    if(length(colors) > 0) p = p + ggplot2 :: scale_fill_manual(values = as.list(colors))
    p = p + ggplot2 :: scale_alpha_continuous(range = c(minAlpha, maxAlpha))
    p = p + ggplot2 :: scale_x_continuous(breaks = 1:nrow(cv), labels = d1$labels)
    p = p + ggplot2 :: scale_y_continuous(labels = function(x) x * 100)
    
    p = p + ggnewscale::new_scale_colour()
    p = p + ggplot2 :: geom_point(data = d2, ggplot2 :: aes(x = series, y = relative), shape = 21, fill = "white", color = "black", position = ggplot2 :: position_jitter(w = 0.4, h = 0), size = pchSize)
    
    p = p + ggplot2 :: xlab("")
    p = p + ggplot2 :: ylab(paste(fg$geneName[i], "(%)"))
    p = p + ggplot2 :: theme(text= ggplot2 :: element_text(family="sans"),
                             legend.position = "none",
                             axis.title.x = ggplot2 :: element_text(size =textSize),
                             axis.title.y = ggplot2 :: element_text(size =textSize), 
                             panel.background = ggplot2 :: element_rect(fill = NA), 
                             panel.grid.major.y = ggplot2 :: element_line(colour = "grey90"),
                             panel.grid.minor.y = ggplot2 :: element_line(colour = "grey90"),
                             panel.grid.major.x = ggplot2 :: element_blank(),
                             axis.text.x = ggplot2 :: element_text(angle = 60, size=textSize, hjust = 1), 
                             axis.text.y = ggplot2 :: element_text(size=textSize))#,
    #p
    
    ggplot2 :: ggsave(filename = paste0(folderStoredPlots, "cumulative_", paste(fg$geneName[i],fg$chr[i],fg$start[i],fg$end[i],"2", formatPlot ,sep="_",collapse = "")), plot = p, width = sizeInchesPlots[1], height = sizeInchesPlots[2],units = "in")
  }
  
  return(invisible(NULL))  
}
