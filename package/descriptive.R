CloneSizesPlot = function(integrSiteObject, filePlot = "~/cloneSizes.png", sizeInchesPlots = c(6, 6)) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  if(length(iSO$data$combined) != 0L) stop("error") # not yet implemented
  
  #iSO = setGroups(iSO, group1 = 1:length(it), group2 = integer(0L))
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  
  tab = iSO$data$integrations
  len = tot = vector("integer", length(tab))
  for(i in seq_along(tab)) {
    len[[i]] = nrow(tab[[i]])
    tot[[i]] = sum(tab[[i]]$count)
  }
  
  colorCell = list(`BM MSC` = "#0054ff", `Ad MSC` = "#008b0b", `HSPC` = "#090909")
  textSize = 8
  pchSize = 2
  minAlpha = .2
  maxAlpha = .8
  minTime = min(cv$time)
  maxTime = max(cv$time)
  
  alp = sapply(cv$time, function(x) (x - minTime) / (maxTime - minTime))
  d = data.frame(series = 1:nrow(cv), length = len, total = tot, type = cv$type, labels = paste(cv$type, cv$idty), alpha = alp)
  
  p = ggplot2 :: ggplot(data = d)
  p = p + ggplot2 :: geom_col(ggplot2 :: aes(x = series, y = total, fill = type, alpha = alpha))
  p = p + ggplot2 :: scale_fill_manual(values = as.list(colorCell))
  p = p + ggplot2 :: scale_alpha_continuous(range = c(minAlpha, maxAlpha))
  p = p + ggplot2 :: scale_x_continuous(breaks = 1:nrow(cv), labels = d$labels)
  
  p = p + ggnewscale::new_scale_colour()
  p = p + ggplot2 :: geom_point(ggplot2 :: aes(x = series, y = length), shape = 21, fill = "white", color = "black", size = pchSize)
  
  p = p + ggplot2 :: xlab("")
  p = p + ggplot2 :: ylab("Clone size")
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
  
  ggplot2 :: ggsave(filename = filePlot, plot = p, width = sizeInchesPlots[1], height = sizeInchesPlots[2], units = "in")
  return(invisible(NULL))
}

DiversityIndexPlot = function(integrSiteObject, filePlot = "~/diversityIndex.png", plottedCovariate = "time", index = c("shannon", "simpson", "invsimpson")[1], sizeInchesPlots = c(6, 6)) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  stopifnot(length(index) == 1L)
  stopifnot(index %in% c("shannon", "simpson", "invsimpson"))
  legend.y = c(shannon = "H", simpson = "1-D", invsimpson = "1/D")[[index]]
  
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  stopifnot(plottedCovariate %in% names(cv))
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  if(length(iSO$data$combined) != 0L) stop("error") # not yet implemented
  
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  
  tab = iSO$data$integrations
  div = vector("double", length(tab))
  for(i in seq_along(tab)) div[[i]] = vegan :: diversity(tab[[i]]$count)
  
  colorCell = list(`Ad MSC` = "#72bf72", `BM MSC` = "#4267c6", `HSPC` = "#878787")
  #alpha = .5
  #alphaCell = list(`Ad MSC` = alpha, `BM MSC` = alpha, `HSPC` = alpha)
  textSize = 8
  
  d = data.frame(x = cv[[plottedCovariate]], y = div, type = cv$type)
  dm = dplyr::summarize(dplyr::group_by(d, x, type), y = mean(y), .groups = 'drop')
  di = split(dm, dm$type)
  
  p = ggplot2 :: ggplot(d, ggplot2 :: aes(x = factor(x), y = y, color = type, group = type))
  p = p + ggplot2 :: geom_point(size=3) #+ ggplot2 :: geom_line(linewidth=1)
  for(i in 1:length(di)) p = p + ggplot2 :: geom_path(data = di[[i]], ggplot2 :: aes(x = factor(x), y = y, color = type))
  p = p + ggplot2 :: scale_color_manual(values = unlist(colorCell), breaks=names(colorCell))
  #p = p + ggplot2 :: scale_alpha_manual(values = unlist(alphaCell), breaks=names(alphaCell))
  p = p + ggplot2 :: xlab(plottedCovariate) + ggplot2 :: ylab(legend.y)
  p = p + ggplot2 :: theme(text= ggplot2 :: element_text(family="sans"),
                           legend.position = "right",
                           axis.title.x = ggplot2 :: element_text(size =textSize),
                           axis.title.y = ggplot2 :: element_text(size =textSize), 
                           panel.background = ggplot2 :: element_rect(fill = NA), 
                           panel.grid.major.y = ggplot2 :: element_line(colour = "grey90"),
                           panel.grid.minor.y = ggplot2 :: element_line(colour = "grey90"),
                           panel.grid.major.x = ggplot2 :: element_line(colour = "grey90"),
                           axis.text.x = ggplot2 :: element_text(size=textSize), 
                           axis.text.y = ggplot2 :: element_text(size=textSize))

  ggplot2 :: ggsave(filename = filePlot, plot = p, width = sizeInchesPlots[1], height = sizeInchesPlots[2], units = "in")
  return(invisible(NULL))
}

AnnotationDistributionPlot = function(integrSiteObject, filePlot = "~/annotationDistribution.png", sizeInchesPlots = c(6, 6)) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  if(length(iSO$data$combined) != 0L) stop("error") # not yet implemented
  
  annotations = annotatr :: build_annotations(genome = 'hg38', annotations = c('hg38_basicgenes', 'hg38_genes_intergenic'))
  annFactor = data.table :: data.table(annot.type = unique(annotations$type), numAnn = c(5,6,4,1,2,3,7))
  
  iSO = setGroups(iSO, group1 = 1:length(it), group2 = integer(0L))
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  #iSO = loadLocationsGenes(iSO)
  #iSO = mergeData(iSO)
  
  tab = iSO$data$integrations
  #lcg = iSO$data$locationsGenes
  
  annDTAll=NULL
  for(i in seq_along(tab)) {
    gr = GenomicRanges :: GRanges(
      seqnames = tab[[i]]$chr,
      strand = tab[[i]]$strand,
      ranges = IRanges :: IRanges(start = tab[[i]]$start, end = tab[[i]]$end),
      #gene = tab[[i]]$geneName,
      #cellCount= vv$readsCount
      cellCount = nrow(tab[[i]])
    )
    
    dm_annotated = annotatr :: annotate_regions(
      regions = gr,
      annotations = annotations,
      ignore.strand = TRUE,
      quiet = FALSE)
    
    annDT = data.table :: as.data.table(dm_annotated)
    annDT = merge(annDT, annFactor, by = "annot.type", all.x = TRUE)
    annDT = data.table :: setorder(annDT, seqnames,start,numAnn)
    annDT = unique(annDT, by=c("seqnames","start"))
    annDT$fileName = it[[i]]
    annDT$series = i
    if(is.null(annDTAll)) annDTAll = annDT
    else annDTAll = rbind(annDTAll, annDT)
  }
  
  #annDTAll = merge(annDTAll, dataInfo, by="fileName")
  
  p = ggplot2 :: ggplot(data=annDTAll, ggplot2 :: aes(
    x = factor(series, labels = paste(cv$type, cv$idty)), 
    fill = factor(annot.type,levels = rev(c("hg38_genes_introns","hg38_genes_exons","hg38_genes_promoters", "hg38_genes_intergenic","hg38_genes_1to5kb")), labels = rev(c("intron","exon","promoters", "intergenic","1to5kb"))) 
    ))
  
  p = p + ggplot2 :: geom_bar(position = "fill")
  p = p + ggplot2 :: scale_fill_manual(values = wesanderson :: wes_palette("Darjeeling1",5),name = "Regions")
  p = p + ggplot2 :: xlab("") + ggplot2 :: ylab("Proportion")
  p = p + ggplot2 :: theme(
    panel.background = ggplot2 :: element_rect(fill = NA),
    panel.grid.major.x = ggplot2 :: element_blank(),
    panel.grid.minor.x = ggplot2 :: element_blank(),
    panel.grid.major.y = ggplot2 :: element_blank(),
    panel.grid.minor.y = ggplot2 :: element_blank(),
    axis.text.x = ggplot2 :: element_text(angle = 60, hjust = 1)
  )

  ggplot2 :: ggsave(filename = filePlot, plot = p, width = sizeInchesPlots[1], height = sizeInchesPlots[2], units = "in")
  return(invisible(NULL))
}

