epigeneticPlot = function(integrSiteObject,
                          inputFiles,
                          chainFile,
                          lengthWindow,
                          filePlot = "~/epigenetic.png", sizeInchesPlots = c(5, 1)) {
  
  iSO = integrSiteObject
  it = iSO$design$filesIntegrations
  cv = iSO$design$covariatesIntegrations
  if(length(it) == 0L | length(cv) == 0L) stop("error")
  
  iSO = setGroups(iSO, group1 = 1:length(it), group2 = integer(0L))
  iSO = loadIntegrations(iSO, excludeIntegrations = list())
  tab = iSO$data$integrations
  groupsList = unique(cv$type)
  
  bbreak = ceiling(seq(- lengthWindow / 2, lengthWindow / 2, length.out = 300))
  #bbreak = (seq(- lengthWindow / 2, lengthWindow / 2, length.out = 30))
  chain = rtracklayer :: import.chain(chainFile)
  
  
  distDTall3=NULL
  for(i in seq_along(inputFiles)){
    distDTall=NULL
    
    cat("reading file", i, "/", length(inputFiles), "\n")
    bedMeth = data.table :: fread(inputFiles[i])
    bedMeth = bedMeth[, c(1:3)]
    colnames(bedMeth) = c("seqnames", "start", "end")
    
    cat("compute genomic ranges\n")
    bed_gr = GenomicRanges :: GRanges(seqnames = bedMeth$seqnames, ranges = IRanges :: IRanges(start = bedMeth$start, end = bedMeth$end))
    hg38_bed = rtracklayer :: liftOver(bed_gr, chain)
    hg38_bed = GenomicRanges :: GRanges(hg38_bed@unlistData)
    
    bedMethmid = hg38_bed@ranges@start + (hg38_bed@ranges@width / 2)
    
    for(j in 1:length(groupsList)){
      ISdata = do.call("rbind", tab[which(cv$type == groupsList[j])])
      ISdata = unique(ISdata, by = c("chr", "start"))
      genomic_interval = GenomicRanges :: GRanges(seqnames = ISdata$chr, ranges = IRanges :: IRanges(start = ISdata$start - (lengthWindow / 2), end = ISdata$end + (lengthWindow / 2)))
      overlapping_rows = IRanges :: findOverlaps(hg38_bed, genomic_interval)
      
      dists = bedMethmid[overlapping_rows@from] - ISdata$start[overlapping_rows@to]
      dists = dists[dists > -(lengthWindow / 2) & dists < (lengthWindow / 2)]
      hh = hist(dists, breaks = bbreak, plot = F)
      
      distDT = data.table :: data.table(breaks = hh$mids, den = hh$density)
      if(is.null(distDTall)) distDTall = distDT
      else distDTall = cbind(distDTall, distDT$den)
      colnames(distDTall)[ncol(distDTall)] = paste(groupsList[j])
    }
    
    if(is.null(distDTall3)) distDTall3 = distDTall
    else distDTall3 = rbind(distDTall3, distDTall)
  }
  distDTall3 = as.data.frame(distDTall3)
  distDTall3 = distDTall3[,-1]
  distPlot = data.frame(distDTall3[1:299,1], distDTall3[300:598,2], distDTall3[599:897,3])
  colnames(distPlot) = colnames(distDTall3)
  
  png(filename = filePlot, width = sizeInchesPlots[1], height = sizeInchesPlots[2], units = "in", pointsize = 12, res = 300)
  #pheatmap :: pheatmap(t(distDTall3), color = colorRampPalette(c("#f3e51c", "#20918c", "#450d56"))(100), 
  pheatmap :: pheatmap(t(distPlot), color = colorRampPalette(c("#f3e51c", "#20918c", "#450d56"))(100), 
                       cluster_rows = F, show_rownames = T, cluster_cols = F, show_colnames = F, legend = F,
                       angle_col = 90, 
                       #gaps_col = c(299, 598)[seq_len(length(inputFiles)-1L)],
                       gaps_row = seq_len(length(inputFiles)-1L),
                       na_col = "grey80", border_color = "white")
  dev.off()
  return(invisible(NULL))
}
