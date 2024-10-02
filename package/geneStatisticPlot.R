setGeneTable = function(fileGeneList, removeDuplicates = TRUE, removeDots = TRUE, alphabeticalOrder = FALSE) {
  geneTable = as.data.table(read.table(fileGeneList, header = TRUE))
  if(removeDuplicates) geneTable = unique(geneTable, by = "geneName")
  if(removeDots) {
    whichExclude = grep('\\.', geneTable$geneName)
    if(length(whichExclude) > 0L) geneTable = geneTable[-whichExclude,] 
    rm(whichExclude)
  }
  if(alphabeticalOrder) geneTable = geneTable[order(geneTable$geneName)]
  geneTable
}

setGenomeInfo = function() {
  length  = as.integer(c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569))
  start   = c(0, cumsum(as.numeric(head(length, -1))))
  bounds  = c(0, cumsum(as.numeric(length)))
  return(list(lengthChr = length, startChr = start, boundsChr = bounds))
}

genesStatisticPlot = function(fileGeneList, fileResults, idSamples, plotHighest = Inf, thresholdPvalue = 1.0, namePlottedGenes = character(), nameSelectedGenes = character(), typePlot = c("segments", "dots")[1], lineWidth = 0.5, transparency = 1, maxLineLength = Inf, inBasePairs = FALSE) {
  geneTable = setGeneTable(fileGeneList)
  if(missing(idSamples)) idSamples = as.character(1:length(fileResults))
  
  startChr = setGenomeInfo()$startChr
  
  tableResults = geneTable[,c("geneName", "chr_int")]
  tableResults$genePos = trunc(geneTable$start + (geneTable$end - geneTable$start) / 2)
  tableResults$genePos = purrr::map2_vec(tableResults$chr_int, tableResults$genePos, ~ startChr[.x] + .y)
  
  tableResults$gene_int = 1:nrow(tableResults)
  i=1
  idZscores = c()
  for(i in seq_along(fileResults)) {
    res_i = data.table::as.data.table(read.table(fileResults[i], header = TRUE))
    #`stopifnot(all(sapply(res_i$geneName, function(x) x %in% geneList)))
    res_i = res_i[, c("geneName", "Zscore", "adjpvalue")]
    res_i = res_i[1 : min(plotHighest, nrow(res_i)),]
    res_i$Zscore[which(res_i$adjpvalue >= thresholdPvalue)] = 0.0
    res_i$Zscore[which(res_i$Zscore < 0)] = 0.0
    
    new_col_name <- paste0("Z", i)
    tableResults <- merge(tableResults, res_i[, .(geneName, Zscore)], by = "geneName", all.x = TRUE)
    tableResults[, (new_col_name) := data.table::fcoalesce(Zscore, 0)]
    tableResults[, Zscore := NULL]
    idZscores = c(idZscores, new_col_name)
  }
  if(length(nameSelectedGenes) == 0) tableResults = tableResults[which(apply(tableResults[,..idZscores], 1, function(x) any(x > 0))),]
  else tableResults = tableResults[geneName %in% nameSelectedGenes,]
  tableResults = tableResults[order(gene_int)]
  tableResults$gene_int = 1:nrow(tableResults)
  maxZscore = ceiling(max(tableResults[, ..idZscores])) + 0
  nameSamples = idZscores
  
  #if(length(namePlottedGenes) > 0 & (!inBasePairs)) {
  if(length(namePlottedGenes) > 0) {  
    namePlottedGenes = unique(namePlottedGenes)
    namePlottedGenes = namePlottedGenes[which(sapply(namePlottedGenes, function(x) x %in% tableResults$geneName))]
    idPlottedGenes = sapply(namePlottedGenes, function(x) which(x == tableResults$geneName))
    if(!inBasePairs) plotted_gene_int = tableResults[idPlottedGenes]$gene_int
    else plotted_gene_int = tableResults[idPlottedGenes]$genePos
    chr_gene_int = tableResults[idPlottedGenes]$chr_int
    labelPlottedGenes = character(length(namePlottedGenes))
    for(j in seq_along(namePlottedGenes)) labelPlottedGenes[j] = paste(chr_gene_int[j], namePlottedGenes[j])
  }
  else { 
    labelPlottedGenes = character()
    plotted_gene_int = chr_gene_int = integer()
  }
  
  
  cat("number plotted genes :", nrow(tableResults),"\n")
  cat("number plotted non-zero in each dataset:\n")
  for(i in seq_along(idZscores)) cat(idSamples[i], ":", sum(tableResults[[idZscores[i]]] > 0), "\n")
  
  #tablePlot = data.table::data.table(geneName = character(), geneIndex = integer(), chrInt = integer(), zScore = numeric(), lower = numeric(), upper = numeric())
  tablePlot = data.table::data.table()
  i = 1L
  for(i in seq_along(fileResults)) {
    #chcol = colourChromosomes[[i]][1 + (tableResults$chr_int %% 2L)] 
    #chcol = as.character(2 * (i - 1) + 1:2)
    chcol = as.integer(2 * (i - 1) + 1:2)[1 + (tableResults$chr_int %% 2L)] 
    posgene = tableResults[[c("gene_int", "genePos")[1+inBasePairs]]]
    tablePlot = rbind(tablePlot, data.table::data.table(geneName = tableResults$geneName, 
                                                        geneIndex = posgene, 
                                                        chrInt = tableResults$chr_int,
                                                        nameSample = nameSamples[i],
                                                        zScore = tableResults[[idZscores[i]]],
                                                        lower = i * maxZscore - tableResults[[idZscores[i]]],
                                                        upper = sapply(tableResults[[idZscores[i]]], function(x) min(i * maxZscore , i * maxZscore - x + maxLineLength)),
                                                        chrColor = chcol))
  }
  tablePlot$chrColor = as.factor(tablePlot$chrColor)
  tablePlot$upHor = tablePlot$upper + tablePlot$upper - tablePlot$lower
  
  #x_major_breaks = unique(tablePlot$upper)
  x_major_breaks = (1 : length(fileResults)) * maxZscore
  #x_minor_breaks = as.numeric(sapply(x_major_breaks, function(x) x - c(6, 3)))
  x_minor_breaks = as.numeric(sapply(x_major_breaks, function(x) x + c(3, 6)))
  
  p = ggplot2::ggplot(data = tablePlot)
  if(typePlot == "segments") {
    #p = p + ggstance :: geom_linerangeh(mapping = ggplot2::aes(y = geneIndex, xmin = lower, xmax = upper, colour = chrColor), size = lineWidth, alpha = transparency)
    #p = p + scale_x_continuous(breaks = x_major_breaks, minor_breaks = x_minor_breaks, labels = idSamples)
    p = p + ggplot2 :: geom_linerange(mapping = ggplot2::aes(x = geneIndex, ymin = upper, ymax = upHor, colour = chrColor), size = lineWidth, alpha = transparency)
    p = p + scale_y_continuous(breaks = x_major_breaks, minor_breaks = x_minor_breaks, labels = idSamples)
  } 
  if(typePlot == "dots") {
    p = p + ggplot2 :: geom_point(mapping = ggplot2::aes(y = geneIndex, x = upper, size = zScore, colour = chrColor), alpha = transparency, shape = 1)
    p = p + scale_x_continuous(breaks = x_major_breaks, minor_breaks = NULL, labels = idSamples)
  } 
  #if(typePlot == "path") p = p + ggplot2 :: geom_path(mapping = ggplot2::aes(y = geneIndex, x = lower))
  p = p + scale_color_manual(values = unlist(colourChromosomes))
  
  #p = p + scale_y_continuous(name = NULL, breaks = plotted_gene_int, labels = labelPlottedGenes)
  p = p + scale_x_continuous(name = NULL, breaks = plotted_gene_int, labels = labelPlottedGenes)
  
  #if(!inBasePairs) p = p + scale_y_continuous(name = NULL, breaks = plotted_gene_int, labels = labelPlottedGenes)
  #else p = p + scale_y_continuous(name = NULL, minor_breaks = startChr[-1])
  p = p + ggplot2::theme(#axis.title.y = ggplot2::element_text(size = 7), 
                         panel.background = ggplot2::element_rect(fill = NA),
                         legend.position = "none",
                         panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
                         panel.grid.minor.x = ggplot2::element_line(colour = "grey90"),
                         #panel.grid.major.y = list(ggplot2::element_blank(), ggplot2::element_line(colour = "grey90"))[[1 +!inBasePairs]],
                         #panel.grid.minor.y = list(ggplot2::element_blank(), ggplot2::element_line(colour = "grey90"))[[1 + inBasePairs]],
                         panel.grid.major.y = list(ggplot2::element_blank(), ggplot2::element_line(colour = "grey90"))[[2]],
                         panel.grid.minor.y = list(ggplot2::element_blank(), ggplot2::element_line(colour = "grey90"))[[2]],
                         axis.text.x = ggplot2::element_text(size = 10, hjust = 1), 
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_text(size = 9))
  
  p
}


