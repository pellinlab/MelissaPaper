VennPlot = function(filesResult, nameCellTypes, filePlot = "~/Venn.svg", 
                    colorCellTypes = c("#6c9cff","#72ba78","#7c7c7c")[1:length(nameCellTypes)]) {
  stopifnot(length(filesResult) == length(nameCellTypes))
  listResults = vector("list", length(filesResult))
  for(i in seq_along(listResults)) {
    ls = readr :: read_tsv(filesResult[[i]])
    whichExclude = grep('\\.', ls$geneName)
    if(length(whichExclude) > 0L) ls = ls[-whichExclude,]
    ls$adjpvalue = p.adjust(ls$pvalue, method = "fdr")
    lowAdj = ls$adjpvalue < 0.05
    higZsc = ls$Zscore > 0
    listResults[[i]] = ls$geneName[lowAdj & higZsc]
    cat("\nfor", nameCellTypes[i], ",", length(listResults[[i]]), "/", length(ls$geneName), "are over targeted\n\n")
  }
  
  p = nVennR :: plotVenn(listResults, 
                         sNames = nameCellTypes,
                         outFile = filePlot,
                         setColors = colorCellTypes,
                         opacity = 0.5,
                         borderWidth = 2,
                         showLegend = T,
                         systemShow = FALSE,
                         labelRegions = F,
                         showNumbers = T,
                         fontScale = 3)
  return(invisible(NULL))
}

