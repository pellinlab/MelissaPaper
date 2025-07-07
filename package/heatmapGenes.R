HeatmapGenes = function(imputedFiles, filePlot = "~/heatmap.png",
                        includeGenes = "",
                        cancerGeneNameList, typeCells,
                        numberTopGenes = 10L, cloneFitness = FALSE,
                        plotOnlyInList = FALSE) {
  
  cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  
  topGenes=c()
  for(i in imputedFiles){
    dataScore = data.table :: fread(i, sep="\t",header = T)
    
    whichExclude = grep('\\.', dataScore$geneName)
    if(length(whichExclude) > 0L) dataScore = dataScore[-whichExclude,]
    
    if(!cloneFitness) dataScore$signTestStat = ifelse(dataScore$Zscore>0, dataScore$testStat, -1 * dataScore$testStat)
    else dataScore$signTestStat = ifelse(dataScore$interaction>0, dataScore$testStat, -1 * dataScore$testStat)
    
    if(plotOnlyInList) dataScore = dataScore[dataScore$geneName %in% cancerList,]
    dataScore= data.table :: setorder(dataScore,-signTestStat)
    
    topGenes=c(topGenes,dataScore$geneName[1:numberTopGenes])
  }
  topGenes = unique(c(topGenes, includeGenes))
  
  topGenesDT = data.table :: data.table(geneName=topGenes)
  topGenesDT$index=1:nrow(topGenesDT)
  topGenesDT=unique(topGenesDT,by="geneName")
  
  heatPlot=topGenesDT
  for(i in 1:length(imputedFiles)){
    dataScore = data.table :: fread(imputedFiles[i], sep="\t",header = T)
    if(!cloneFitness) dataScore=dataScore[, c("geneName", "Zscore")]
    else dataScore=dataScore[, c("geneName", "interaction")]
    #heatPlot=merge(heatPlot,dataScore,all.x=T)
    heatPlot = merge(heatPlot,dataScore, by = "geneName", all.x=T)
  }  
  
  #cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  heatPlot$geneName=ifelse(heatPlot$geneName %in% cancerList, paste(heatPlot$geneName,"*",sep = ""), heatPlot$geneName)
  ########
   heatPlot = unique(heatPlot, by = "geneName")
  ########
  heatPlot = data.table :: setorder(heatPlot,index)
  
  
  heatPlot=as.data.frame(heatPlot)
  rownames(heatPlot)=heatPlot$geneName
  heatPlot=heatPlot[,-c(1:2)]
  #heatPlot[is.na(heatPlot)] <- 0
  #colnames(heatPlot)=c("BM MSC", "adip MSC", "HSPC")
  #colnames(heatPlot)=c("BM MSC", "adip MSC")
  colnames(heatPlot) = typeCells
  #rownames(heatPlot) = paste0("X", 1:657)
  
  #png(filename = paste(paste(imputedFiles,collapse = "_",sep="_"),".png",sep=""),    width =9, height = 1.5, units = "in", pointsize = 12,res = 600)
  png(filename = filePlot, width =9, height = 1.5, units = "in", pointsize = 12,res = 600)
  pheatmap :: pheatmap(t(heatPlot), color = colorRampPalette(c("white","yellow", "orange","firebrick3"))(100),cluster_rows=T, show_rownames=T,cluster_cols=F,angle_col=90,na_col="grey80",border_color="white")
  dev.off()

  return(invisible())
}