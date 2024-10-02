MiamiOvertargetPlot = function(fileResult, 
                               filePlot = "~/MiamiPlot",
                               formatPlot = c("png", "svg", "pdf")[1],
                               cloneFitness = FALSE,
                               nPrintedGenes = 20L,
                               highlightRobust = FALSE,
                               robust = FALSE,
                               keepRobustStat = TRUE,
                               removeRobustDifferentSign = TRUE,
                               typePlot = c("1", "2")[1], 
                               #typeCell = c("MSC", "adipMSC", "CD34LV")[1],
                               thresholdPvalue = 0.05,
                               annotateByFdrPvalue = FALSE, 
                               cancerGeneNameList = NULL,
                               sizePlot = c(12, 6),
                               colors = c("#090909", "#7c7c7c")) {
  #stopifnot(length(typeCell) == 1L)
  #stopifnot(typeCell %in% c("MSC", "adipMSC", "CD34LV"))
  stopifnot(length(typePlot) == 1L)
  stopifnot(typePlot %in% c("1", "2"))
  
  # MSC c("#6c9cff", "#0054ff")
  # adipMSC c("#72ba78", "#008b0b")
  # CD34LV c("#7c7c7c","#090909")
  
  #colorPoints = function(typeCell) {
  #  if(typeCell == "MSC") return(c("#6c9cff", "#0054ff"))
  #  if(typeCell == "adipMSC") return(c("#72ba78", "#008b0b"))
  #  if(typeCell == "CD34LV") return(c("#7c7c7c","#090909"))
  #}
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)
  #dataScore = data.table::fread(fileResult, sep=",", header = T)
  #cellT = typeCell
  
  if(highlightRobust) stopifnot("robustAdjPvalue" %in% names(dataScore))
  
  
  #print(summary(dataScore))
  
  dataScoreU = unique(dataScore,by = "geneName")
  dataScoreU = data.table::setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex = 1:nrow(dataScoreU)
  
  #dataScoreU$signTestStat=ifelse(dataScoreU$Zscore>0, dataScoreU$Chiscore,-1*dataScoreU$Chiscore)
  #dataScoreU$FDRpvalue=p.adjust(dataScoreU$pvalue,method = "fdr")
  #dataScoreU=setorder(dataScoreU,-Chiscore)
  #dataScoreU$annotate=ifelse((dataScoreU$FDRpvalue<0.05 & dataScoreU$signTestStat>0) ,1,0)
  
  #if(!cloneFitness) dataScoreU$signTestStat = ifelse(dataScoreU$Zscore>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  #else dataScoreU$signTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  
  if(robust) {
    if(keepRobustStat) whichNonRobust = which(dataScoreU$robustAdjPvalue > thresholdPvalue)
    else whichNonRobust = 1:nrow(dataScoreU)
    if(length(whichNonRobust) > 0L) {
      #dataScoreU$interaction [whichNonRobust] = dataScoreU$robustInteraction[whichNonRobust]
      dataScoreU$testStat    [whichNonRobust] = dataScoreU$robustTestStat   [whichNonRobust]
      dataScoreU$adjpvalue   [whichNonRobust] = dataScoreU$robustAdjPvalue  [whichNonRobust]
      if(cloneFitness) dataScoreU$signedTestStat[whichNonRobust] = dataScoreU$robustSignedTestStat[whichNonRobust]
    }
    if(removeRobustDifferentSign) {
      sameSign = sign(dataScoreU$interaction) == sign(dataScoreU$robustInteraction)
      dataScoreU$testStat = dataScoreU$testStat * as.double(sameSign)
      dataScoreU$signedTestStat = dataScoreU$signedTestStat * as.double(sameSign)
      dataScoreU$adjpvalue = ifelse(sameSign, dataScoreU$adjpvalue, 1.0)
    } 
  } 
  
  #if(cloneFitness) dataScoreU$signedTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  
  dataScoreU = data.table::setorder(dataScoreU, - testStat)
  #dataScoreU = data.table::setorder(dataScoreU, - signedTestStat)
  if(annotateByFdrPvalue) {
    dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
    dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue & dataScoreU$signedTestStat > 0, yes = 1, no = 0)
  }
  else {
    dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue & dataScoreU$signedTestStat > 0, yes = 1, no = 0)
  }
  
  whichExclude = grep('\\.', dataScoreU$geneName)
  if(length(whichExclude) > 0L) dataScoreU = dataScoreU[-whichExclude,] 
  #dataScoreU$npt = 0
  #dataScoreU$nIS = 0
  
  #dataScoreU$geneNameOnco=ifelse(dataScoreU$geneName %in% bushmanCancerGeneList$V1,paste(dataScoreU$geneName,"*",sep = ""),dataScoreU$geneName)
  if(!is.null(cancerGeneNameList)) cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  if(!is.null(cancerGeneNameList)) dataScoreU$geneNameOnco = ifelse(dataScoreU$geneName %in% cancerList, yes = paste(dataScoreU$geneName,"*",sep = ""), no = dataScoreU$geneName)
  else dataScoreU$geneNameOnco = dataScoreU$geneName
  
  dateYear = dataScoreU[annotate == 1]
  dateYear$labelling = 1
  if(nrow(dateYear) > 20) {
    dateYear$labelling = 0
    if(nPrintedGenes > 0) dateYear$labelling[1:nPrintedGenes] = 1
  }
  if(nPrintedGenes == 0) dateYear$labelling = 0
  
  
  if(highlightRobust) {
    highlight = which(dataScoreU$robustAdjPvalue <= thresholdPvalue & dataScoreU$adjpvalue <= thresholdPvalue)
    dateYear$labelRob = dateYear$labelling * (dateYear$robustAdjPvalue <= thresholdPvalue)
  } 
  
  textSize = 14
  textSizeRepel = 4
  
  p = ggplot2::ggplot(dataScoreU, ggplot2::aes(x=testIndex,y=signedTestStat,color = as.factor(chrInt)))
  p = p + ggplot2::theme(axis.title.y = ggplot2::element_text(size = textSize), 
                panel.background = ggplot2::element_rect(fill = NA),
                legend.position = "none",
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                axis.text.x = ggplot2::element_blank(), 
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = textSize))
  #p = p + geom_hline(yintercept = -log2(0.05), color = "grey40", linetype = "dashed")
  p = p + ggplot2::geom_point(alpha = 0.75,size=1)
  p = p + ggplot2::scale_color_manual(values = rep(colors, times = unique(length(dataScoreU$chr))))
  if(typePlot == "1") p = p + ggplot2::geom_point(data=subset(dateYear, annotate == 1), color = "#f40000", size = 3, shape = 1) 
  if(highlightRobust & (!robust)) p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(label=geneNameOnco), size = textSizeRepel, color=c("grey20", "#f40000")[1+subset(dateYear, labelling==1)$labelRob])
  else p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(label=geneNameOnco), size = textSizeRepel, color=c("grey20"))
  p = p + ggplot2::geom_hline(yintercept = 0,color="grey20",linewidth=0.3)
  if(typePlot == "2") p = p + ggplot2::geom_hline(yintercept = c(min(dateYear$testStat),-min(dateYear$testStat))[1:(2-cloneFitness)], color="#f40000",linetype = "dashed")
  # ylim(c(-max(dateYear$testStat),max(dateYear$testStat)))
  p = p + ggplot2::xlab("Chromosomes")
  p = p + ggplot2::ylab(c("Gene targeting score", "Gene/growth association score")[1 + cloneFitness])
  if(highlightRobust) p = p + ggplot2::annotate("point", x = dataScoreU$testIndex[highlight], y = dataScoreU$signedTestStat[highlight], pch = 1, color = "#f40000", size = 2.5)
  #ggsave(filename =  paste(unlist(strsplit(fileNameOut,"\\."))[1],"MiamiPlot2.png",sep="",collapse = ""), width = 12,height = 6,dpi = 300)
  ggplot2::ggsave(filename = paste0(filePlot, ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)
  
  return(invisible(NULL))
}

MiamiDifferentialPlot = function(fileResult, 
                                 filePlot = "~/MiamiDifferentialPlot",
                                 formatPlot = c("png", "svg", "pdf")[1],
                                 cloneFitness = FALSE, 
                                 #colPlot = c(1:3)[1], 
                                 thresholdPvalue = 0.05,
                                 annotateByFdrPvalue = FALSE,
                                 cancerGeneNameList = NULL,
                                 nPrintedGenes = 20,
                                 highlightRobust = FALSE,
                                 robust = FALSE,
                                 keepRobustStat = TRUE,
                                 removeRobustDifferentSign = TRUE,
                                 sizePlot = c(12, 6),
                                 colors1 = c("#0054ff", "#6c9cff"),
                                 colors2 = c("#090909", "#7c7c7c")) {
  
  #colorPoints = function(i, j) {
  #  stopifnot(length(i) == 1L & length(j) == 1L)
  #  stopifnot(i %in% 1:3 & j %in% c("col1L", "col1H", "col2L", "col2H"))
  #  cdt = data.table::data.table(col1L = c("#6c9cff", "#6c9cff", "#72ba78"),
  #                               col1H = c("#0054ff", "#0054ff", "#008b0b"),
  #                               col2L = c("#7c7c7c", "#72ba78", "#7c7c7c"),
  #                               col2H = c("#090909", "#008b0b", "#090909"))
  #  return(as.character(cdt[i, ..j]))
  #}
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)
  #dataScore = data.table::fread(fileResult, sep=",", header = T)
  
  if(highlightRobust) stopifnot("robustAdjPvalue" %in% names(dataScore))
  
  dataScoreU=unique(dataScore,by = "geneName")
  dataScoreU=data.table::setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex=1:nrow(dataScoreU)
  
  if(robust) {
    if(keepRobustStat) whichNonRobust = which(dataScoreU$robustAdjPvalue > thresholdPvalue)
    else whichNonRobust = 1:nrow(dataScoreU)
    #dataScoreU$interaction[whichNonRobust] = dataScoreU$robustInteraction[whichNonRobust]
    if(length(whichNonRobust) > 0L) {
      dataScoreU$testStat   [whichNonRobust] = dataScoreU$robustTestStat   [whichNonRobust]
      dataScoreU$adjpvalue  [whichNonRobust] = dataScoreU$robustAdjPvalue  [whichNonRobust]
      if(cloneFitness) dataScoreU$signedTestStat[whichNonRobust] = dataScoreU$robustSignedTestStat[whichNonRobust]
    }
    if(removeRobustDifferentSign) {
      sameSign = sign(dataScoreU$interaction) == sign(dataScoreU$robustInteraction)
      dataScoreU$testStat = dataScoreU$testStat * as.double(sameSign)
      dataScoreU$signedTestStat = dataScoreU$signedTestStat * as.double(sameSign)
      dataScoreU$adjpvalue = ifelse(sameSign, dataScoreU$adjpvalue, 1.0)
    } 
  } 
  
  #if(cloneFitness) dataScoreU$signedTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  #if(!cloneFitness) dataScoreU$signTestStat = ifelse(dataScoreU$Zscore>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  #else dataScoreU$signTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  
  dataScoreU = data.table::setorder(dataScoreU, - testStat)
  if(annotateByFdrPvalue) {
    dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
    dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  else {
    dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  
  whichExclude = grep('\\.', dataScoreU$geneName)
  if(length(whichExclude) > 0L) dataScoreU = dataScoreU[-whichExclude,]
  #dataScoreU$npt = 0
  #dataScoreU$nIS = 0
  
  if(!is.null(cancerGeneNameList)) cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  if(!is.null(cancerGeneNameList)) dataScoreU$geneNameOnco = ifelse(dataScoreU$geneName %in% cancerList, yes = paste(dataScoreU$geneName,"*",sep = ""), no = dataScoreU$geneName)
  else dataScoreU$geneNameOnco = dataScoreU$geneName
  
  dateYear = dataScoreU[dataScoreU$annotate == 1]
  dateYear$labelling = 1
  if(nrow(dateYear) > 20) {
    dateYear$labelling = 0
    if(nPrintedGenes > 0) dateYear$labelling[1:nPrintedGenes] = 1
  }
  if(nPrintedGenes == 0) dateYear$labelling = 0
  dateYearSwitch=ifelse(nrow(dateYear) > 0, yes = 1, no = 0)
  
  if(highlightRobust) {
    highlight = which(dataScoreU$robustAdjPvalue <= thresholdPvalue & dataScoreU$adjpvalue <= thresholdPvalue)
    dateYear$labelRob = dateYear$labelling * (dateYear$robustAdjPvalue <= thresholdPvalue)
  } 
  
  #textSize = c(14, 24)[1 + cloneFitness]
  #textSizeRepel = c(4, 6)[1 + cloneFitness]
  textSize = c(14)
  textSizeRepel = c(4)
  
  p = ggplot2::ggplot(dataScoreU[(chrInt %% 2L) == 0L,], ggplot2::aes(x = testIndex, y = signedTestStat, color = signedTestStat))
  p = p + ggplot2::theme(axis.title.y = ggplot2::element_text(size = textSize), 
                panel.background = ggplot2::element_rect(fill = NA),
                legend.position = "none",
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                axis.text.x = ggplot2::element_blank(), 
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = textSize))
  p = p + ggplot2::geom_point(alpha = 0.75,size=1)
  #p = p + ggplot2::scale_color_gradient2(low = colorPoints(colPlot, "col2H"), midpoint = 0, mid = "white", high = colorPoints(colPlot, "col1H"), limits = c(-0.5, 0.5), oob = scales :: squish)
  p = p + ggplot2::scale_color_gradient2(low = colors2[1], midpoint = 0, mid = "white", high = colors1[1], limits = c(-0.5, 0.5), oob = scales :: squish)
  p = p + ggnewscale::new_scale_colour()
  p = p + ggplot2::geom_point(data=dataScoreU[(chrInt %% 2L) == 1L,], ggplot2::aes(x = testIndex, y = signedTestStat, color = signedTestStat), alpha = 0.75, size = 1)
  #p = p + ggplot2::scale_color_gradient2(low = colorPoints(colPlot, "col2L"), midpoint = 0, mid = "white", high = colorPoints(colPlot, "col1L"), limits = c(-0.1, 0.1), oob = scales :: squish)
  p = p + ggplot2::scale_color_gradient2(low = colors2[2], midpoint = 0, mid = "white", high = colors1[2], limits = c(-0.1, 0.1), oob = scales :: squish)
  if(highlightRobust & (!robust)) p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(label=geneNameOnco), size = textSizeRepel, color=c("grey20", "#f40000")[1+subset(dateYear, labelling==1)$labelRob])
  else p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(label=geneNameOnco), size = textSizeRepel, color=c("grey20"))
  p = p + ggplot2::geom_hline(yintercept = 0, color="grey20", linewidth = 0.3)
  if(dateYearSwitch) p = p + ggplot2::geom_hline(yintercept = c(min(dateYear$testStat), - min(dateYear$testStat)), color= "#f40000",linetype = "dashed")
  p = p + ggplot2::ylim(c(-max(dataScoreU$testStat),max(dataScoreU$testStat)))
  p = p + ggplot2::xlab("Chromosomes")
  p = p + ggplot2::ylab(c("Gene targeting score", "Gene/growth association score")[1 + cloneFitness])
  if(highlightRobust) p = p + ggplot2::annotate("point", x = dataScoreU$testIndex[highlight], y = dataScoreU$signedTestStat[highlight], pch = 1, color = "#f40000", size = 2.5)
  #ggsave(filename =  paste(unlist(strsplit(fileNameOut,"\\."))[1],"DifferentialMiami.png",sep="",collapse = ""), width = 12,height = 6,dpi = 300)
  ggplot2::ggsave(filename = paste0(filePlot, ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)

  return(invisible(NULL))
}

BushmanDifferentialPlot = function(fileResult, 
                                   filePlot = "~/BushmanDifferentialPlot.png", 
                                   cloneFitness = FALSE, 
                                   colPlot = c(1:3)[1],
                                   thresholdPvalue = 0.05,
                                   annotateByFdrPvalue = FALSE, 
                                   cancerGeneNameList = NULL,
                                   maxOverlapPlot = 10) {
  colorPoints = function(i, ncol) {
    if(i == 1) colfunc = colorRampPalette(c("#090909", "white","#0054ff"))  
    if(i == 2) colfunc = colorRampPalette(c("#008b0b", "white","#0054ff"))
    if(i == 3) colfunc = colorRampPalette(c("#090909", "white","#008b0b"))  
    return(colfunc(ncol))
  }
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)
  #dataScore = data.table::fread(fileResult, sep=",", header = T)
  
  dataScoreU = unique(dataScore,by = "geneName")
  dataScoreU = data.table::setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex = 1:nrow(dataScoreU)
  
  if(cloneFitness) dataScoreU$signedTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  
  dataScoreU = data.table::setorder(dataScoreU, - testStat)
  if(annotateByFdrPvalue) {
    dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
    #dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue & dataScoreU$signedTestStat>0, yes = 1, no = 0)
    dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  else {
    #dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue & dataScoreU$signedTestStat>0, yes = 1, no = 0)
    dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  
  if(!is.null(cancerGeneNameList)) {
    cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
    dataScoreU = dataScoreU[dataScoreU$geneName %in% cancerList, ]
    dataScoreU = unique(dataScoreU, by = "geneName")
    dataScoreU = data.table::setorder(dataScoreU, - signedTestStat)
    dataScoreU$testIndex = 1:nrow(dataScoreU)
  }
  
  dateYear = dataScoreU[dataScoreU$annotate == 1]
  dateYear$labelling = 1
  if(nrow(dateYear) > 30) {
    dateYear$labelling = 0
    dateYear$labelling[1:15] = 1
    dateYear$labelling[(nrow(dateYear) - 15L) : nrow(dateYear)] = 1
  }
  dateYearSwitch=ifelse(nrow(dateYear) > 0, yes = 1, no = 0)
  
  textSize = 24
  textSizeRepel = 6
  
  qsts = list(identity, unique)[[1 + cloneFitness]](quantile(dataScoreU$signedTestStat, probs = seq(0, 1, 0.01)))
  ncol = c(101, length(qsts))[1 + cloneFitness]
  #ylmp = ifelse(cloneFitness, yes = max(dataScoreU$testStat), no = max(dataScoreU$signedTestStat))
  ylmp = max(dataScoreU$testStat)
  
  p = ggplot2 :: ggplot(dataScoreU, ggplot2 :: aes(x = testIndex, y = signedTestStat, color = cut(signedTestStat, qsts)))
  p = p + ggplot2 :: theme(text = ggplot2 :: element_text(size = textSize),
                axis.title.y = ggplot2 :: element_text(size = textSize), 
                panel.background = ggplot2 :: element_rect(fill = NA),
                legend.position = "none",
                panel.grid.major.x = ggplot2 :: element_blank(),
                panel.grid.minor.x = ggplot2 :: element_blank(),
                panel.grid.major.y = ggplot2 :: element_line(colour = "grey90"),
                panel.grid.minor.y = ggplot2 :: element_line(colour = "grey90"),
                axis.text.x = ggplot2 :: element_blank(), 
                axis.ticks.x = ggplot2 :: element_blank(),
                axis.text.y = ggplot2 :: element_text(size = textSize))
  p = p + ggplot2 :: geom_point(alpha = 0.75,size=1)
  p = p + ggplot2 :: scale_colour_manual(values = colorPoints(colPlot, ncol))
  p = p + ggrepel::geom_text_repel(data=subset(dateYear, labelling == 1), ggplot2 :: aes(label=geneName), size = textSizeRepel, color = "grey20",
                                   max.overlaps = getOption("ggrepel.max.overlaps", default = maxOverlapPlot))
  p = p + ggplot2 :: geom_hline(yintercept = 0, color = "grey20", linewidth = 0.3)
  if(dateYearSwitch) p = p + ggplot2 :: geom_hline(yintercept = c(min(dateYear$testStat), - min(dateYear$testStat)), color= "#f40000", linetype = "dashed")
  p = p + ggplot2 :: geom_hline(yintercept = mean(dataScoreU$signedTestStat), color = "#000000", linetype = "dotted", size = 1)
  p = p + ggplot2 :: ylim(c(-ylmp, ylmp))
  p = p + ggplot2 :: xlab("Genes") + ggplot2 :: ylab(c("Gene targeting score", "Gene/growth association score")[1 + cloneFitness])
  ggplot2 :: ggsave(filename = filePlot, width = 9, height = 9, dpi = 300)
  
  return(invisible(NULL))
}

