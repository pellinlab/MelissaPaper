MiamiOvertargetPlot = function(fileResult, 
                               storePlot = F,
                               filePlot = "~/MiamiPlot",
                               formatPlot = c("png", "svg", "pdf")[1],
                               cloneFitness = FALSE,
                               robust = FALSE,
                               cancerGeneNameList = NULL,
                               bpsDomain = FALSE,
                               log10Scale = FALSE,
                               
                               nPrintedGenes = 20L,
                               highlightAnnotated = FALSE,
                               #highlightRobust = FALSE,
                               keepRobustStat = FALSE,
                               emptyRobust = T,
                               removeRobustDifferentSign = TRUE,
                               
                               thresholdPvalue = 0.05,
                               annotateByFdrPvalue = FALSE, 
                               keepAboveThreshold = 10,

                               colors = c("#090909", "#7c7c7c"),
                               scaleByPvalue = FALSE,
                               sizePlot = c(12, 6)) {
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)
  #if(highlightRobust) stopifnot("robustAdjPvalue" %in% names(dataScore))
  
  dataScoreU = unique(dataScore,by = "geneName")
  dataScoreU = data.table::setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex = NULL
  dataScoreU$testIndex = 1:nrow(dataScoreU)
  if(bpsDomain) {
    chrs_domain = cumsum(c(0,248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
                           138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
                           83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569))
    dataScoreU$bpIndex = purrr::map2_dbl(dataScoreU$chrInt, dataScoreU$start, ~ chrs_domain[.x] + .y)
  } else {
    chrs_domain = c(1L, unique(dataScoreU$testIndex[diff(dataScoreU$chrInt) > 0L], nrow(dataScoreU)))
    dataScoreU$bpIndex = dataScoreU$testIndex
  } 
  
  
  if(robust) {
    if(keepRobustStat) whichNonRobust = which(dataScoreU$robustAdjPvalue > thresholdPvalue)
    else whichNonRobust = 1:nrow(dataScoreU)
    if(length(whichNonRobust) > 0L) {
      dataScoreU$testStat [whichNonRobust] = dataScoreU$robustTestStat [whichNonRobust]
      dataScoreU$adjpvalue[whichNonRobust] = dataScoreU$robustAdjPvalue[whichNonRobust]
      if(cloneFitness) dataScoreU$signedTestStat[whichNonRobust] = dataScoreU$robustSignedTestStat[whichNonRobust]
    }
    if(removeRobustDifferentSign) {
      sameSign = sign(dataScoreU$interaction) == sign(dataScoreU$robustInteraction)
      dataScoreU$testStat = dataScoreU$testStat * as.double(sameSign)
      dataScoreU$signedTestStat = dataScoreU$signedTestStat * as.double(sameSign)
      dataScoreU$adjpvalue = ifelse(sameSign, dataScoreU$adjpvalue, 1.0)
    } 
  } else {
    if(cloneFitness & emptyRobust) {
      dataScoreU$robustGenes = dataScoreU$robustAdjPvalue <= thresholdPvalue
      dataScoreU$robustGenes = dataScoreU$robustGenes & (sign(dataScoreU$robustInteraction) > 0)
      dataScoreU$robustGenes = dataScoreU$robustGenes & (sign(dataScoreU$interaction) > 0)
      # dataScoreU$robustGenes = dataScoreU$robustGenes & (sign(dataScoreU$interaction) == sign(dataScoreU$robustInteraction))
    }
  }
  
  # dataScoreU = data.table::setorder(dataScoreU, - testStat) # for differential target
  dataScoreU = data.table::setorder(dataScoreU, - signedTestStat)
  if(annotateByFdrPvalue) {
    dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
    dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue & dataScoreU$signedTestStat > 0, yes = T, no = F)
  } else {
    dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue & dataScoreU$signedTestStat > 0, yes = T, no = F)
  }
  
  if(cloneFitness & emptyRobust) dataScoreU$annotate = dataScoreU$annotate & dataScoreU$robustGenes
  
  whichExclude = grep('\\.', dataScoreU$geneName)
  if(length(whichExclude) > 0L) dataScoreU = dataScoreU[-whichExclude,] 
  
  if(!is.null(cancerGeneNameList)) cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  if(!is.null(cancerGeneNameList)) dataScoreU$geneNameOnco = ifelse(dataScoreU$geneName %in% cancerList, yes = paste(dataScoreU$geneName,"*",sep = ""), no = dataScoreU$geneName)
  else dataScoreU$geneNameOnco = dataScoreU$geneName
  
  dateYear = dataScoreU[annotate == T]
  data.table::setorder(dateYear, -signedTestStat)
  dateYear$labelling = 1
  if(nrow(dateYear) > 20) {
    dateYear$labelling = 0
    if(nPrintedGenes > 0) dateYear$labelling[1:min(nPrintedGenes, nrow(dateYear))] = 1
  }
  if(nPrintedGenes == 0) dateYear$labelling = 0
  
  #if(highlightRobust) {
  #  highlight = which(dataScoreU$robustAdjPvalue <= thresholdPvalue & dataScoreU$adjpvalue <= thresholdPvalue)
  #  dateYear$labelRob = dateYear$labelling * (dateYear$robustAdjPvalue <= thresholdPvalue)
  #} 
  
  if(scaleByPvalue) {
    #pv = -log10(dataScoreU$adjpvalue)
    #wt = purrr::map_int(pv, ~ which.max(.x > 5:0))
    #dataScoreU$alp = (5 / seq(from = 5, to = 13, length.out = 6))[wt]
    #print(head(dataScoreU$alp))
    #dataScoreU$alp[dataScore$annotate == T] = 1
    #dataScoreU$alp[dataScore$annotate == F] = .2
  } else { dataScoreU$alp = 1 }
  
  #print(dataScoreU)
  
  if(log10Scale) {
    dataScoreU = dataScoreU[dataScoreU$adjpvalue < keepAboveThreshold * thresholdPvalue & dataScoreU$signedTestStat > 0]
  }
  
  textSize = 14
  textSizeRepel = 4
  
  #p = ggplot2::ggplot(dataScoreU, ggplot2::aes(x=Zscore,y=log10adjpv))
  p = ggplot2::ggplot(dataScoreU)
  #if(scaleByPvalue) p = p + ggplot2::geom_point(ggplot2::aes(x = bpIndex, y = signedTestStat, colour = factor(2 - chrInt%%2), alpha = as.factor(annotate)), size = 1)
  if(cloneFitness & emptyRobust) {
    p = p + ggplot2::geom_point(ggplot2::aes(x = bpIndex, y = signedTestStat, colour = factor(2 - chrInt%%2), shape = robustGenes), size = 1)
    p = p + ggplot2::scale_shape_manual(values = c(1, 19))
  } else {
    p = p + ggplot2::geom_point(ggplot2::aes(x = bpIndex, y = signedTestStat, colour = factor(2 - chrInt%%2)), size = 1, shape = 19)
  }
  p = p + ggplot2::scale_color_manual(values = colors)
  if(scaleByPvalue) p = p + ggplot2::scale_alpha_manual(values = c(.1, 1))
  if(highlightAnnotated) p = p + ggplot2::geom_point(data=subset(dateYear, annotate == T), color = "#f40000", size = 3, shape = 1)
  #if(highlightRobust) p = p + ggplot2::annotate("point", x = dataScoreU$testIndex[highlight], y = dataScoreU$signedTestStat[highlight], pch = 1, color = "#f40000", size = 2.5)
  p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(x = bpIndex, y = signedTestStat, label=geneNameOnco), size = textSizeRepel, color=c("grey20"))
  if(!log10Scale) p = p + ggplot2::geom_hline(yintercept = 0,color="grey20",linewidth=0.3)
  #p = p + ggplot2::geom_hline(yintercept = c(min(dateYear$testStat),-min(dateYear$testStat))[1:(2-cloneFitness)], color="#f40000",linetype = "dashed")
  #p = p + ggplot2::geom_hline(yintercept = min(dateYear$testStat), color="#f40000",linetype = "dashed")
  p = p + ggplot2::geom_hline(yintercept = min(dateYear$testStat[dateYear$testStat > 0]), color="#f40000",linetype = "dashed")
  p = p + ggplot2::scale_x_continuous(breaks = chrs_domain, labels = NULL)
  if(log10Scale) {
    p = p + ggplot2::scale_y_log10()
    #p = p + ggplot2::scale_y_continuous(breaks = 10^(1:10), labels = 10^(1:10))
  } 
  #p = p + ggplot2::scale_x_continuous(breaks = chrs_domain, labels = NULL, minor_breaks = (head(chrs_domain, -1) + tail(chrs_domain, -1)) / 2, minor_labels = c(1:22, "X", "Y", ""))#c(1:22, "X", "Y", ""))
  p = p + ggplot2::xlab(paste0("Chromosomes (", c("ordinal gene coordinates", "base pair gene coordinates")[1+bpsDomain], ")"))
  p = p + ggplot2::ylab(paste(c("Gene targeting score", "Gene/growth association score")[1 + cloneFitness], c("", "(log10 scale)")[1 + log10Scale]))
  p = p + ggplot2::theme(axis.title.y = ggplot2::element_text(size = textSize),
                         axis.title.x = ggplot2::element_text(size = textSize),
                panel.background = ggplot2::element_rect(fill = NA),
                legend.position = "none",
                panel.grid.major.x = ggplot2::element_line(colour = c("white","grey90")[1+log10Scale]),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                axis.text.x = ggplot2::element_text(size = textSize), 
                #axis.ticks.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = textSize))
  #p
  if(storePlot) {
    ggplot2::ggsave(filename = paste0(filePlot, ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)
    return(invisible(NULL))
  } else return(p)
}

MiamiDifferentialPlot = function(fileResult, 
                                 storePlot = F,
                                 filePlot = "~/MiamiDifferentialPlot",
                                 formatPlot = c("png", "svg", "pdf")[1],
                                 cloneFitness = FALSE,
                                 robust = FALSE,
                                 cancerGeneNameList = NULL,
                                 bpsDomain = FALSE,
                                 log10Scale = FALSE,
                                 
                                 keepAboveThreshold = 20,
                                 nPrintedGenes = 20L,
                                 highlightAnnotated = FALSE,
                                 #highlightRobust = FALSE,
                                 keepRobustStat = FALSE,
                                 removeRobustDifferentSign = TRUE,
                                 emptyRobust = T,
                                 
                                 thresholdPvalue = 0.05,
                                 annotateByFdrPvalue = FALSE, 
                                 
                                 colors1 = c("#0054ff", "#6c9cff"),
                                 colors2 = c("#090909", "#7c7c7c"),
                                 scaleByPvalue = FALSE,
                                 sizePlot = c(12, 6)) {
  
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)

  dataScoreU = unique(dataScore,by = "geneName")
  dataScoreU = data.table::setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex = NULL
  dataScoreU$testIndex = 1:nrow(dataScoreU)
  if(bpsDomain) {
    chrs_domain = cumsum(c(0,248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
                           138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
                           83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569))
    dataScoreU$bpIndex = purrr::map2_dbl(dataScoreU$chrInt, dataScoreU$start, ~ chrs_domain[.x] + .y)
  } else {
    chrs_domain = c(1L, unique(dataScoreU$testIndex[diff(dataScoreU$chrInt) > 0L], nrow(dataScoreU)))
    dataScoreU$bpIndex = dataScoreU$testIndex
  } 
  
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
  } else {
    if(cloneFitness & emptyRobust) {
      dataScoreU$robustGenes = dataScoreU$robustAdjPvalue <= thresholdPvalue
      dataScoreU$robustGenes = dataScoreU$robustGenes & (sign(dataScoreU$robustInteraction) == sign(dataScoreU$interaction))
    }
  }
  
  dataScoreU = data.table::setorder(dataScoreU, - testStat)
  if(annotateByFdrPvalue) {
    dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
    dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  else {
    dataScoreU$annotate = ifelse(dataScoreU$adjpvalue < thresholdPvalue, yes = 1, no = 0)
  }
  
  if(cloneFitness & emptyRobust) dataScoreU$annotate = dataScoreU$annotate & dataScoreU$robustGenes
  
  whichExclude = grep('\\.', dataScoreU$geneName)
  if(length(whichExclude) > 0L) dataScoreU = dataScoreU[-whichExclude,]
  
  if(!is.null(cancerGeneNameList)) cancerList = data.table::fread(cancerGeneNameList, sep="\t", header = F)$V1
  if(!is.null(cancerGeneNameList)) dataScoreU$geneNameOnco = ifelse(dataScoreU$geneName %in% cancerList, yes = paste(dataScoreU$geneName,"*",sep = ""), no = dataScoreU$geneName)
  else dataScoreU$geneNameOnco = dataScoreU$geneName
  
  dateYear = dataScoreU[annotate == 1]
  dateYear$labelling = 1
  if(nrow(dateYear) > 20) {
    dateYear$labelling = 0
    if(nPrintedGenes > 0)
    {
      maxSTS = tail(sort(dateYear$signedTestStat), nPrintedGenes)
      minSTS = head(sort(dateYear$signedTestStat), nPrintedGenes)
      dateYear$labelling[purrr::map_lgl(dateYear$signedTestStat, ~ .x %in% c(minSTS,maxSTS))] = 1
    }
  }
  if(nPrintedGenes == 0) dateYear$labelling = 0
  dateYearSwitch=ifelse(nrow(dateYear) > 0, yes = 1, no = 0)
  
  #if(highlightRobust) {
  #  highlight = which(dataScoreU$robustAdjPvalue <= thresholdPvalue & dataScoreU$adjpvalue <= thresholdPvalue)
  #  dateYear$labelRob = dateYear$labelling * (dateYear$robustAdjPvalue <= thresholdPvalue)
  #} 
  
  # if(scaleByPvalue) {
  #   pv = -log10(dataScoreU$adjpvalue)
  #   wt = purrr::map_int(pv, ~ which.max(.x >= 5:0))
  #   #dataScoreU$alpha = 1/c(1, 1.5, 3, 5, 10, 20)[wt]
  #   dataScoreU$alpha = 1/c(1, 1, 1.5, 3, 5, 10)[wt]
  # } else dataScoreU$alpha = 1 #.75
  # 
  
  
  textSize = c(14)
  textSizeRepel = c(4)
  
  if(!log10Scale) {
    
    p = ggplot2::ggplot()
    #if(scaleByPvalue) p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 0L,], mapping = ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat, alpha = alpha), size = 1)
    #else p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 0L,], mapping = ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat), size = 1)
    if(cloneFitness & emptyRobust) {
      p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 0L,], ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat, shape = robustGenes), size = 1)
      p = p + ggplot2::scale_shape_manual(values = c(1, 19))
    } else {
      p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 0L,], ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat), size = 1, shape = 19)
    }
    p = p + ggplot2::scale_color_gradient2(low = colors2[1], midpoint = 0, mid = "white", high = colors1[1], limits = c(-0.5, 0.5), oob = scales :: squish)
    p = p + ggnewscale::new_scale_colour()
    #if(scaleByPvalue) p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 1L,], mapping = ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat, alpha = alpha), size = 1)
    #else p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 1L,], mapping = ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat), size = 1)
    if(cloneFitness & emptyRobust) {
      p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 1L,], ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat, shape = robustGenes), size = 1)
      p = p + ggplot2::scale_shape_manual(values = c(1, 19))
    } else {
      p = p + ggplot2::geom_point(data = dataScoreU[(chrInt %% 2L) == 1L,], ggplot2::aes(x = bpIndex, y = signedTestStat, colour = signedTestStat), size = 1, shape = 19)
    }
    p = p + ggplot2::scale_color_gradient2(low = colors2[2], midpoint = 0, mid = "white", high = colors1[2], limits = c(-0.5, 0.5), oob = scales :: squish)
    if(highlightAnnotated) p = p + ggplot2::geom_point(data=subset(dateYear, annotate == 1), color = "#f40000", size = 3, shape = 1)
    p = p + ggrepel::geom_text_repel( data=subset(dateYear, labelling==1), ggplot2::aes(x = bpIndex, y = signedTestStat, label=geneNameOnco), size = textSizeRepel, color=c("grey20"))
    p = p + ggplot2::geom_hline(yintercept = 0,color="grey20",linewidth=0.3)
    p = p + ggplot2::geom_hline(yintercept = c(-1,1) * min(dateYear$testStat[dateYear$testStat > 0]), color="#f40000",linetype = "dashed")
    p = p + ggplot2::scale_x_continuous(breaks = chrs_domain, labels = NULL)
    p = p + ggplot2::ylim(c(-max(dataScoreU$testStat),max(dataScoreU$testStat)))
    p = p + ggplot2::xlab("Chromosomes")
    p = p + ggplot2::ylab(paste(c("Gene targeting score", "Gene/growth association score")[1 + cloneFitness], c("", "(log10 scale)")[1 + log10Scale]))
    p = p + ggplot2::theme(axis.title.y = ggplot2::element_text(size = textSize),
                           axis.title.x = ggplot2::element_text(size = textSize),
                           panel.background = ggplot2::element_rect(fill = NA),
                           legend.position = "none",
                           panel.grid.major.x = ggplot2::element_line(colour = "white"),
                           panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                           panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
                           axis.text.x = ggplot2::element_text(size = textSize), 
                           #axis.ticks.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_text(size = textSize))
    
    if(storePlot) {
      ggplot2::ggsave(filename = paste0(filePlot, ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)
      return(invisible(NULL))
    } else return(p)
  }
  if(log10Scale) {
    dataScoreU = dataScoreU[dataScoreU$adjpvalue < keepAboveThreshold * thresholdPvalue]
    dataScoreU$chrEven = factor(dataScoreU$chrInt %% 2L == 0L)
    rangeY = c(1, max(abs(dataScoreU$signedTestStat)))
    slines = min(dataScoreU$testStat[dataScoreU$adjpvalue <= thresholdPvalue])
    
    d1 = dataScoreU[signedTestStat >=  1]
    d2 = dataScoreU[signedTestStat <= -1]
    d2$signedTestStat = -1 * d2$signedTestStat
    
    p1 = ggplot2::ggplot(data = d1, mapping = aes(x = bpIndex, y =  signedTestStat, colour = chrEven)) +
      geom_point(aes(alpha = factor(annotate))) + scale_color_manual(values = colors1) +
      scale_alpha_manual(values = list(c(1, 1), c(.1,1))[[1+scaleByPvalue]]) +
      geom_hline(yintercept = 1,  linewidth = .1, color = colors1[1]) +
      ggrepel::geom_text_repel(data=subset(d1, annotate==1), ggplot2::aes(x = bpIndex, y = signedTestStat, label=geneNameOnco), size = textSizeRepel, color=c("grey20"), direction = "both") +
      ggplot2::geom_hline(yintercept = slines, color="#f40000", linetype = "dashed") +
      scale_y_continuous(trans = c("log10")) +
      scale_x_continuous(breaks = chrs_domain, labels = NULL, minor_breaks = NULL) +
      coord_cartesian(xlim = range(chrs_domain), ylim = rangeY[1:2]) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        #axis.text.x = ggplot2::element_blank(), 
        axis.ticks.x = ggplot2::element_blank(),
        #axis.text.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = NA),
        panel.grid.major.x = ggplot2::element_line(colour = c("white","grey90")[2]),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor.y = ggplot2::element_line(colour = "grey90")
      )
    
    p2 = ggplot2::ggplot(data = d2, aes(x = bpIndex, y = signedTestStat, colour = chrEven)) +
      geom_point(aes(alpha = factor(annotate))) + scale_color_manual(values = colors2) +
      scale_alpha_manual(values = list(c(1, 1), c(.1,1))[[1+scaleByPvalue]]) +
      geom_hline(yintercept = 1, linewidth = .1, color = colors2[1]) +
      ggplot2::geom_hline(yintercept = slines, color="#f40000", linetype = "dashed") +
      #scale_y_continuous(trans = c("log10", "reverse")) +
      scale_y_continuous(trans = c("log10")) +
      ggrepel::geom_text_repel(data=subset(d2, annotate==1), ggplot2::aes(x = bpIndex, y = signedTestStat, label=geneNameOnco), size = textSizeRepel, color=c("grey20"), direction = "both") +
      scale_x_continuous(breaks = chrs_domain, labels = NULL, minor_breaks = NULL) +
      coord_cartesian(xlim = range(chrs_domain), ylim = rangeY[1:2]) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        #axis.text.x = ggplot2::element_blank(), 
        axis.ticks.x = ggplot2::element_blank(),
        #axis.text.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = NA),
        panel.grid.major.x = ggplot2::element_line(colour = c("white","grey90")[2]),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor.y = ggplot2::element_line(colour = "grey90")
      )
    
    if(storePlot) {
      ggplot2::ggsave(filename = paste0(filePlot, "_1", ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)
      ggplot2::ggsave(filename = paste0(filePlot, "_2", ".", formatPlot), width = sizePlot[1], height = sizePlot[2], dpi = 300)
      return(invisible(NULL))
    } else return(list(p1, p2))
  }
}

BushmanDifferentialPlot = function(fileResult, 
                                   filePlot = "~/BushmanDifferentialPlot.png", 
                                   cloneFitness = FALSE, 
                                   colPlot = c(1:3)[1],
                                   thresholdPvalue = 0.05,
                                   annotateByFdrPvalue = FALSE, 
                                   cancerGeneNameList = NULL,
                                   maxOverlapPlot = 10,
                                   sizePlot = c(6, 6)) {
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
  
  textSize = 14
  textSizeRepel = 4
  
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
  ggplot2 :: ggsave(filename = filePlot, width = sizePlot[1], height = sizePlot[2], dpi = 300)
  
  return(invisible(NULL))
}

