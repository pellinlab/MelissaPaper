GoKeggPlots = function(fileResult, folderStoredPlots = "~", 
                       namePlots = c(ks = "keggScore.png", kh = "keggHeat.png", re = "reactome.png", gs = "goScore.png", gh = "goHeat.png"), 
                       cloneFitness = FALSE, annotateByFdrPvalue = FALSE, 
                       sizeInchesPlots = list(ks = c(9, 4), kh = c(9, 2), re = c(9, 3), gs = c(9,4), gh = c(9, 2))) {
  cat("select BioMart database\n")
  mart = biomaRt :: useMart('ENSEMBL_MART_ENSEMBL')
  cat("select BioMart dataset\n")
  mart = biomaRt :: useDataset('hsapiens_gene_ensembl', mart)
  
  cat("retrieve attributes from BioMart database\n")
  annotLookup = biomaRt :: getBM(
    mart = mart,
    attributes = c(
      'hgnc_symbol',
      'ensembl_gene_id',
      'gene_biotype',
      "ensembl_transcript_id",
      "entrezgene_id"),
    uniqueRows = TRUE)
  
  entrezRetrieval = function(genel){
    unlist(lapply(genel,function(xn) {
      xv=unlist(strsplit(xn, "\\."))[1]
      symb=unique(annotLookup$entrezgene_id[annotLookup$hgnc_symbol==xv])[1]
      if(is.na(symb)){
        symb=""
      }
      if(symb!=""){
        return(symb)
      }else{
        return("NA")
      }
    } ))
  }
  
  dataScore = data.table::fread(fileResult, sep="\t", header = T)
  dataScoreU = unique(dataScore,by = "geneName")
  whichExclude = grep('\\.', dataScoreU$geneName)
  if(length(whichExclude) > 0L) dataScoreU = dataScoreU[-whichExclude,]
  dataScoreU = data.table :: setorder(dataScoreU,chrInt,start)
  dataScoreU$testIndex = 1:nrow(dataScoreU)
  
  if(!cloneFitness) dataScoreU$signTestStat = ifelse(dataScoreU$Zscore>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  else dataScoreU$signTestStat = ifelse(dataScoreU$interaction>0, dataScoreU$testStat, -1 * dataScoreU$testStat)
  
  if(annotateByFdrPvalue) dataScoreU$FDRpvalue = p.adjust(dataScoreU$pvalue,method = "fdr")
  else dataScoreU$FDRpvalue = dataScoreU$adjpvalue
  dataScoreU$annotate = ifelse(dataScoreU$FDRpvalue < 0.05 & dataScoreU$signTestStat > 0, yes = 1, no = 0)
  
  dataScoreU$entrezid = entrezRetrieval(dataScoreU$geneName)
  
  dataScoreU=dataScoreU[!(is.na(dataScoreU$FDRpvalue)) & ((dataScoreU$entrezid!="NA")),]
  
  allResALL=NULL
  
  dataScoreU = dataScoreU[order(abs(dataScoreU$signTestStat), decreasing = T),]
  aG = data.table :: as.data.table(dataScoreU[,c("entrezid","signTestStat")])
  aGU = unique(aG,by="entrezid")
  aG = as.vector(aGU$signTestStat)
  names(aG) = as.character(aGU$entrezid)
  
  gene_list = aG
  gene_list = sort(gene_list, decreasing = TRUE)
  
  cat("compute gseKEGG\n")
  gse = clusterProfiler :: gseKEGG(geneList = gene_list,
                                   organism = 'hsa',
                                   #keyType = "ncbi-geneid", 
                                   minGSSize = 5, 
                                   maxGSSize = 800, 
                                   pvalueCutoff = 1, 
                                   pAdjustMethod = "fdr")
  
  gseaTable = data.table :: as.data.table(gse@result)
  gseaTable = data.table :: setorder(gseaTable, by = NES)
  gseaTable$sig_pos = paste(ifelse(gseaTable$p.adjust < 0.05, "Sig", "NoSig"), ifelse(gseaTable$NES < 0,"Neg", "Pos"), sep = "")
  
  topN=10
  
  if(length(sizeInchesPlots[[1]]) > 0L) {
    filePlot = paste0(folderStoredPlots, "/", namePlots[[1]])
    ks = ggplot2 :: ggplot(rbind(tail(gseaTable[NES>0,],topN),head(gseaTable[NES<0,],topN)), ggplot2 :: aes(reorder(Description, NES), NES)) 
    ks = ks + ggplot2 :: geom_col(ggplot2 :: aes(fill=sig_pos)) 
    ks = ks + ggplot2 :: scale_fill_manual("", 
                                           breaks = c("SigNeg", "NoSigNeg", "SigPos","NoSigPos"), 
                                           values = c("royalblue", "lightblue","red2", "coral1"), 
                                           labels = c("NES<0; FDR<0.05", "NES<0; FDR>0.05","NES>0; FDR<0.05", "NES>0; FDR>0.05"))
    ks = ks + ggplot2 :: coord_flip()
    ks = ks + ggplot2 :: labs(x = "KEGG term", y = "Normalized Enrichment Score", title = "")
    ks = ks + ggplot2 :: theme_minimal(base_size = 12)
    ggplot2 :: ggsave(filePlot, plot = ks, width = sizeInchesPlots[[1]][1], height = sizeInchesPlots[[1]][2], units = "in", dpi=300, bg = "white")
    #ggplot2 :: ggsave(paste(substring(fileNameOut,1,nchar(fileNameOut)-4), "__padj__KEGGGSEAPlot.png",sep = ""), plot=p1, width = 9 ,height = 4,units = "in",dpi=300,bg = "white")
  }
  if(length(sizeInchesPlots[[2]]) > 0L) {
    filePlot = paste0(folderStoredPlots, "/", namePlots[[2]])
    kh = enrichplot :: heatplot(gse, foldChange = gene_list, showCategory = 5, label_format = 30)
    ggplot2 :: ggsave(filePlot, plot = kh, width = sizeInchesPlots[[2]][1], height = sizeInchesPlots[[2]][2], units = "in", dpi=300, bg = "white")
    #ggsave(paste(substring(fileNameOut,1,nchar(fileNameOut)-4), "__padj__KEGGGSEAHeatPlot.png",sep = ""), plot=p2, width = 9 ,height = 2,units = "in",dpi=300,bg = "white")
  }
  # ?
  #fwrite(as.data.table(gse@result),file = "~/try.R,sep = "\t") 
  # ?
  
  cat("compute enrichments pathway\n")
  gse = ReactomePA :: enrichPathway(names(gene_list),
                                    organism = "human",
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH",
                                    minGSSize = 5,
                                    maxGSSize = 800,
                                    qvalueCutoff = 0.2)
  
  gseaTable = data.table :: as.data.table(gse@result)
  gseaTable = data.table :: setorder(gseaTable, by = p.adjust)
  #gseaTable=setorder(gseaTable,by=GeneRatio )
  
  gseaTable$sig_pos = ifelse(gseaTable$p.adjust < 0.05, "Sig", "NoSig")
  gseaTable$log2adjpvalue = -1 * log2(gseaTable$p.adjust)
  topN=10
  
  if(length(sizeInchesPlots[[3]]) > 0L) {
    filePlot = paste0(folderStoredPlots, "/", namePlots[[3]])
    re = ggplot2 :: ggplot(head(gseaTable, topN), ggplot2 :: aes(reorder(Description, log2adjpvalue), log2adjpvalue))
    re = re + ggplot2 :: geom_col(ggplot2 :: aes(fill = sig_pos))
    re = re + ggplot2 :: scale_fill_manual("", breaks = c("Sig", "NoSig"), values = c("red2", "coral1"), labels = c("FDR<0.05", "FDR>0.05"))
    re = re + ggplot2 :: coord_flip() 
    re = re + ggplot2 :: labs(x="Reactome Pathway", y="-Log2(Pvalue)", title="") 
    re = re + ggplot2 :: theme_minimal(base_size = 12)
    ggplot2 :: ggsave(filePlot, plot = re, width = sizeInchesPlots[[3]][1], height = sizeInchesPlots[[3]][2], units = "in", dpi = 300, bg = "white")
    #ggsave(paste(substring(fileNameOut,1,nchar(fileNameOut)-4), "__padj__ReactomeEnrichPlot.png",sep = ""), plot=p2, width = 9 ,height = 3,units = "in",dpi=300,bg = "white")
  }
  
  cat("compute gseGO\n")
  gse = clusterProfiler :: gseGO(geneList=gene_list, 
                                  ont ="BP", 
                                  keyType = "ENTREZID", 
               minGSSize = 5, 
               maxGSSize = 800, 
               pvalueCutoff = 1, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH")
  
  gseaTable = data.table :: as.data.table(gse@result)
  gseaTable = data.table :: setorder(gseaTable, by = NES)
  gseaTable$sig_pos = paste(ifelse(gseaTable$p.adjust<0.05, "Sig", "NoSig"), ifelse(gseaTable$NES<0, "Neg", "Pos"),sep = "")
  topN=10
  
  if(length(sizeInchesPlots[[4]]) > 0L) {
    filePlot = paste0(folderStoredPlots, "/", namePlots[[4]])
    gs = ggplot2 :: ggplot(rbind(tail(gseaTable[NES>0,], topN), head(gseaTable[NES<0,], topN)), ggplot2 :: aes(reorder(Description, NES), NES))
    gs = gs + ggplot2 :: geom_col(ggplot2 :: aes(fill=sig_pos)) 
    gs = gs + ggplot2 :: scale_fill_manual("", 
                                           breaks = c("SigNeg", "NoSigNeg", "SigPos","NoSigPos"),
                                           values = c("royalblue", "lightblue","red2", "coral1"), 
                                           labels = c("NES<0; FDR<0.05", "NES<0; FDR>0.05", "NES>0; FDR<0.05", "NES>0; FDR>0.05"))
    gs = gs + ggplot2 :: coord_flip() 
    gs = gs + ggplot2 :: labs(x="BP GO term", y="Normalized Enrichment Score",title="") 
    gs = gs + ggplot2 :: theme_minimal(base_size = 12)
    ggplot2 :: ggsave(filePlot, plot = gs, width = sizeInchesPlots[[4]][1], height = sizeInchesPlots[[4]][2], units = "in", dpi=300, bg = "white")
  }
  if(length(sizeInchesPlots[[5]]) > 0L) {
    filePlot = paste0(folderStoredPlots, "/", namePlots[[5]])
    gh = enrichplot :: heatplot(gse, foldChange = gene_list, showCategory = 5, label_format = 1)
    ggplot2 :: ggsave(filePlot, plot = gh, width = sizeInchesPlots[[5]][1], height = sizeInchesPlots[[5]][2], units = "in", dpi=300, bg = "white")
  }
  
  return(invisible(NULL))
}

