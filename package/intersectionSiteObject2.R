createEmptyIntegrSiteObject = function() {
  integrSiteObject = list(
    organism = list(type = character(),
                    infoChromosomes = data.frame(),
                    genome = character(),
                    fileInfoGenes = character()),
    design = list(id = integer(),
                  filesIntegrations = character(),
                  covariatesIntegrations = data.frame(),
                  selectedCovariates = character(),
                  selCovAreFactors = logical()),
    data = list(groups = list(),
                integrations = list(),
                combined = data.table :: data.table(),
                locationsGenes = data.table :: data.table(),
                splitOverGenes = list()),
    model = list(twoGroups = logical(),
                 groupedByCovariate = logical(),
                 contrastFactors = list(),
                 typeRegressionModel = character(),
                 estimatedH0 = list(),
                 design = NULL,
                 response = NULL,
                 namesParameters = character(),
                 lengthInvestigatedGenome = integer(),
                 sumsOutcome = integer(),
                 lengthsOutcome = integer()),
    infoStorage = list(pathWorkingDirectory = "~/")
  )
  class(integrSiteObject) = "integrSiteObject"
  return(integrSiteObject)
}

setWorkingDirectory = function(integrSiteObject, pathWorkingDirectory = "~/") {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  integrSiteObject$infoStorage$pathWorkingDirectory = pathWorkingDirectory
  return(integrSiteObject)
}

includeOrganismHuman = function() {
  infoChromosomes = data.frame(
    idInteger = as.integer(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)),
    length = as.integer(c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569)),
    idCharacter = as.character(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")),
    idChar2digits = as.character(c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"))
  )
  rownames(infoChromosomes) = infoChromosomes$idCharacter
  genome = "hg38"
  return(list(infoChromosomes = infoChromosomes, genome = genome))
}

includeOrganismMouse = function() {
  infoChromosomes = data.frame()
  genome = character(0)
  return(list(infoChromosomes = infoChromosomes, genome = genome))
}

includeOrganismSimulation = function() {
  infoChromosomes = data.frame(idInteger = 15L, length = 101991189L, idCharacter = "chr15", idChar2digits = "chr15")
  genome = "hg38chr15"
  return(list(infoChromosomes = infoChromosomes, genome = genome))
}

setFilesIntegration = function(integrSiteObject, filesIntegrations, id) {
  integrSiteObject$design$id = id
  integrSiteObject$design$filesIntegrations = filesIntegrations
  names(integrSiteObject$design$filesIntegrations) = id
  return(integrSiteObject)
}

setCovariatesIntegrations = function(integrSiteObject, covariatesIntegrations, id) {
  integrSiteObject$design$covariatesIntegrations = covariatesIntegrations
  rownames(integrSiteObject$design$covariatesIntegrations) = id
  return(integrSiteObject)
}

setInfoSelectedCovariates = function(integrSiteObject, selectedCovariates = character(), selCovAreFactors = logical()) {
  integrSiteObject$design$selectedCovariates = selectedCovariates
  integrSiteObject$design$selCovAreFactors   = selCovAreFactors
  return(integrSiteObject)
}

setGroups = function(integrSiteObject, group1 = integer(), group2 = integer()) {
  iSO = integrSiteObject
  iSO$model$twoGroups = (length(group1) > 0L) & (length(group2) > 0L)
  iSO$data$groups = list(`1` = group1, `2` = group2)
  sc = iSO$design$selectedCovariates
  sf = iSO$design$selCovAreFactors
  gbc = FALSE
  if(length(sc)) if(any(sf)) for(i in which(sf)) {
    cov = iSO$design$covariatesIntegrations[,i, drop = T]
    n0 = length(unique(cov))
    n1 = length(unique(cov[group1]))
    n2 = length(unique(cov[group2]))
    if((n1 + n2) == n0) gbc = TRUE
  }
  iSO$model$groupedByCovariate = gbc
  return(iSO)
}

loadIntegrations = function(integrSiteObject, excludeIntegrations = list(), header = FALSE, variableNames = c("chr", "start",	"end", "strand", "count")) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  
  wd = integrSiteObject$infoStorage$pathWorkingDirectory
  fi  = integrSiteObject$design$filesIntegrations
  integrSiteObject$data$integrations = vector("list", length(fi))
  names(integrSiteObject$data$integrations) = names(fi)
  for(i in seq_along(fi)) {
    int_i = data.table :: as.data.table(read.table(file.path(wd, fi[[i]]), header = header))
    names(int_i) = variableNames
    if(length(excludeIntegrations) > 0L) {
      whichExcluded = integer()
      for(j in seq_along(excludeIntegrations)) {
        whichExcluded = c(whichExcluded, which(int_i$chr == excludeIntegrations[[j]][1] & int_i$start >= excludeIntegrations[[j]][2] & int_i$start < excludeIntegrations[[j]][3]))
      }
      if(length(whichExcluded) > 0L) int_i = int_i[-whichExcluded,]
    }
    integrSiteObject$data$integrations[[i]] = int_i
  } 
  
  return(integrSiteObject)
}

loadLocationsGenes = function(integrSiteObject, addPromoter = 0L, header = TRUE, exclDotInName = F, 
                              variableNames = c("chr", "start",	"end", "strand", "geneName", "chrInt", "testIndex")) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  
  pwd = integrSiteObject$infoStorage$pathWorkingDirectory
  flg = integrSiteObject$organism$fileInfoGenes
  loc = integrSiteObject$data$locationsGenes
  
  loc = data.table :: as.data.table(read.table(file.path(pwd, flg), header = header))
  loc$start[loc$strand == "+"] = loc$start[loc$strand == "+"] - addPromoter
  loc$end[loc$strand == "-"] = loc$end[loc$strand == "-"] + addPromoter
  
  if(exclDotInName) {
    withDotInName = grep("\\.", loc$geneName)
    if(length(withDotInName) > 0L) loc = loc[-withDotInName]
  }

  loc$testIndex = 1 : nrow(loc)
  colnames(loc) = variableNames
  data.table :: setkey(loc, chrInt, start, end)
  
  integrSiteObject$data$locationsGenes = loc
  return(integrSiteObject)
}

excludeLargestClone3 = function(integrSiteObject, selectedCovariates = character()) {
  cov = integrSiteObject$design$covariatesIntegrations
  int = integrSiteObject$data$integrations
  lge = integrSiteObject$data$locationsGenes
  
  eventsInGene = function(events) {
    eig = vector("list", nrow(lge))
    for(i in seq_len(nrow(lge))) {
      gene = lge[i]
      eig[[i]] = events[chr == gene$chr & start >= gene$start & start < gene$end]
    }
    return(eig)
  }
  
  if(length(selectedCovariates) > 0L) {
    unc = unique(cov[,selectedCovariates], by = selectedCovariates)
    idd = as.integer(apply(cov[,selectedCovariates], 1, function(x) base :: match(TRUE, apply(unc, 1, function(y) all(y == x)))))
    spi = base :: split(lapply(int, function(x) x[,c("chr", "start", "count")]), idd)
    spl = lapply(spi, function(x) do.call("rbind", x))
    spg = lapply(spl, eventsInGene)
  }
  else {
    spl = list(`1` = do.call("rbind", lapply(int, function(x) x[,c("chr", "start", "count")])))
  }
  return(integrSiteObject)
}

mergeData = function(integrSiteObject, reorder = FALSE) {
  convertChromosomeToInteger = function(chromosome) {
    label = substr(chromosome, start = 4L, stop = 5L)
    if(label == "X") return(23L)
    if(label == "Y") return(24L)
    if(label == "M") return(25L)
    return(as.integer(label))
  }

  chromosomes = integrSiteObject$organism$infoChromosomes$idCharacter
  id = integrSiteObject$design$id
  covariates = integrSiteObject$design$covariatesIntegrations
  integrations = integrSiteObject$data$integrations
  selcov = integrSiteObject$design$selectedCovariates
  
  group1 = integrSiteObject$data$group[["1"]]
  group2 = integrSiteObject$data$group[["2"]]
  #includedIntegrations = sort(c(group1, group2))
  whichGroup1 = sapply(id, function(x) x %in% group1)
  #isGroup2 = vector("integer", length(includedIntegrations)) 
  #isGroup2[whichGroup2] = 1L
  #isGroup2 = rep(1L, length(includedIntegrations)) 
  isGroup1 = rep(0L, length(id)) 
  isGroup1[whichGroup1] = 1L
  
  selectedId = id[sort(c(group1, group2))]
  
  tables = vector("list", length = length(selectedId))
  for(i in seq_along(tables)) {
    sid = selectedId[i]
    int_i = integrations[[sid]]
    integrationsInChromosomes = which(int_i$chr %in% chromosomes)
    tables[[i]] = data.table :: data.table(
      chr = int_i$chr[integrationsInChromosomes],
      #position = int_i$start[integrationsInChromosomes],
      start = int_i$start[integrationsInChromosomes],
      end = int_i$end[integrationsInChromosomes],
      strand = int_i$strand[integrationsInChromosomes],
      #strand = sapply(int_i$strand[integrationsInChromosomes], function(x) ifelse(x == "+", yes = 1, no = 0)),
      chrInt = sapply(int_i$chr[integrationsInChromosomes], convertChromosomeToInteger, USE.NAMES = FALSE),
      count = int_i$count[integrationsInChromosomes]
    )
    tables[[i]]$group = isGroup1[sid]
    tables[[i]]$series = sid
    
    relativeCounts = vector("integer", length(chromosomes))
    for(j in seq_along(relativeCounts)) relativeCounts[j] = sum(int_i$count[which(int_i$chr %in% chromosomes[j])])
    #tables[[i]]$relativeCounts = relativeCounts[tables[[i]]$chrInt]
    tables[[i]]$totalCounts = sum(relativeCounts)
    
    for(j in seq_along(selcov)) tables[[i]][[selcov[j]]] = covariates[sid, selcov[j]]
    
    #print(tables[[i]][which(chr == "chr6"  & start >= 35573584 & end < 35688915),])
    #print(tables[[i]][which(chr == "chr15" & start >= 60488283 & end < 61229302),])
  }
  combinedTable = do.call("rbind", tables) 
  
  if(reorder) combinedTable = combinedTable[order(chrInt, start, -group)]
  
  #setkey(combinedTable, chrInt, position)
  data.table :: setkey(combinedTable, chrInt, start, end)
  
  integrSiteObject$data$combined = combinedTable
  return(integrSiteObject)
}

excludeLargestClone = function(integrSiteObject, excludeInBothGroups = FALSE) {
  tab = integrSiteObject$data$combined
  gen = integrSiteObject$data$locationsGenes
  exl = data.table :: data.table(chr = character(), start = integer())
  for(i in 1:nrow(gen)) {
    ge = gen[i, ]
    rig = which((tab$chr == ge$chr) & (tab$start >= ge$start) & (tab$end < ge$end))
    iig = tab[rig, ] # integrations in gene ge
    va = c("start", "group"[excludeInBothGroups])
    
    if(nrow(iig) > 0L) {
      #df = iig[, c(va, "count")] # not working
      sel = c(va, "count")
      df = data.table::`[.data.table`(iig, , ..sel)
      
      ct = plyr::count(df = df, vars = va, wt_var = "count")
      if(excludeInBothGroups) {
        c1 = ct[which(ct$group == 1L), ]
        c2 = ct[which(ct$group == 2L), ]
        if(nrow(c1) > 0L) {
          e1 = which.max(c1$freq)
          exl = rbind(exl, data.table :: data.table(chr = ge$chr, start = c1$start[e1]))
        }
        if(nrow(c2) > 0L) {
          e2 = which.max(c2$freq)
          exl = rbind(exl, data.table :: data.table(chr = ge$chr, start = c2$start[e2]))
        }
      } else {
        if(nrow(ct) > 0L) {
          e0 = which.max(ct$freq)
          exl = rbind(exl, data.table :: data.table(chr = ge$chr, start = ct$start[e0]))
        }
      }
    }
  }
  #print(tab)
  #print(exl)
  #integrSiteObject$data$combined = tab[!exl, on = .(chr, start)]
  integrSiteObject$data$combined <- data.table :: `[.data.table`(tab, !exl, on = .(chr, start))
  return(integrSiteObject)
}



setFactorsCombinedData = function(integrSiteObject, ordered = FALSE) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  sc = integrSiteObject$design$selectedCovariates
  sf = integrSiteObject$design$selCovAreFactors
  cf = vector("list", length(sc))
  #va = c(c("chr", "strand", "group"), sc[which(sf)])
  va = c(c("chr", "strand"), sc[which(sf)])
  if(!all(va %in% colnames(integrSiteObject$data$combined))) stop("error")
  for(i in 1:length(va)) integrSiteObject$data$combined[[va[i]]] = factor(integrSiteObject$data$combined[[va[i]]], ordered = ordered)  
  if(length(sc) > 0L) for(j in 1:length(sc)) if(sf[j]) cf[[j]] = contrasts(integrSiteObject$data$combined[[sc[j]]])
  integrSiteObject$model$contrastFactors = cf
  return(integrSiteObject)
}

reorderCombinedData = function(integrSiteObject, groups = FALSE) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  if(groups) integrSiteObject$data$combined = integrSiteObject$data$combined[order(chrInt, start, group)]
  else integrSiteObject$data$combined = integrSiteObject$data$combined[order(chrInt, start)]
  return(integrSiteObject)
}


computeSplitIntegrationsGenes = function(integrSiteObject, nminSplit = 2L, removeDuplicates = TRUE) {
  interactionsGene_xIntegration_y = data.table :: foverlaps(x = integrSiteObject$data$locationsGenes, y = integrSiteObject$data$combined, type = "any", which = TRUE)
  interactionsGene_xIntegration_y = interactionsGene_xIntegration_y[!is.na(yid),]
  splitOverGenes = split(interactionsGene_xIntegration_y, by = "xid")
  #genesForIntegrations = lapply(splitOverGenes, function(x) x$yid)
  if(removeDuplicates) bb1 = !(duplicated(lapply(splitOverGenes, function(x) x$yid)))
  else bb1 = TRUE
  bb = bb1 & (lapply(splitOverGenes, nrow) > nminSplit)
  integrSiteObject$data$splitOverGenes = splitOverGenes[bb]
  
  return(integrSiteObject)
}


computeResponse = function(integrSiteObject) {
  integrSiteObject$model$response = cbind(integrSiteObject$data$combined$count, integrSiteObject$data$combined$totalCounts - integrSiteObject$data$combined$count)
  return(integrSiteObject)
}

computeDesignLogistic = function(integrSiteObject, testedCovariate = character(), otherCovariates = character(), differential = FALSE) {
  #if(length(covariates) > 0L) if(!all(covariates %in% colnames(integrSiteObject$data$combined))) stop("error")
  gbc = integrSiteObject$model$groupedByCovariate
  covariates = c(otherCovariates[which(otherCovariates != testedCovariate)], testedCovariate, "group"[(!gbc) & differential])
  f = as.formula(paste0("~", paste0(covariates, collapse = "+")))
  X = model.matrix(f, data = integrSiteObject$data$combined)
  integrSiteObject$model$design = X
  integrSiteObject$model$namesParameters = colnames(X)
  return(integrSiteObject)
}

fitLogisticRegressionUnderNull = function(integrSiteObject) {
  y = integrSiteObject$model$response
  X = integrSiteObject$model$design
  m = stats :: glm.fit(X, y, family = binomial())
  integrSiteObject$model$estimatedH0 = m
  return(integrSiteObject)
}

computeLengthInvestigatedGenome = function(integrSiteObject) {
  chromosomes = integrSiteObject$organism$infoChromosomes$idCharacter
  geneDet = integrSiteObject$data$locationsGenes[integrSiteObject$data$locationsGenes[["chr"]] %in% chromosomes,]
  geneDet = unique(geneDet, by = c("chr", "start", "end"))
  geneDet$length = geneDet$end - geneDet$start
  geneDet = data.table :: setorder(geneDet, by = -length) 
  geneDet = unique(geneDet, by = "geneName")
  integrSiteObject$model$lengthInvestigatedGenome = sum(geneDet$length)
  #integrSiteObject$model$lengthInvestigatedGenome = sum(integrSiteObject$organism$infoChromosomes$length)
  return(integrSiteObject)
}

computeSumLengthOutcomeVariable = function(integrSiteObject) {
  #outcomeVariable = "count"
  combined = integrSiteObject$data$combined
  groups = integrSiteObject$data$groups
  seriesTotIS = c()
  seriesTotCIS = c()
  for(i in c(groups[["1"]], groups[["2"]])) {
    seriesTotIS = c(seriesTotIS, sum(combined[["count"]][combined$series == i])) 
    seriesTotCIS = c(seriesTotCIS, length(combined[["count"]][combined$series == i]))
  }
  integrSiteObject$model$sumsOutcome = seriesTotIS
  integrSiteObject$model$lengthsOutcome = seriesTotCIS
  return(integrSiteObject)
}

computeDesignLogLinear = function(integrSiteObject, dummyForDatasets = FALSE, otherCovariates = character(), typeRegressionModel = "logistic") {
  group1 = integrSiteObject$data$groups[["1"]]
  group2 = integrSiteObject$data$groups[["2"]]
  
  differential = (length(group2) != 0L)
  grpByCov = integrSiteObject$model$groupedByCovariate
  
  
  if(length(unique(integrSiteObject$data$combined$series)) == 1L) X = data.frame(event = c(1,0,1,0), inZ = c(1,1,0,0), intercept = c(1,1,1,1))
  else {
    computeSeriesMat = function(start, increment, nrow, ncol, names) {
      seriesMat = matrix(0, nrow = nrow, ncol = ncol)
      #print(seriesMat)
      for(ii in 1:ncol) {
        i = ii + start
        tt = rep(0, nrow)
        tt[(increment * (i-1L) + 1L) : (increment * i)] = 1
        #print(tt)
        seriesMat[,ii] = tt
      }
      colnames(seriesMat) = names[1:ncol]
      return(seriesMat)
    }
    
    if(!differential) {
      nc = length(group1)
      if(typeRegressionModel == "logistic") {
        event = rep(c(1,0,1,0), nc)
        intercept = rep(rep(1,4), nc)
        inZ = rep(c(1,1,0,0), nc)
        if((nc > 1L) & dummyForDatasets) seriesMat = computeSeriesMat(start = 0L, increment = 4L, nrow = 4L * nc, ncol = nc - 1L, names = group1)
        else seriesMat = NULL
      }
      if(typeRegressionModel == "poisson") {
        event = rep(c(1,1), nc)
        intercept = rep(rep(1,2), nc)
        inZ = rep(c(1,0), nc)
        if((nc > 1L) & dummyForDatasets) seriesMat = computeSeriesMat(start = 0L, increment = 2L, nrow = 4L * nc, ncol = nc - 1L, names = group1)
        else seriesMat = NULL
      }
      if(dummyForDatasets) X = cbind(data.frame(event = event, inZ = inZ, intercept = intercept), seriesMat)
      else X = data.frame(event = event, inZ = inZ, intercept = intercept)
    }
    else {
      if(typeRegressionModel == "logistic") {
        event = c(rep(c(1,0,1,0), length(group1)), rep(c(1,0,1,0), length(group2)))
        
        nc = length(group1)
        nc1 = nc
        #if(nc >= 2L) seriesMat1 = computeSeriesMat(start = 0L,  increment = 4L, nrow = 4L * nc, ncol = nc - 1L, names = group1)
        if((nc >= 2L) & dummyForDatasets) seriesMat1 = computeSeriesMat(start = 0L,  increment = 4L, nrow = length(event), ncol = nc - 1L, names = group1)
        else seriesMat1 = NULL
        nc = length(group2)
        #if(nc >= 2L) seriesMat2 = computeSeriesMat(start = nc1, increment = 4L, nrow = 4L * nc, ncol = nc - 1L, names = group2)
        if((nc >= 2L) & dummyForDatasets) seriesMat2 = computeSeriesMat(start = nc1, increment = 4L, nrow = length(event), ncol = nc - 1L, names = group2)
        else seriesMat2 = NULL
        
        intercept = c(rep(rep(1,4), length(group1)), rep(rep(1,4), length(group2)))
        inZ = c(rep(c(1,1,0,0), length(group1)), rep(c(1,1,0,0), length(group2)))
        group = c(rep(rep(1,4), length(group1)), rep(rep(0,4), length(group2)))
      }
      if(typeRegressionModel == "poisson") {
        event = c(rep(c(1,1), length(group1)), rep(c(1,1), length(group2)))
        
        nc = length(group1)
        nc1 = nc
        #if(nc >= 2L) seriesMat1 = computeSeriesMat(start = 0L,  increment = 2L, nrow = 4L * nc, ncol = nc - 1L, names = group1)
        if((nc >= 2L) & dummyForDatasets) seriesMat1 = computeSeriesMat(start = 0L,  increment = 2L, nrow = length(event), ncol = nc - 1L, names = group1)
        else seriesMat1 = NULL
        nc = length(group2)
        #if(nc >= 2L) seriesMat2 = computeSeriesMat(start = nc1, increment = 2L, nrow = 4L * nc, ncol = nc - 1L, names = group2)
        if((nc >= 2L) & dummyForDatasets) seriesMat2 = computeSeriesMat(start = nc1, increment = 2L, nrow = length(event), ncol = nc - 1L, names = group2)
        else seriesMat2 = NULL
        
        intercept = c(rep(rep(1,2), length(group1)), rep(rep(1,2), length(group2)))
        inZ = c(rep(c(1,0), length(group1)), rep(c(1,0), length(group2)))
        group = c(rep(rep(1,2), length(group1)), rep(rep(0,2), length(group2)))
      }
      inZ_group = inZ * group
      sm = cbind(seriesMat1, seriesMat2)
      if((nc >= 2L) & dummyForDatasets) X = cbind(data.frame(event = event, inZ_group = inZ_group, intercept = intercept, inZ = inZ), sm)
      else X = data.frame(event = event, inZ_group = inZ_group, intercept = intercept, inZ = inZ)
      #X = cbind(data.frame(event = event, inZ_group = inZ_group, intercept = intercept, inZ = inZ, group = group), sm)
    }
  }
  
  if(length(otherCovariates) > 0L) {
    ci = integrSiteObject$design$covariatesIntegrations
    sc = integrSiteObject$design$selectedCovariates # argument not used
    sf = integrSiteObject$design$selCovAreFactors
    cf = integrSiteObject$model$contrastFactors
    
    nc = colnames(ci)
    
    desCovMatrix = function(wc) {
      cnFac = cf[[wc]]
      nameCov = nc[wc]
      isFac = sf[wc]
      
      if(isFac) {
        df = data.frame()
        for(i in 1:nrow(ci)) {
          lev_i = ci[i,wc]
          
          for(j in 1:4) df = rbind(df, cnFac[as.character(lev_i),])
          if(i == 1L) nam_1 = colnames(cnFac)
        }
        nameFac = paste0(nameCov, " : ", nam_1)
        colnames(df) = nameFac
        return(df)
      }
      else {
        df = data.frame(x = as.numeric(rep(ci[[wc]], each = 4L)))
        colnames(df) = nameCov
        return(df)
      }
    }
    
    colBlocks = lapply(1:length(sc), desCovMatrix)
    X = cbind(X, do.call("cbind", colBlocks))
    
    if(differential & (!grpByCov)) X$group = group
    
  }
  else if(differential) X$group = group
  
  integrSiteObject$model$typeRegressionModel = typeRegressionModel
  integrSiteObject$model$design = X
  
  return(integrSiteObject)
}


computeTestStatistics = function(integrSiteObject, testedCovariate = character(0), testType = "likRt", numberCores = 1L, sizeChunks = 1L) {
  if(class(integrSiteObject) != "integrSiteObject") stop("error")
  if(!is.numeric(numberCores) | !is.numeric(sizeChunks)) stop("error")
  #if(!(testType %in% c("likRt", "score"))) stop("error")
  
  splitOverGenes = integrSiteObject$data$splitOverGenes
  integrations  = integrSiteObject$data$combined
  response = integrSiteObject$model$response
  design = integrSiteObject$model$design
  modelH0 = integrSiteObject$model$estimatedH0
  
  idSplitIntegr = unique(c(seq(1, length(splitOverGenes), by = sizeChunks), length(splitOverGenes)))
  idSplitCores = 2L
  testInd = 0L
  tests = array(list(), length(idSplitIntegr))
  
  while(idSplitCores <= length(idSplitIntegr)){
    testInd = testInd + 1L
    print(paste(idSplitCores, "...", length(idSplitIntegr)))
    splitOverCores = array(list(), numberCores)
    for(jk in c(1:numberCores)){
      if(idSplitCores <= length(idSplitIntegr)) {
        splitOverCores[[jk]] = splitOverGenes[idSplitIntegr[idSplitCores - 1L] : (idSplitIntegr[idSplitCores] - 1L)]  
        idSplitCores=idSplitCores + 1L
      }
      else {
        splitOverCores = splitOverCores[!(unlist(lapply(splitOverCores,is.null)))]
      } 
    }
    if(testType=="likRt"){
      if(integrSiteObject$model$twoGroups) test = parallel :: mclapply(splitOverCores, function(X) computeLikRtTest2gr(X, integrations, testedCovariate, response, design), mc.cores = numberCores)
      else test = parallel :: mclapply(splitOverCores, function(X) computeLikRtTest1gr(X, integrations, testedCovariate, response, design, modelH0), mc.cores = numberCores)
      tests[[testInd]] = t(do.call("cbind", test))
    }
    if(testType=="score"){
      if(integrSiteObject$model$twoGroups) stop("score test can be computed only for one group")
      else test = parallel :: mclapply(splitOverCores, function(X) computeScoreTest(X, integrations, testedCovariate, modelH0), mc.cores = numberCores)
      tests[[testInd]] = test
    }
    if(testType=="robst"){
      if(integrSiteObject$model$twoGroups) test = parallel :: mclapply(splitOverCores, function(X) computeRobstTest2gr(X, integrations, testedCovariate, response, design), mc.cores = numberCores)
      else test = parallel :: mclapply(splitOverCores, function(X) computeRobstTest1gr(X, integrations, testedCovariate, response, design, modelH0), mc.cores = numberCores)
      tests[[testInd]] = test
    }
  }
  return(tests)
}

computeLikRtTest2gr = function(included, integrations, testedCovariate, response, design) {
  computeTest = function(x) {
    newCovariate = rep(0.0, length = nrow(integrations))
    newCovariate[x$yid] = integrations[[testedCovariate]][x$yid]
    #groupInteraction = newCovariate * as.numeric(as.numeric(integrations[["group"]]) == 1)
    groupInteraction = newCovariate * integrations[["group"]]
    newDesign = as.matrix(cbind(design, newCovariate, groupInteraction))
    fit1 = fastglm :: fastglmPure(newDesign, response, family = binomial(), method = 3,tol = 1e-5)
    fit0 = fastglm :: fastglmPure(newDesign[,-ncol(newDesign), drop = FALSE], response, family = binomial(), method = 3, tol = 1e-5)
    tstat = deviance(fit0) - deviance(fit1)
    return(as.numeric(c(as.numeric(fit1$coefficients), tstat)))
  }
  return(do.call("cbind", lapply(included, computeTest)))  
}

computeLikRtTest1gr = function(included, integrations, testedCovariate, response, design, modelH0) {
  computeTest = function(x) {
    newCovariate = rep(0, length = nrow(integrations))
    newCovariate[x$yid] = integrations[[testedCovariate]][x$yid]
    newDesign = as.matrix(cbind(design, newCovariate))
    fit1 = fastglm :: fastglmPure(newDesign, response, family=binomial())
    tstat = deviance(modelH0) - deviance(fit1)
    return(as.numeric(c(as.numeric(fit1$coefficients), tstat)))
  }
  return(do.call("cbind", lapply(included, computeTest)))
}

computeRobstTest2gr = function(included, integrations, testedCovariate, response, design) {
  computeTest = function(x) {
    minForRobust = 2L
    robustInteraction = 0.0
    robustAnalysis = FALSE
    robustTestStat = 0.0
    seriesAll = integrations[["series"]]
    if(length(x$yid) > 1L) {
      starts = integrations[["start" ]][x$yid]
      counts = integrations[["count" ]][x$yid]
      series = integrations[["series"]][x$yid]
      ust = unique(starts)
      if(length(ust) >= minForRobust) {
        df = data.table :: data.table(starts = starts, counts = counts)
        ct = plyr::count(df = df, vars = "starts", wt_var = "counts")
        mc = which.max(ct$freq)
        ex = x$yid[which(starts == ct$starts[mc])]
        we = which(x$yid %in% ex)
        ce = counts[we]
        se = series[we]
        
        newCovariate = rep(0, length = nrow(integrations))
        newCovariate[x$yid] = integrations[[testedCovariate]][x$yid]
        #groupInteraction = newCovariate * as.numeric(as.numeric(integrations[["group"]]) == 1)
        groupInteraction = newCovariate * integrations[["group"]]
        newDesign0 = as.matrix(cbind(design, newCovariate))[-ex,]
        newDesign1 = as.matrix(cbind(design, newCovariate, groupInteraction))[-ex,]
        #newDesign0 = as.matrix(design[-ex,])
        #newDesign1 = as.matrix(cbind(design, newCovariate))[-ex,]
        newResponse = response
        for(j in seq_along(se)) {
          pos_j = which(seriesAll == se[j])
          newResponse[pos_j,2] = newResponse[pos_j,2] - ce[j]
        }
        newResponse = newResponse[-ex,]
        
        #fit0 = stats :: glm.fit(newDesign0, newResponse, family = binomial())
        fit0 = fastglm :: fastglmPure(newDesign0, newResponse, family = binomial())
        fit1 = fastglm :: fastglmPure(newDesign1, newResponse, family = binomial())
        
        robustAnalysis = TRUE
        robustInteraction = tail(as.numeric(fit1$coefficients), 1)
        robustTestStat = deviance(fit0) - deviance(fit1)
      }
    }
    return(as.numeric(c(robustAnalysis, robustInteraction, robustTestStat)))
  }
  return(do.call("cbind", lapply(included, computeTest)))
}

computeRobstTest1gr = function(included, integrations, testedCovariate, response, design, modelH0) {
  computeTest = function(x) {
    minForRobust = 2L
    robustAnalysis = FALSE
    robustInteraction = 0.0
    robustTestStat = 0.0
    seriesAll = integrations[["series"]]
    if(length(x$yid) > 1L) {
      starts = integrations[["start" ]][x$yid]
      counts = integrations[["count" ]][x$yid]
      series = integrations[["series"]][x$yid]
      ust = unique(starts)
      if(length(ust) >= minForRobust) {
        df = data.table :: data.table(starts = starts, counts = counts)
        ct = plyr::count(df = df, vars = "starts", wt_var = "counts")
        mc = which.max(ct$freq)
        ex = x$yid[which(starts == ct$starts[mc])]
        we = which(x$yid %in% ex)
        ce = counts[we]
        se = series[we]
        
        newCovariate = rep(0, length = nrow(integrations))
        newCovariate[x$yid] = integrations[[testedCovariate]][x$yid]
        
        newDesign0 = as.matrix(design[-ex,])
        newDesign1 = as.matrix(cbind(design, newCovariate))[-ex,]
        newResponse = response
        for(j in seq_along(se)) {
          pos_j = which(seriesAll == se[j])
          newResponse[pos_j,2] = newResponse[pos_j,2] - ce[j]
        }
        newResponse = newResponse[-ex,]
        
        #fit0 = stats :: glm.fit(newDesign0, newResponse, family = binomial())
        fit0 = fastglm :: fastglmPure(newDesign0, newResponse, family = binomial())
        fit1 = fastglm :: fastglmPure(newDesign1, newResponse, family = binomial())
        
        robustAnalysis = TRUE
        robustInteraction = tail(as.numeric(fit1$coefficients), 1)
        robustTestStat = deviance(fit0) - deviance(fit1)
      }
    }
    return(as.numeric(c(robustAnalysis, robustInteraction, robustTestStat)))
  }
  return(do.call("cbind", lapply(included, computeTest)))
}

computeScoreTest = function(included, integrations, testedCovariate, modelH0) {
  computeNewCovariate = function(x) {
    newCovariate = rep(0, length = nrow(integrations))
    newCovariate[x$yid] = integrations[[testedCovariate]][x$yid]
    return(newCovariate)
  }
  return(statmod :: glm.scoretest(modelH0, do.call("cbind", lapply(included, computeNewCovariate))))
}

computePvalueLikRt = function(tests, correctionMethod = "fdr", namesParameters = character()) {
  tests = do.call("rbind", tests)
  #tests = unique(tests, by = "geneName")
  testIndex = as.numeric(rownames(tests))
  tests = data.table :: as.data.table(tests)
  #colnames(tests)[1:(ncol())] # names design
  #colnames(tests)[ncol(tests) - 1L] # name covariate
  if(length(namesParameters) > 0L) colnames(tests)[1 : length(namesParameters)] = namesParameters
  colnames(tests)[ncol(tests) - 1L] = "interaction"
  colnames(tests)[ncol(tests)] = "testStat"
  tests$signedTestStat = tests$testStat * sign(tests$interaction)
  tests$testIndex = testIndex
  tests$pvalue = (1 - pchisq(abs(tests$testStat), df = 1))
  tests$adjpvalue = p.adjust(tests$pvalue, correctionMethod)
  return(tests)
}

computePvalueRobst = function(tests, correctionMethod = "fdr") {
  tests = tests[sapply(tests, function(x) length(x) > 0L)]
  tests = lapply(tests, function(x) do.call("cbind", x))
  tests = do.call("cbind", tests)
  #tests = unique(tests, by = "geneName")
  testIndex = as.numeric(colnames(tests))
  tests = data.table :: as.data.table(t(tests))
  tests$V1 = as.logical(tests$V1)
  colnames(tests)[1] = "robustAnalysis"
  colnames(tests)[2] = "robustInteraction"
  colnames(tests)[3] = "robustTestStat"
  tests$testIndex = testIndex
  tests$robustSignedTestStat = sign(tests$robustInteraction) * tests$robustTestStat
  tests$robustPvalue = (1 - pchisq(abs(tests$robustTestStat), df = 1))
  tests$robustAdjPvalue = 1.0
  tests$robustAdjPvalue[which(tests$robustAnalysis)] = p.adjust(tests$robustPvalue[which(tests$robustAnalysis)], correctionMethod)
  return(tests)
}

computePvalueScore = function(tests, correctionMethod = "fdr") {
  tests = unlist(tests)
  tests = data.table :: data.table(testIndex = as.numeric(names(tests)), testStat = as.numeric(tests))
  tests$pvalue = 2 * (1 - pnorm(abs(tests$testStat)))
  tests$adjpvalue = p.adjust(tests$pvalue, correctionMethod)
  return(tests)
}


computeTestStatisticLogLinear = function(integrSiteObject, numberCores = 1L) {
  #if(!(outcomeVariable %in% colnames(integrSiteObject$data$combined))) stop("error")
  splitOverGenes = integrSiteObject$data$splitOverGenes
  integrations  = integrSiteObject$data$combined
  locationsGenes = integrSiteObject$data$locationsGenes
  groups = c(integrSiteObject$data$groups[["1"]], integrSiteObject$data$groups[["2"]])
  X = integrSiteObject$model$design
  typeRegressionModel = integrSiteObject$model$typeRegressionModel
  #numIntegrations = integrSiteObject$model$sumsOutcome
  numIntegrations = integrSiteObject$model$lengthsOutcome
  seriesTotCIS = integrSiteObject$model$lengthsOutcome
  lengthInvestigatedGenome = integrSiteObject$model$lengthInvestigatedGenome
  
  computeTest = function(osp) {
    gene = locationsGenes[osp$xid[1],]
    geneLength = gene$end - gene$start
    #integrationsInGene = integrations[osp$yid, c("series", outcomeVariable)]
    integrationsInGene = integrations[osp$yid, c("series", "count")]
    #numIntegrInGene = unlist(lapply(groups, function(x, y) sum(y[[outcomeVariable]][y$series == x]), integrationsInGene))
    #numIntegrInGene = unlist(lapply(groups, function(x, y) sum(y[["count"]][y$series == x]), integrationsInGene))
    #numIntegrInSeries = unlist(lapply(groups, function(x, y) sum(y$series == x), integrationsInGene))
    numIntegrInGene = unlist(lapply(groups, function(x, y) length(y[["count"]][y$series == x]), integrationsInGene))
    numIntegrInSeries = unlist(lapply(groups, function(x, y) length(y$series == x), integrationsInGene))
    if(typeRegressionModel == "logistic") {
      numIntegrOutsideGene = numIntegrations - numIntegrInGene
      outZnIS = lengthInvestigatedGenome - geneLength - numIntegrOutsideGene
      freQ = as.vector(rbind(numIntegrInGene, geneLength - numIntegrInGene, numIntegrOutsideGene, outZnIS))
      tstat = logLinearRegMultLogistic(frequencies = freQ, designMatrix = X)
    }
    if(typeRegressionModel == "poisson") {
      numIntegrOutsideGene = numIntegrations - numIntegrInGene
      numIntegrInGene=numIntegrInGene/numIntegrInSeries
      numIntegrInGene[is.na(numIntegrInGene)]=0
      numIntegrOutsideGene=numIntegrOutsideGene/(seriesTotCIS-(numIntegrInSeries))
      of = as.vector(rbind(numIntegrInSeries,(seriesTotCIS-(numIntegrInSeries))))
      freQ = as.vector(rbind(numIntegrInGene,numIntegrOutsideGene))
      tstat = logLinearRegMultPoisson(frequencies = freQ, designMatrix = X, weightsPoisson = of)
    }
    retV = list(gene = gene, nIS = freQ, Tstat = tstat)
    return(retV)
  }
  
  tests = parallel :: mclapply(splitOverGenes, computeTest, mc.cores = numberCores)
  tests = cbind(do.call("rbind", lapply(tests,function(g) g[["gene"]])),
                do.call("rbind", lapply(tests,function(g) g[["nIS"]])),
                do.call("rbind", lapply(tests,function(g) g[["Tstat"]])))
  
  if(length(tests) > 0L) {
    colnames(tests)[ncol(tests) - 1L] = "Zscore"
    colnames(tests)[ncol(tests)] = "testStat"
    tests$signedTestStat = tests$testStat * sign(tests$Zscore)
  } 
  return(tests)
}

logLinearRegMultLogistic = function(frequencies, designMatrix, weightsPoisson = NULL) {
  fit = stats :: glm.fit(x = as.matrix(designMatrix[,-c(1)]), y = as.matrix(designMatrix[,1], ncol=1), weights = frequencies, family = binomial())
  if(fit$converged == T){
    fit1 = stats :: glm.fit(x = as.matrix(designMatrix[,-c(1,2)]), y = as.matrix(designMatrix[,1], ncol=1), weights = frequencies, family = binomial())
    t = deviance(fit1) - deviance(fit)
  }
  else{
    fit = brglm :: brglm.fit(x = designMatrix[,-1], y = designMatrix[,1], weights = frequencies, family = binomial(), intercept = F)
    fit1 = brglm :: brglm.fit(x = designMatrix[,-c(1,2)], y = designMatrix[,1], weights = frequencies, family = binomial(), intercept = F)
    t = fit1$deviance - fit$deviance
  }
  return(as.numeric(c(as.numeric(fit$coefficients[1]), t)))
}

logLinearRegMultPoisson = function(frequencies, designMatrix, weightsPoisson = NULL) {
  if(fit$converged == T) {
    fit1 = stats :: glm.fit(x = as.matrix(designMatrix[,-c(1,2)]), y = as.matrix(frequencies, ncol=1), weights = weightsPoisson, family = poisson())
    t = deviance(fit1) - deviance(fit)
  }
  else {
    t = NA
  }
  return(as.numeric(c(as.numeric(fit$coefficients[1]), t)))
}

computePvalueLgLin = function(tests, correctionMethod = "fdr") {
  if(length(tests) > 0L) {
    tests$pvalue = 1 - pchisq(abs(tests$testStat), df = 1)
    tests$adjpvalue = p.adjust(tests$pvalue, method = correctionMethod)
  }
  return(tests)
}

processTests = function(integrSiteObject, tests, rmDuplicateGenes = TRUE, isTarget = FALSE, overOnly = FALSE) {
  if(length(tests) > 0L) {
    if(!isTarget) tests = merge(integrSiteObject$data$locationsGenes, tests, by = "testIndex")
    if(overOnly) tests = tests[order(tests$testStat,decreasing = T),]
    else tests = tests[order(tests$testStat,decreasing = T),]
    if(rmDuplicateGenes) tests = unique(tests, by = "geneName")
    if(!isTarget) tests = tests[!is.na(testStat),]
    class(tests) = c("procTestsISObj", class(tests))
  }
  return(tests)
}


