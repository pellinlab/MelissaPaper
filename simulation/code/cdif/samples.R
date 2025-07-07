# cdif

cp = read.csv(paste0(folCode, "../data/chpt.csv"), header = F, sep = ",")$V1
whichSim = read.csv(paste0(folCode, "../data/whsi.csv"), header = F)$V1
prob = read.csv(paste0(folCode, "../data/prob.csv"), header = F, sep = ",")$V1

nrowProb = 10000L
a = 1
b = c(0,2)
igr = 1:5

validSample = FALSE
ks = 1L

while(!validSample & ks <= limIterSamples) {
  baseInt = rDiscreteUnifMix(nrowProb, cp, prob, replace = F, sorted = T)
  
  growth = rep(1, length(baseInt))
  for(h in 1:(length(cp)-1)) {
    whichTarget = which(baseInt >= cp[h] & baseInt < cp[h+1L])
    if(h %in% whichSim) growth[whichTarget] = 1
    else growth[whichTarget] = runif(1, b[1], b[2])
  }
  
  growth1 = growth
  prob1 = matrix(0L, length(baseInt), 6L)
  prob1[,1] = 1L
  for(j in intersect(includedSamples, 2:6)) prob1[,j] = a + growth1 * igr[j-1]
  
  for(j in includedSamples) {
    sam = sample.int(nrow(prob1), sizeSample, replace = TRUE, prob = prob1[,j])
    tab = table(sam)
    locations = baseInt[as.integer(names(tab))]
    counts = as.integer(tab)
    if(j == 1L) counts = rep(1L, length(counts))
    
    lcit = lapply(1:5, function(x) whInTarget(locations, whGene = x))
    itit = sapply(lcit, length)
    itl1[[i]][[j]] = itit
    nl1[i,j] = sum(itit)
    nc1[i,j] = sum(unlist(lapply(lcit, function(x) counts[x])))
    #nl1[i,j] = sum(locations >= cp[whichSim] & locations < cp[whichSim+1L])
    #nc1[i,j] = sum(counts[locations >= cp[whichSim] & locations < cp[whichSim+1L]])
    
    s = data.frame(chr = "chr15", 
                   start = locations, end = locations + 1L, 
                   strand = c("+","-")[1 + rbinom(length(locations), 1, .5)], 
                   count = counts)
    nameFile = paste0("s", c("", "0", "00", "000")[1 + numberDigits - nd], i, "T", j, "G1")
    write.table(s, 
                paste0(folSamples, nameFile,".bed"), 
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  
  growth2 = growth
  for(w in seq_along(whichSim)) growth2[which(baseInt >= cp[whichSim[w]] & baseInt < cp[whichSim[w]+1L])] = heightHotspot
  #growth2[which(baseInt >= cp[whichSim] & baseInt < cp[whichSim+1L])] = heightHotspot
  prob2 = matrix(0L, length(baseInt), 6L)
  prob2[,1] = 1L
  for(j in intersect(includedSamples, 2:6)) prob2[,j] = a + growth2 * igr[j-1]
  
  for(j in includedSamples) {
    sam = sample.int(nrow(prob2), sizeSample, replace = TRUE, prob = prob2[,j])
    tab = table(sam)
    locations = baseInt[as.integer(names(tab))]
    counts = as.integer(tab)
    if(j == 1L) counts = rep(1L, length(counts))
    
    lcit = lapply(1:5, function(x) whInTarget(locations, whGene = x))
    itit = sapply(lcit, length)
    itl2[[i]][[j]] = itit
    nl2[i,j] = sum(itit)
    nc2[i,j] = sum(unlist(lapply(lcit, function(x) counts[x])))
    #nl2[i,j] = sum(locations >= cp[whichSim] & locations < cp[whichSim+1L])
    #nc2[i,j] = sum(counts[locations >= cp[whichSim] & locations < cp[whichSim+1L]])
    
    s = data.frame(chr = "chr15", 
                   start = locations, end = locations + 1L, 
                   strand = c("+","-")[1 + rbinom(length(locations), 1, .5)], 
                   count = counts)
    nameFile = paste0("s", c("", "0", "00", "000")[1 + numberDigits - nd], i, "T", j, "G2")
    write.table(s, 
                paste0(folSamples, nameFile,".bed"), 
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  
  if(sum(nl1[i,]) > 1L & sum(nl2[i,]) > 1L) validSample = TRUE
  ks = ks+1L
}


