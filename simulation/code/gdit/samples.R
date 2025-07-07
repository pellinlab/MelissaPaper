# gdit

cp = read.csv(paste0(folCode, "../data/chpt.csv"), header = F, sep = ",")$V1
whichSim = read.csv(paste0(folCode, "../data/whsi.csv"), header = F)$V1
prob = read.csv(paste0(folCode, "../data/prob.csv"), header = F, sep = ",")$V1

prob2 = prob1 = prob
prob2[whichSim] = prob1[whichSim] * heightHotspot

validSample = FALSE
ks = 1L

while(!validSample & ks <= limIterSamples) {
  for(j in includedSamples) {
    locations = rDiscreteUnifMix(n = sizeSample, cp, 
                                 prob1, 
                                 replace = FALSE, sorted = TRUE)
    sampled = data.frame(chr = "chr15", 
                         start = locations, end = locations + 1L, 
                         strand = c("+","-")[1 + rbinom(length(locations), 1, .5)], 
                         count = 1L)
    nd = setNumberDigits(i)
    nameFile = paste0("s", c("", "0", "00", "000")[1 + numberDigits - nd], i)
    write.table(sampled, 
                paste0(folSamples, nameFile, "s", j, "g1.bed"), 
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    itit = sapply(lapply(1:5, function(x) whInTarget(locations, whGene = x)), length)
    itl1[[i]][[j]] = itit
    nl1[i,j] = sum(itit)
  }
  
  for(j in includedSamples) {
    locations = rDiscreteUnifMix(n = sizeSample, cp, 
                                 prob2, 
                                 replace = FALSE, sorted = TRUE)
    sampled = data.frame(chr = "chr15", 
                         start = locations, end = locations + 1L, 
                         strand = c("+","-")[1 + rbinom(length(locations), 1, .5)], 
                         count = 1L)
    nd = setNumberDigits(i)
    nameFile = paste0("s", c("", "0", "00", "000")[1 + numberDigits - nd], i)
    write.table(sampled, 
                paste0(folSamples, nameFile, "s", j, "g2.bed"), 
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    itit = sapply(lapply(1:5, function(x) whInTarget(locations, whGene = x)), length)
    itl2[[i]][[j]] = itit
    nl2[i,j] = sum(itit)
  }
  
  if(sum(nl1[i,]) >= 1L & sum(nl2[i,]) >= 1L) validSample = TRUE
  ks = ks+1L
}



