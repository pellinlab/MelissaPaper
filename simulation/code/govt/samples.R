# govt

cp = c(1L, targetGenes(), 101973591L + 100000L)
whichSim = 2L * (1:5)
prob = c(1, rep(c(heightHotspot, 1), times = 5))

validSample = FALSE
ks = 1L

while(!validSample & ks <= limIterSamples) {
  for(j in includedSamples) {
    integrations = rDiscreteUnifMix(n = sizeSample, cp, prob, 
                                    replace = FALSE, sorted = TRUE)
    sampled = data.frame(chr = "chr15", 
                         start = integrations, end = integrations + 1L, 
                         strand = c("+","-")[1 + rbinom(length(integrations), 1, .5)], 
                         count = 1L)
    nd = setNumberDigits(i)
    nameFile = paste0("s", c("", "0", "00", "000")[1 + numberDigits - nd], i,"s",j)
    write.table(sampled, 
                paste0(folSamples, nameFile,".bed"), 
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    itit = sapply(lapply(1:5, function(x) whInTarget(integrations, whGene = x)), length)
    itl1[[i]][[j]] = itit
    nl1[i,j] = sum(itit)
  }
  
  if(sum(nl1[i,]) >= 1L) validSample = TRUE
  ks = ks+1L
}
