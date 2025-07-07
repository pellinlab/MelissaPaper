mainFolder = "~/Downloads/MELISSApaper"

source(paste0(mainFolder, "/simulation/code/simutils.R"))

typeAnalysis = 1
isRecursive = T

idType = c("govt", "gdit", "clfi", "cdif")
idSize = c("s01", "s02", "s04", "s08")
idHeig = c("h01", "h02", "h04", "h08", "h16", "h32", "h64")
idGene = c("g20", "g10", "g80", "g40", "g05")

folAnalysis = paste0(mainFolder, "/simulation/simu/", idType[typeAnalysis], "/")

numberDigits = 4L
numberSimulatedValues = 1000L
thresholdPv = 0.05

ppv = vector("list", length(idSize))
cov = vector("list", length(idGene))
names(ppv) = idSize
names(cov) = idGene
fwt = cov

i = j = h = 1L
for(i in seq_along(idSize)) {
  print(i)
  ppv[[i]] = vector("list", length(idHeig))
  names(ppv[[i]]) = idHeig
  for(j in seq_along(idHeig)) {
    ppv[[i]][[j]] = vector("numeric", numberSimulatedValues)
    
    for(h in 1:numberSimulatedValues) {
      nd = setNumberDigits(h)
      idFile = paste0(c("", "0", "00", "000")[1 + numberDigits - nd], h)
      re = read.csv(paste0(folAnalysis, idSize[i], "/", idHeig[j], "/", c("results/", "recress/")[1 + isRecursive], "r", idFile, ".csv"), header = TRUE, sep = "\t")
      gen = re$geneName
      if(typeAnalysis == 1L) sig = re$Zscore
      apv = re$adjpvalue
      
      whp = grep("SIMULA", gen)
      if(length(whp) > 0) {
        if(typeAnalysis == 1L) {
          nTruPos = sum((apv[whp] <  thresholdPv) & (sig[whp] > 0))
          nFalNeg = sum((apv[whp] >= thresholdPv) | (sig[whp] < 0))
          nAllPos = sum((apv < thresholdPv) & (sig > 0))
        }
        else {
          nTruPos = sum(apv[whp] <  thresholdPv)
          nFalNeg = sum(apv[whp] >= thresholdPv)
          nAllPos = sum(apv < thresholdPv)
        }
        if(nAllPos > 0) ppv[[i]][[j]][[h]] = nTruPos / nAllPos
        else ppv[[i]][[j]][[h]] = 0.0
      }
      else ppv[[i]][[j]][[h]] = 0.0 
      
      
      whp = sapply(paste0("SIMULA", c("20K-", "10K-", "80K+", "40K+", "05K-")), function(x) grep(x, gen))
      if(length(whp) > 0) {
        if(any(sapply(whp, length) > 1)) stop("error")
        for(g in seq_along(whp)) {
          if(length(whp[[g]]) > 0L) {
            fwt[[g]][[i]][[j]][h] = T
            if(typeAnalysis == 1L) cov[[g]][[i]][[j]][h] = (apv[whp[[g]]] <  thresholdPv) & (sig[whp[[g]]] > 0)
            else cov[[g]][[i]][[j]][h] = apv[whp[[g]]] <  thresholdPv
          }
        }
      }
    }
  }
}



for(i in 1:length(idSize)) {
  idFile = paste0(mainFolder, "/simulation/resu/ppv/", idType[typeAnalysis], "_", idSize[i], ".csv")
  write.csv(do.call("cbind", ppv[[i]]), idFile, row.names = F)
} 

idCov = paste0("g", 1:length(idGene))
idFwt = paste0("f", 1:length(idGene))
for(i in 1:length(idSize)) {
  for(j in 1:length(idGene)) {
    idFile = paste0(mainFolder, "/simulation/resu/cov/", idType[typeAnalysis], "_", idCov[j], "_", idSize[i], ".csv")
    write.csv(do.call("cbind", cov[[j]][[i]]), idFile, row.names = F)
    idFile = paste0(mainFolder, "/simulation/resu/cov/", idType[typeAnalysis], "_", idFwt[j], "_", idSize[i], ".csv")
    write.csv(do.call("cbind", fwt[[j]][[i]]), idFile, row.names = F)
  } 
} 

