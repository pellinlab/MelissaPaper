folAnalysis = paste0(folCode, typeAnalysis, "/")
folSamples = paste0(folSim, idSimulation, "/samples/")
folResults = paste0(folSim, idSimulation, "/results/")
folRecRess = paste0(folSim, idSimulation, "/recress/")

timesSim = vector("list",    nSimulations)
sampleIt = vector("integer", nSimulations)
nl1 = nl2 = nc1 = nc2 = matrix(NA_integer_, nSimulations, length(includedSamples))
itl1 = itl2 = itc1 = itc2 = replicate(nSimulations, vector("list", length(includedSamples)), simplify = FALSE)

i = 1L
for(i in startLoop:nSimulations) {
  cat("simulation", i, "/", nSimulations, "\n")
  nd = setNumberDigits(i)
  
  pt0 = proc.time()["elapsed"]
  cat("draw samples:\n")
  source(paste0(folAnalysis, "samples.R"))
  pt1 = proc.time()["elapsed"]
  cat("analysis:\n")
  source(paste0(folAnalysis, "analysis.R"))
  pt2 = proc.time()["elapsed"]
  cat("recursive analysis:\n")
  source(paste0(folAnalysis, "recursive.R"))
  pt3 = proc.time()["elapsed"]
  
  sampleIt[[i]] = ks - as.integer(validSample)
  timesSim[[i]] = c(pt0, pt1, pt2, pt3)
}

