mainFolder = "~/Downloads/MELISSApaper"

type = 1
size = 3 
heig = 5 
### t : [1:"govt", 2:"gdit", 3:"clfi", 4:"cdif"]
### s : 100[1000] * s
### h : 2 ^ (h-1)

nSimulations = 2L

limIterSamples = 20L
addPromoter = 1000L
dummyForDatasets = c(T, F)[1+type %in% c(3,4)]
pvalueCorrectionMethod = "fdr"
nminSplit = 1L
nCores = 8L
sizeChunks = 100L

startLoop = 1L
numberDigits = 4L

typeAnalysis = c("govt", "gdit", "clfi", "cdif")[type]
isCloneFitness = typeAnalysis %in% c("clfi", "cdif")
isDiffAnalysis = typeAnalysis %in% c("gdit", "cdif")
sizeSample = as.integer(c(1,2,4,8)[size] * c(100, 1000)[1 + isCloneFitness])
heightHotspot = as.integer(c(1,2,4,8,16,32,64)[heig])
includedSamples = list(c(1:4), c(1:6))[[1 + isCloneFitness]]

folCode = paste0(mainFolder, "/simulation/code/")
folSim  = paste0(mainFolder, "/simulation/simu/")

idsTyA = c("govt", "gdit", "clfi", "cdif")
idsHeH = c("h01", "h02", "h04", "h08", "h16", "h32", "h64")
idsSiS = c("s01", "s02", "s04", "s08")
idsGen = c("20K", "10K", "80K", "40K", "05K")
idSimulation = paste0(idsTyA[type], "/", idsSiS[size], "/", idsHeH[heig])

source(paste0(folCode, "../../package/intersectionSiteObject2.R"))
source(paste0(folCode, "../../package/intersectionSiteObject4.R"))
source(paste0(folCode, "sampleUniform.R"))
source(paste0(folCode, "simutils.R"))

source(paste0(folCode, "main.R"))
source(paste0(folCode, "store.R"))
