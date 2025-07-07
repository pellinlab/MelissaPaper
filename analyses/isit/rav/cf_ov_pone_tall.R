#mainFolder = "~/Downloads/MELISSApaper"
#mainFolder = "~/Work/Boston/Melissa/MELISSApaper"
dataFolder = paste0(mainFolder, "/analyses/data/ravi")
source(paste0(mainFolder,"/package/intersectionSiteObject2.R"))
source(paste0(mainFolder,"/package/intersectionSiteObject4.R"))

integrSiteObject = createEmptyIntegrSiteObject()
integrSiteObject = setWorkingDirectory(integrSiteObject, mainFolder)
integrSiteObject = setOrganism(integrSiteObject, "human", "analyses/data/othr/hg38_genesLongestMinlen500.bed")

#library(data.table)
#library(stringr)

cell_typid = c("B", "T", "NK", "Neut", "Macr")#, "HSPC")
info_data = data.table::data.table(fileName = character(), fileId = character(), patient = character(), cellType = character(), timeChr = character(), timeMts = integer(), timeYrs = numeric())

i=j=1
for(i in seq_along(patients)) for(j in seq_along(cell_typid))
{
  pattern_id = paste0("^", patients[i], "_", cell_typid[j], "_",  ".*\\.bed$")
  files_id = list.files(path = dataFolder, pattern = pattern_id, full.names = F)
  files_pt = paste0("analyses/data/ravi/", files_id)
  
  times_id = paste0(stringr::str_pad(sub(".*_(\\d+)m.*", "\\1", files_id), width = 2, pad = "0"), "m")
  times_mt = as.integer(sub(".*_(\\d+)m.*", "\\1", files_id))
  times_yr = times_mt / 12
  
  id = data.table::data.table(fileName = files_pt, fileId = files_id, patient = patients[i], cellType = cell_typid[j], timeChr = times_id, timeMts = times_mt, timeYrs = times_yr)
  data.table::setorder(id, timeMts)
  
  info_data = rbind(info_data, id)
  info_data = info_data[timeMts <= time_limit]
}
rm(pattern_id, files_id, files_pt, times_id, times_mt, times_yr, id, i, j)

integrSiteObject = setDesign(integrSiteObject, filesIntegrations = info_data$fileName, covariatesIntegrations = as.data.frame(info_data[, .(time = timeYrs, type = cellType)]))

result = iSiteCloneFitness(integrSiteObject, 
                           testedCovariate = "time", 
                           idFiles = 1:nrow(info_data), 
                           addPromoter = 1000L,
                           selectedCovariates = c("time", "type"), 
                           selCovAreFactors = c(F, T), 
                           nminSplit = 1L,
                           numberCores = 15L,
                           robust = T)

#result$result
StoreResults(integrSiteObject, result, paste0(nameFileResult = "/analyses/results/singlePatient/cf_ov_", patients, "_tmax", time_limit,".csv"), storeTable = F)


