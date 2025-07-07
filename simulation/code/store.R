write.csv(sampleIt, paste0(folSim, idSimulation, "/itersPerSample.csv"))

if(isDiffAnalysis) {
  events = as.data.frame(cbind(nl1, nc1, nl2, nc2))
  colnames(events) = c(paste0("nl1.", includedSamples), 
                       paste0("nc1.", includedSamples), 
                       paste0("nl2.", includedSamples), 
                       paste0("nc2.", includedSamples))
  itl2 = as.data.frame(do.call("rbind", lapply(itl2, function(x) do.call("rbind", x))))
  names(itl2) = idsGen
  write.csv(itl2, paste0(folSim, idSimulation, "/intInTarget2.csv"))
} else {
  events = as.data.frame(cbind(nl1, nc1))
  colnames(events) = c(paste0("nl1.", includedSamples), 
                       paste0("nc1.", includedSamples))
}
itl1 = as.data.frame(do.call("rbind", lapply(itl1, function(x) do.call("rbind", x))))
names(itl1) = idsGen
write.csv(itl1, paste0(folSim, idSimulation, "/intInTarget1.csv"))
write.csv(events, paste0(folSim, idSimulation, "/eventsInTarget.csv"))

matTimes = do.call("rbind", timesSim)
colnames(matTimes) = paste0("pt", 0:3)
write.csv(matTimes, paste0(folSim, idSimulation, "/times.csv"))


