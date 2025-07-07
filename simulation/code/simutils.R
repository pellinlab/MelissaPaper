setNumberDigits = function(numberSamples) {
  q = numberSamples; i = 1L
  while(q > 0L) {q = numberSamples %/% i; i = i * 10L}
  return(as.integer(log(i - 1L, 10)))
}

targetGenes = function(wh) {
  all = as.integer(c(22300000, 22320000, 37500000, 37510000, 54900000, 54980000, 97000000,	97040000, 100000000, 100005000))
  if(missing(wh)) return(all)
  else return(all[2*(wh-1)+(1:2)])
}

whInTarget = function(int, whGene) which(int >= targetGenes(whGene)[1] & int < targetGenes(whGene)[2])


