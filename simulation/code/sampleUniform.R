# n = 5L 
# cp = c(2L, 17L, 23L); prob = c(2,6)
# replace = FALSE; sorted = FALSE
rDiscreteUnifMix = function(n, cp, prob, right = FALSE, replace = TRUE, sorted = FALSE) {
  stopifwrong1 = function(cp, prob) {
    stopifnot(is.vector(cp))
    stopifnot(length(cp) >= 2L)
    stopifnot(is.numeric(cp))
    stopifnot(all(cp == as.integer(cp)))
    stopifnot(all(cp == sort(cp)))
    if(length(cp) > 2L) {
      stopifnot(is.vector(prob))
      stopifnot(is.numeric(prob))
      stopifnot(length(prob) == (length(cp) - 1L))
    } 
    return(invisible(NULL))
  }
  stopifwrong2 = function(replace, sorted) return(invisible(NULL))
  stopifwrong1(cp, prob)
  stopifwrong2(replace, sorted)
  
  if(length(n) == 0L) return(integer())
  else n = as.integer(n[1L])
  if(right) cp = cp + 1L
  
  if(length(cp) > 2L) {
    min = cp[-length(cp)]
    max = cp[-1L] #+ c(rep(0L, length(min) - 1L), 1L)
    probBlocks = prob * (max - min)
    tb = table(sample.int(length(prob), n, replace = TRUE, prob = probBlocks))
    blocks = as.integer(names(tb))
    nsamples = as.integer(tb)
    rl = vector("list", length(blocks))
    #rl = lapply(1:length(rl), function(i) sample.int(n = max[blocks[i]] - min[blocks[i]], size = nsamples[i], replace = replace) - 1L + min[blocks[i]])
    for(i in 1:length(rl)) rl[[i]] = sample.int(n = max[blocks[i]] - min[blocks[i]], size = nsamples[i], replace = replace) - 1L + min[blocks[i]]
    if(sorted) rl = lapply(rl, sort)
    r = unlist(rl)
  }
  else {
    r = sample.int(n = cp[2] + 1L - cp[1], size = n, replace = replace) - 1L + cp[1]
    if(sorted) r = sort(r)
  }
  
  return(r)
}

dDiscreteUnifMix = function(x, cp, prob, right = FALSE, log = FALSE) {
  stopifwrong1 = function(cp, prob) return(invisible(NULL))
  stopifnot(is.numeric(x))
  stopifwrong1(cp, prob)
  
  if(length(x) == 0L) return(numeric())
  if(right) cp = cp + 1L
  
  p = prob * (cp[-1L] - cp[-length(cp)])
  p = prob / sum(p)
  f = x == as.integer(x)
  g = list(base::identity, base::log)[[1L + log]]
  
  return(g(sapply(x, function(i) c(0, p, 0)[(1L + sum(i >= cp))]) * f))
}
