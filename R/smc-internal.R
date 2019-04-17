
newSMC = function(G, S, L=3) {

  .generateNames = function(n, label) {
    if(is.character(n)) return(list(n, n))
    sp_code = sprintf("%s_%%0%dd", label, ceiling(log10(n)))
    species_names = sprintf(sp_code, seq_len(n))
    return(list(species_names, species_names))
  }

  if(length(L)==1) L = rep(L, S)
  if(length(L)!=S)
    stop("Number of size groups (L) must be equal to species number (S).")

  species_names = .generateNames(S, "species")
  group_names   = .generateNames(G, "group")

  if(is.character(S)) S = length(S)
  if(is.character(G)) G = length(G)

  out = list()

  out$groups = list()
  out$groups$prop = setNames(numeric(G), nm = group_names[[1]])
  out$groups$jump = matrix(0, nrow=G, ncol=G, dimnames = group_names)

  out$species = list()
  out$size    = list()

  for(i in seq_len(G)) {
    out$species[[i]] = matrix(0, nrow=S, ncol=S, dimnames = species_names)
  }

  names(out$species) = group_names[[1]]

  for(i in seq_len(S)) {
    size_names = .generateNames(L[i], "size")
    if(L[i]==2) size_names = .generateNames(c("small", "large"))
    if(L[i]==3) size_names = .generateNames(c("small", "medium", "large"))
    out$size[[i]]    = matrix(0, nrow=L[i], ncol=L[i], dimnames = size_names)
  }
  names(out$size) = species_names[[1]]

  class(out) = "smc"

  return(out)

}


# auxiliar ----------------------------------------------------------------

# simulates one markov chain trajectory
.simMC = function(n, A, initial=NULL) {

  out = numeric(n)
  initialNull = which(rowSums(A) != 0)[1]
  out[1] = if(is.null(initial)) initialNull else initial # deterministic now
  states = seq_len(nrow(A))
  for(i in seq_len(n-1)) {
    out[i+1] = sample(states, size=1, prob = A[out[i], ])
  }
  return(out)
}

# simulates one trajectory conditional to previous state
.simSMC = function(x, A) {
  Ns = table(x)
  groups = as.numeric(names(Ns))
  y1 = numeric(length(x))
  for(i in groups) {
    y1[x==i] = .simMC(Ns[as.character(i)], A[[i]])
  }
  return(y1)
}


# Estimation --------------------------------------------------------------

.completeTails = function(x) {
  nn = length(x)
  nona = which(!is.na(x))
  h = head(nona, 1)
  t = tail(nona, 1)
  x[seq_len(h-1)] = x[h]
  x[t+seq_len(nn-t)] = x[t]
  return(x)
}

.getBlocks = function(x, n) {
  x1 = 0 + (.getRates(x=x, n=n) > 0.5)
  x1 = .completeTails(x1)
  return(x1)
}

.getRates = function(x, n) {
  x0 = as.numeric(filter(x, filter=rep(1,n)/n))
  return(x0)
}


.removeBumps = function(x, thr) {

  x = rle(x)

  px = x$lengths/sum(x$lengths)
  x$values[px<thr] = NA
  validGroups = unique(na.omit(x$values))
  validGroups = setNames(seq_along(validGroups), nm=validGroups)
  validGroups = validGroups[as.character(x$values)]
  names(validGroups) = NULL
  x$values = validGroups
  x = inverse.rle(x)
  x = round(approx(x=x, xout = seq_along(x), rule=2)$y, 0)

  x = rle(x)
  if(length(x$values)<3) {
    x$values = seq_along(x$values)
    x = inverse.rle(x)
    return(x)
  }

  while(any(table(x$values)>1)) {
    bumps = numeric(length(x$values))
    for(i in seq_len(length(x$values)-2)) {
      # if(x$values[i]==x$values[i+2]) x$values[i+1] = x$values[i]
      if(x$values[i]==x$values[i+2]) bumps[i+1] = 1
    }
    thisFirst = which.min(x$lengths/bumps)
    x$values[thisFirst] = x$values[thisFirst+1]
    x = inverse.rle(x)
    x = rle(x)
  }
  # x = inverse.rle(x)
  # x = rle(x)
  x$values = seq_along(x$values)
  x = inverse.rle(x)

  return(x)
}


