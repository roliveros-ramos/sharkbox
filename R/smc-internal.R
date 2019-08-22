
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
    if(sum(A[out[i],])==0) {
      print(i)
      print(out[i])
      print(rownames(A)[i])
      print(A[out[i],])
      print(out)
      print(A)
    }
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

.getSimpleTM = function(x, levels, simplify=FALSE, diagonal2zero=FALSE) {
  x = as.factor(x)
  levels(x) = levels
  out = table(head(x, -1), tail(x, -1))
  if(isTRUE(simplify)) {
    cs = colSums(out)
    rs = rowSums(out)
    ind = which((cs==0) & (rs==0))
    if(length(ind)>0) {
      out = out[-ind, ]
      out = out[, -ind]
    }
  }
  if(isTRUE(diagonal2zero)) diag(out) = 0
  return(out)
}

estimateTM = function(x, INDEX, ...) {
  x0 = tapply(x, INDEX = INDEX, FUN = .getSimpleTM, levels=levels(x), ...)
  out = 0
  for(i in seq_along(x0)) out = out + x0[[i]]
  return(out)
}


.norm = function(x) {
  s = sum(x, na.rm=TRUE)
  if(s==0) return(x)
  return(x/s)
}

.normMatrix = function(x) {
  s = rowSums(x, na.rm=TRUE)
  if(all(s==0)) return(x)
  s[s==0] = 1
  return(x/s)
}

.normSMC = function(x) {

  x$groups$prop = .norm(x$groups$prop)
  x$groups$jump = .normMatrix(x$groups$jump)
  for(j in seq_along(x$species)) x$species[[j]] = .normMatrix(x$species[[j]])
  for(k in seq_along(x$size)) x$size[[k]] = .normMatrix(x$size[[k]])
  return(x)

}

.verifySMC = function(x, thr) {

  for(j in seq_along(x$species)) {

    cs = colSums(x$species[[j]])
    rs = rowSums(x$species[[j]])
    ind0 = which((cs == 0) & (rs == 0))
    ind1 = which((cs > 0) & (rs == 0))

    x$species[[j]][ind1, ] = 1
    x$species[[j]][, ind0] = 0

  }

  for(k in seq_along(x$size)) {
    n = sum(x$size[[k]], na.rm=TRUE)
    if(n < thr) {
      warning("Insuficient data to estimate size-class transition matrix, using uniform prior.")
      x$size[[k]] = x$size[[k]] + 1
    }
  }

  return(x)

}

