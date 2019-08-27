
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

