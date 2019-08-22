
#' Title
#'
#' @param data
#' @param S
#' @param L
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
xsmc = function(data, S, L, control=NULL, ...) {

  if(is.null(control)) control=list(thr=0.03, min_count=10)

  x = lapply(data, FUN=.countTraj_old, S=S, L=L, thr=control$thr) # check defaults
  out = newSMC(G=S, S=S, L=L)

  for(i in seq_along(x)) {
    out$groups$prop = out$groups$prop + x[[i]]$groups$prop
    out$groups$jump = out$groups$jump + x[[i]]$groups$jump
    for(j in seq_along(out$species)) out$species[[j]] = out$species[[j]] + x[[i]]$species[[j]]
    for(k in seq_along(out$size)) out$size[[k]] = out$size[[k]] + x[[i]]$size[[k]]
  }

  # verify counts: if lower than NN, assume random? give a warning!
  out = .verifySMC(out, thr=control$min_count)

  out = .normSMC(out) # transform to probabilities

  kg = head(order(out$groups$prop, decreasing = TRUE), sum(out$groups$prop>control$thr))

  out$groups$prop = out$groups$prop[kg]
  out$groups$jump = out$groups$jump[kg, kg]
  out$species = out$species[kg]

  out = .normSMC(out) # re-normalize after removing non-key groups

  M = out$groups$jump
  prop = out$groups$prop

  out$groups$TM = function(N) {
    A = M
    diag(A) = -1
    A = A/(pmax(N*prop, 2))
    diag(A) = diag(A) + 1
    A = A/rowSums(A)
    return(A)
  }

  return(out)

}


#' Simulate a transition matrix for a Sequential Markov Chaing
#'
#' @param G
#' @param S
#' @param L
#' @param N
#'
#' @return
#' @export
#'
#' @examples
smcSim = function(G, S, L) {

  x = newSMC(G=G, S=S, L=L) # species greater than groups

  G = nrow(x$groups$jump)
  S = nrow(x$species[[1]])
  L = sapply(x$size, nrow)

  # Unloading group matrix
  M = x$groups$jump
  # diag(M) = -1
  M[row(M) == (col(M) - 1)] = 1 # can vary! rows sums one.
  prop = as.numeric(rbrokenbar(n=1, S=G))
  x$groups$jump = M
  x$groups$prop[] = prop

  x$groups$TM = function(N) {
    A = M
    diag(A) = -1
    A = A/(pmax(N*prop, 2))
    diag(A) = diag(A) + 1
    A = A/rowSums(A)
    return(A)
  }

  # Species matrix
  main_species = sample(x=S, size=G)

  .sortSpecies = function(first, size) {
    c(first, sample(setdiff(seq_len(size), first)))
  }

  .coverage = function(A) {
    out = sapply(seq_along(A), FUN=function(x) length(unique(A[seq_len(x)])))
    return(which.max(out))
  }

  A = t(sapply(main_species, .sortSpecies, size=S))
  A[tail(seq_along(A), -.coverage(A))] = NA
  A = unlist(apply(A, 1, FUN=function(x) list(c(na.omit(x)))), recursive = FALSE)

  for(i in seq_len(G)) x$species[[i]][A[[i]], A[[i]]] =
    rbrokenbar(n=length(A[[i]]), S=length(A[[i]]), sequential=TRUE)

  # Size matrix
  for(i in seq_len(S)) x$size[[i]][,] =
    rbrokenbar(n=nrow(x$size[[i]]), S=nrow(x$size[[i]]), shuffle = TRUE)

  return(x)

}


# Auxiliar functions ------------------------------------------------------

#' Title Steady states for a Markov chain
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
steadyStates = function(x) {

  .eigen = function (x, transpose = TRUE) {
    eigenResults = eigen(x = t(x), symmetric = FALSE)
    ind = which(round(eigenResults$values, 3) == 1)
    if (length(ind) == 0) {
      warning("No eigenvalue = 1 found - the embedded Markov Chain must be irreducible, recurrent")
      return(NULL)
    }
    if (length(ind) > 1) {
      warning("The Markov Chain is not irreducible nor recurrent.")
      return(NULL)
    }
    out = as.matrix(t(eigenResults$vectors[, ind]))
    if (rowSums(Im(out)) != 0) {
      warning("The Markov Chain is not irreducible nor recurrent.")
      return(NULL)
    }
    return(out)
  }

  out = .eigen(x=x, transpose=TRUE)
  if(is.null(out)) {
    warning("Warning! No steady state")
    return(NULL)
  }
  out = out / rowSums(out)
  out = suppressWarnings(as.double(out))

  return(out)
}

findLandingGroups = function(x, S=NULL, thr=0.03) {

  if(is.null(S)) S = max(x)
  spp = seq_len(S)

  n = max(10, thr*length(x))

  out = matrix(0, ncol=S, nrow=length(x))

  for(iSp in spp) {
    out[, iSp] = 0 + (x == iSp)
  }
  colnames(out) = seq_len(S)

  out0 = apply(out, 2, .getRates, n=n)
  out1 = apply(out, 2, .getBlocks, n=n)
  out2 = rle(rowSums(out1))
  ones = out2$values == 1
  out2$values[ones] =  seq_len(sum(ones))
  out2$values[!ones] = NA
  out2 = inverse.rle(out2)
  out3 = round(approx(x=out2, xout = seq_along(out2), rule=2)$y, 0)

  spp = table(out3, x)
  sp0 = apply(spp, 1, which.max)
  names(sp0) = NULL
  gg0 = setNames(seq_along(unique(sp0)), nm = unique(sp0))
  out4 = sp0[out3]
  out4 = gg0[as.character(out4)]
  names(out4) = NULL
  out4 = .removeBumps(out4, thr=thr)

  return(out4)

}
