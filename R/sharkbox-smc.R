



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
    A = A/(N*prop)
    diag(A) = diag(A) + 1
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
  for(i in seq_len(S)) x$size[[i]] =
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

