
#' Exact Multinomial Test
#'
#' @param observed Vector (or matrix) of observed values. The length
#' (or number of columns) defines the number of categories of the multinomial
#' distribution. The sum of the vector (or rows) the sample size.
#' @param prob The theoretical probability of the multinomial distribution.
#'
#' @return
#' @export
#'
#' @examples
ExactMultinomialTest = function (observed, prob) {

  .pvalue = function(p, mass) sum(mass[mass <= p])

  if(is.matrix(observed)) {

    groups = ncol(observed)
    sizes  = rowSums(observed)
    size   = unique(sizes)
    if(length(size)==1) {
      pObs = apply(observed, 1, function(x) dmultinom(x, size = size, prob = prob))
      all_outcomes = multinomial_support(groups, size)
      mass = apply(all_outcomes, 1, function(x) dmultinom(x, size = size, prob = prob))
      p.value = sapply(pObs, FUN = .pvalue, mass=mass)
    } else {
      stop("All observed samples must have the same size.")
      pObs = p.value = rep(NA, nrow(observed))
      for(iSize in size) {
        ind = which(sizes==iSize)
        pObs[ind] = apply(observed[ind, ], 1, function(x) dmultinom(x, size = iSize, prob = prob))
        all_outcomes = multinomial_support(groups, iSize)
        mass = apply(all_outcomes, 1, function(x) dmultinom(x, size = iSize, prob = prob))
        p.value[ind] = sapply(pObs[ind], FUN = .pvalue, mass=mass)
      }
    }

  } else {

    groups = length(observed)
    size   = sum(observed)
    pObs = dmultinom(observed, size = size, prob)
    all_outcomes = multinomial_support(groups, size)
    mass = apply(all_outcomes, 1, function(x) dmultinom(x, size = size, prob = prob))
    p.value = .pvalue(pObs, mass)

  }

  output = list(statistic = pObs, parameter = list(size = size, groups = groups),
                p.value = p.value, method="Exact Multinomial Test",
                observed = observed, expected = size*prob, residuals = NULL,
                pmf = sort(mass, decreasing = TRUE), n=nrow(all_outcomes))
  class(output) = "EMtest"

  return(output)

}

#' @export
print.EMtest = function(x, ...) {

  tab = data.frame(events=x$n, statistic=round(x$statistic, digits = 4),
                   p.value = round(x$p.value, digits = 4))
  cat("\n Exact Multinomial Test\n\n")
  print(tab, row.names = FALSE)
  return(invisible(tab))

}



# Internal ----------------------------------------------------------------


#' Support of a multinomial distribution
#'
#' @param k The number of categories.
#' @param n The sample size.
#'
#' @return
#' @export
#'
#' @examples
multinomial_support = function(k, n) {
  if(k == 1) return(n)
  size = choose(n+k-1,k-1)
  mat = matrix(0, nrow=size, ncol=k)
  j = 1
  for(i in seq_len(n)) {
    p = choose(i+k-2, k-2)
    mat[j + 1:p, seq_len(k-1)] = multinomial_support(k - 1, i)
    j = j + p
  }
  mat[, k] = n - rowSums(mat)
  return(mat)
}
