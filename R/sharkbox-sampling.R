
#' Sampling simulator
#'
#' @param object Population to be sampled, usually a data.frame
#' @param skip Number of individuals to skip before next observation.
#' The default is zero.
#' @param start Start of the sampling as a fraction of the total observations.
#' The default is zero (sampling start at the beginning).
#' @param fraction Fraction to be sampled (after start). The default is one,
#' all the population is sampled.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @return The sample from the population.
#' @export
#'
#' @examples
sampling = function(object, skip=0, start=0, fraction=1, ...) {
  UseMethod("sampling")
}

#' @export
sampling.data.frame = function(object, skip=0, start=0, fraction=1, ...) {
  n = nrow(object)
  ind = sampling(seq_len(n), skip=skip, start=start, fraction=fraction, ...)
  return(object[ind, ])
}

#' @export
sampling.default = function(object, skip=0, start=0, fraction=1, ...) {

  n = length(object)
  pos  = (seq_len(n)-1)/n
  ind = which(pos>=start & pos<(start+fraction))
  sub = c(TRUE, rep(FALSE, skip))
  ind = suppressWarnings(ind[sub])
  return(object[ind])

}

#' @export
sampling.list = function(object, skip=0, start=0, fraction=1, ...) {
  lapply(object, FUN=sampling, skip=skip, start=start, fraction=fraction, ...)
}
