
#' Title
#'
#' @param formula
#' @param data
#' @param thr
#' @param verbose
#' @param na.action
#' @param method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
smc = function(formula, data, thr=0.03, verbose=TRUE, na.action,
                  method="default", ...) {

  if(missing(na.action)) na.action = na.pass
  gp = interpret.smc(formula)
  # check only one block term is included
  cl = match.call()
  # full model frame
  mf = match.call(expand.dots = FALSE)
  mf$formula = gp$fake.formula
  mf[[1]] = quote(stats::model.frame)
  mf$na.action = na.pass # check for NAs manually
  mf = eval(mf, data, parent.frame())
  ind = .sort(mf, gp)
  sf = if(gp$doSplit) mf[ind , gp$split.names, drop=FALSE] else rep(1, length=nrow(mf))
  # check each split has complete data

  x0 = lapply(gp$smooth.spec, FUN=estimateTM, mf=mf[ind, , drop=FALSE], INDEX=sf)
  terms = lapply(gp$smooth.spec, FUN="[[", i="term")

  output = list(coefficients = x0, residuals = NULL,
                fitted.values=NULL, model=NULL, na.action=na.action,
                call=call, formula=formula, terms=terms, data=mf, control=NULL,
                method=method)

  class(output) = "smc"

  return(output)
}


# Methods -----------------------------------------------------------------


#' Title
#'
#' @param object
#' @param newdata
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.smc = function(object, newdata, ...) {

  if(missing(newdata)) newdata = NULL

  simulator = function(N) {

    par = coef(object)
    term = terms(object)
    exec = seq_along(par)

    .simulator = function(N) {
      mf = data.frame('__index__'=seq_len(N))
      for(i in seq_along(par)) {
        mf[, unlist(term[exec[i]])] = predict(par[[exec[i]]], mf=mf)
      }
      mf[, 1] = NULL
      return(mf)
    }

    if(length(N)==1) return(.simulator(N))

    out = vector("list", length = length(N))
    for(i in seq_along(N)) out[[i]] = .simulator(N=N[i])

    return(out)

  }

  if(is.null(newdata)) return(simulator)

}

