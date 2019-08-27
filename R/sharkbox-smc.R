
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

  smooths = sapply(gp$smooth.spec, class)

  out = list(mf=mf, sf=sf, gp=gp, ind=ind)
  return(out)

  output = list(coefficients = x0, residuals = NULL,
                fitted.values=NULL, model=NULL, na.action=na.action,
                call=call, formula=formula, terms=terms, data=mf, control=NULL,
                method=method)

  class(output) = "smc"

  return(output)
}

