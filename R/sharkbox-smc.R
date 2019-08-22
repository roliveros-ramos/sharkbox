
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

  fac = .getFactorLevels(mf=mf, vars=gp$smooth.spec[[1]]$term, ind=ind)
  x0 = estimateTM(fac, INDEX=sf, diagonal2zero=TRUE)
  x0 = .normMatrix(x0)

  output = list(coefficients = x0, residuals = NULL,
                fitted.values=NULL, model=NULL, na.action=na.action,
                call=call, formula=formula, terms, data=mf, control=NULL,
                method=method)

  class(output) = "smc"

  return(output)
}


.countTraj = function(x, S, L=3, thr=0.03) {
  x0 = tapply(x[,1], INDEX = y, FUN = .getTM, S=S)
  x1 = tapply(x[,2], INDEX = x[,1], FUN = .getTM, S=L)
  mainSp = apply(table(y, factor(x[,1], levels=seq_len(S))), 1, which.max)
  nn = newSMC(G=S, S=S, L=L)
  nn$species[mainSp] = x0
  nn$groups$prop[mainSp] = table(y)
  nn$groups$jump[mainSp, mainSp] = .getTM(x=y, simplify=TRUE)
  diag(nn$groups$jump) = 0
  nn$size[as.numeric(names(x1))] = x1
  return(nn)
}

