
#' Estimate transition matrix for SMC model components
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
estimateTM = function(x, ...) {
  UseMethod("estimateTM")
}

#' @export
estimateTM.block.spec = function(x, mf, INDEX) {
  fac_x  = .getFactorLevels(mf=mf, vars=x$term, ind=ind)
  fac_by = rep(1, nrow(mf))
  levels(fac_by) = "1"
  x = data.frame(x=fac_x, by=fac_by)
  y = .estimateTM(x, INDEX=INDEX)
  class(y) = c("block.TM", class(y))
  return(y)
}

#' @export
estimateTM.mc.spec = function(x, mf, INDEX) {
  fac_x  = .getFactorLevels(mf=mf, vars=x$term, ind=ind)
  fac_by = .getFactorLevels(mf=mf, vars=x$by, ind=ind)
  x = data.frame(x=fac_x, by=fac_by)
  y = .estimateTM(x, INDEX=INDEX)
  class(y) = c("mc.TM", class(y))
  return(y)
}

#' @export
estimateTM.multinomial.spec = function(x, mf, INDEX) {
  fac_x  = .getFactorLevels(mf=mf, vars=x$term, ind=ind, asFactor = TRUE)
  fac_by = .getFactorLevels(mf=mf, vars=x$by, ind=ind, asFactor = TRUE)
  y = table('X(t)'=fac_x, by=fac_by)
  class(y) = c("multinomial.TM", class(y))
  return(y)
}

estimateTM.default = function(x, ...) {}

# Internal ----------------------------------------------------------------

.getFactorLevels = function(mf, vars, ind=NULL, asFactor=FALSE) {
  if(is.null(vars)) return(rep(1, length=nrow(mf)))
  mf = mf[, vars, drop=FALSE]
  fac = factor(apply(mf, 1, paste, collapse="_"))
  cod = as.numeric(fac)
  if(!is.null(ind)) cod = cod[ind]
  levels(cod) = levels(fac)
  if(isTRUE(asFactor)) cod = .num2Factor(cod)
  return(cod)
}

.num2Factor = function(x, .levels=NULL) {
  if(!is.numeric(x)) stop("x must be numeric!")
  lev = if(is.null(.levels)) levels(x) else .levels
  if(is.null(lev)) lev = seq_len(max(x, na.rm=TRUE))
  x = factor(x, levels=seq_along(lev))
  levels(x) = lev
  return(x)
}

.getSimpleTM = function(ind, x) {

  if(!is.data.frame(x)) stop("Argument must be a data.frame!")
  checkVars = all(c("x", "by") %in% names(x))
  if(!checkVars) stop("Incorrect 'data.frame' format.")
  by = .num2Factor(x$by[ind], .levels=levels(x$by))
  x  = .num2Factor(x$x[ind], .levels=levels(x$x))
  out = table('X(t)'=head(x, -1), 'X(t+1)'=tail(x, -1), by=head(by, -1))
  out = drop(out)

  return(out)

}

.estimateTM = function(x, INDEX, diagonal2zero=FALSE, ...) {
  x0 = tapply(seq_len(nrow(x)), INDEX = INDEX, FUN = .getSimpleTM, x=x)
  out = 0
  for(i in seq_along(x0)) out = out + x0[[i]]
  if(isTRUE(diagonal2zero)) diag(out) = 0
  return(out)
}

.simplify = function(out) {
  cs = colSums(out)
  rs = rowSums(out)
  ind = which((cs==0) & (rs==0))
  if(length(ind)>0) {
    out = out[-ind, ]
    out = out[, -ind]
  }
  return(out)
}

.normMatrix = function(x) {
  s = rowSums(x, na.rm=TRUE)
  if(all(s==0)) return(x)
  s[s==0] = 1
  return(x/s)
}

.norm = function(x) {
  s = sum(x, na.rm=TRUE)
  if(s==0) return(x)
  return(x/s)
}
