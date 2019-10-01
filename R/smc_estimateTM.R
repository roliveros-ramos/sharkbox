
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
  fac_x  = .getFactorLevels(mf=mf, vars=x$term)
  fac_by = rep(1, nrow(mf))
  levels(fac_by) = "1"
  x = data.frame(x=fac_x, by=fac_by)
  y = .estimateTM(x, INDEX=INDEX)
  diag(y) = 0
  JUMP = .norm(y) # jump matrix
  PROP = .norm(table(.num2Factor(fac_x)))

  TM = function(N) {
    A = JUMP
    diag(A) = -1
    A = A/(pmax(N*PROP, 2))
    diag(A) = diag(A) + 1
    A = A/rowSums(A)
    attr(A, "initial") = PROP
    return(A)
  }

  y = list(jump=JUMP, prop=PROP, TM=TM)

  class(y) = c("block.TM", class(y))

  return(y)
}

#' @export
estimateTM.mc.spec = function(x, mf, INDEX) {

  fac_x  = .getFactorLevels(mf=mf, vars=x$term)
  fac_by = .getFactorLevels(mf=mf, vars=x$by)
  by.var = x$by
  x = data.frame(x=fac_x, by=fac_by)
  y = .estimateTM(x, INDEX=INDEX)
  attr(y, "initial") = .norm(table(.num2Factor(fac_x)))
  attr(y, "by") = by.var
  class(y) = c("mc.TM", class(y))
  return(y)
}

#' @export
estimateTM.multinomial.spec = function(x, mf, INDEX) {

  fac_x  = .getFactorLevels(mf=mf, vars=x$term, asFactor = TRUE)
  fac_by = .getFactorLevels(mf=mf, vars=x$by, asFactor = TRUE)
  by.var = x$by
  y = table('X(t)'=fac_x, by=fac_by)
  y = t(t(y)/colSums(y, na.rm=TRUE))

  attr(y, "by") = by.var
  class(y) = c("multinomial.TM", class(y))
  return(y)
}

estimateTM.default = function(x, ...) {}

# Internal ----------------------------------------------------------------

.getFactorLevels = function(mf, vars, ind=NULL, asFactor=FALSE, asCharacter=FALSE) {
  if(isTRUE(asCharacter)) asFactor = FALSE
  if(is.null(vars) | identical(vars, "NA")) return(rep(1, length=nrow(mf)))
  mf = mf[, vars, drop=FALSE]
  fac = factor(apply(mf, 1, paste, collapse="_"))
  cod = as.numeric(fac)
  if(!is.null(ind)) cod = cod[ind]
  levels(cod) = levels(fac)
  if(isTRUE(asFactor)) cod = .num2Factor(cod)
  if(isTRUE(asCharacter)) cod = as.character(.num2Factor(cod))
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

.estimateTM = function(x, INDEX, ...) {
  x0 = tapply(seq_len(nrow(x)), INDEX = INDEX, FUN = .getSimpleTM, x=x)
  out = 0
  for(i in seq_along(x0)) out = out + x0[[i]]
  attr(out, "n") = .count(out) # total transitions
  out = .norm(out)
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

# Normalization -----------------------------------------------------------

.normMatrix = function(x) {
  s = rowSums(x, na.rm=TRUE)
  if(all(s==0)) return(x)
  s[s==0] = 1
  y = x/s
  attributes(y) = attributes(x)
  return(y)
}

.normArray = function(x) {
  y = apply(x, 3, .normMatrix)
  dim(y) = dim(x)
  attributes(y) = attributes(x)
  return(y)
}

.normVector = function(x) {
  s = sum(x, na.rm=TRUE)
  if(s==0) return(x)
  y = x/s
  dim(y) = NULL
  names(y) = names(x)
  return(y)
}

.norm = function(x) {
  if(length(dim(x)) == 1) return(.normVector(x))
  if(length(dim(x)) == 2) return(.normMatrix(x))
  if(length(dim(x)) == 3) return(.normArray(x))
  return(.normVector(x))
}

.count = function(x) {
  if(is.matrix(x)) return(sum(x, na.rm=TRUE))
  if(is.array(x)) {
    out = apply(x, 3, sum, na.rm=TRUE)
    names(out) = dimnames(x)[[3]]
  return(out)
  }
  return(invisible(NULL))
}

