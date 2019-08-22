#' Estimating a block model for sequential data
#'
#' @description \code{block} is used to estimate a 'block model' (REFERENCE)
#' based on the dominant level of a variable in a given ordered sequence.
#' @param formula an object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' The details of model specification are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible by
#'  as.data.frame to a data frame) containing the variables in the model.
#'  If not found in data, the variables are taken from environment(formula),
#'  typically the environment from which glm is called.
#' @param thr Threshold for the length of a block in order to be identified.
#' @param verbose a logical indicating if some "progress report" should be given.
#' @param na.action a function which indicates what should happen when the data contain NAs.
#'  The default is set by to na.pass.
#' @param method The method used for reduction of groups. See 'Details'.
#' @param ... Not currently used, just in case.
#' @details We define an "unloading block" as a group of fish that were
#'  predominantly of the same type (e.g. species, quality, size class) within
#'  each individual unloading. Unloading blocks are allowed to contain some minor
#'  amount of fish with different characteristics for practical reasons.
#'  To estimate the unloading blocks in each super-sample, the n-running
#'  proportion for each type of fish in the unloading is calculated and the
#'  dominant type (more than 50%) identified for each group of n consecutive
#'  fishes. An unloading block is then computed as the union of contiguous
#'  groups with the same dominant type. Small blocks (less than n fish unloaded)
#'  dividing two blocks of the same dominant type are absorved to generate
#'  uninterrupted unloading blocks of the same type. Since the expected length
#'  of a unloading block is expected to change as function of the total
#'  unloading length (i.e. the longer the unloading the longer the unloading
#'  blocks), the value n to compute the running proportions should be taken as
#'  a fraction (e.g. 3\%, thr=0.03) of the total number fish unloaded N.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data_file = system.file("data/supersamples_demo.csv", package="sharkbox")
#' dat = read.csv(data_file)
#' mod_block0 = blocks(order|trip ~ species, data=dat)
#' mod_block1 = blocks(order|trip ~ group, data=dat)
#' }
blocks = function(formula, data, thr=0.03, verbose=TRUE, na.action,
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

  fac = .getFactorLevels(mf=mf, vars=all.vars(gp$pred.formula), ind=ind)

  dat = split(fac, f = sf)
  thisBlocks = levels(fac)[unsplit(lapply(dat, .findBlocks, thr=thr), f=sf)]

  # tb = fac$levels[unlist(tapply(fac$f, INDEX = sf, FUN=.findBlocks, thr=thr))]

  check = .checkBlocks(thisBlocks, fac, thr, method)
  while(!attr(check, "ok")) {
    thisBlocks[thisBlocks==attr(check, "old")] = attr(check, "new")
    if(isTRUE(verbose))
      message(sprintf("Merging '%s' with '%s'", attr(check, "old"), attr(check, "new")))
    check = .checkBlocks(thisBlocks, fac, thr, method)
  }

  thisBlocks = thisBlocks[order(ind)] # back to original order

  bb = table(thisBlocks)

  output = list(coefficients = names(bb), residuals = NULL,
                fitted.values=thisBlocks, model=NULL, na.action=na.action,
                call=call, formula=formula, terms, data=mf, control=NULL,
                method=method)

  class(output) = "block"

  return(output)
}



# Internals ---------------------------------------------------------------

.sort = function(xf, gp) {
  if(!is.null(gp$split.names))
    return(order(xf[, gp$split.names], xf[, gp$response]))
  return(order(xf[, gp$response]))
}

.getFactorLevels = function(mf, vars, ind=NULL) {
  if(is.null(vars)) stop("vars cannot be NULL.")
  mf = mf[, vars, drop=FALSE]
  fac = factor(apply(mf, 1, paste, collapse="_"))
  cod = as.numeric(fac)
  if(!is.null(ind)) cod = cod[ind]
  levels(cod) = levels(fac)
  return(cod)
}


#' Internal for identify blocks in a sequence
#'
#' @param x A vector of integers
#' @param thr The threshold to identify blocks, as fraction of total length
#'
#' @return A vector with the blocks labeled by the dominant integer in the block
#'
#' @examples
.findBlocks = function(x, thr=0.03) {

  S = max(x)
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
  if(sum(!is.na(out2))<n) {
    warning("No blocks identified, returning NA")
    return(out2)
  }
  out3 = round(approx(x=out2, xout = seq_along(out2), rule=2)$y, 0)

  spx = table(out3, x)
  sp0 = setNames(apply(spx, 1, which.max), nm=NULL)
  gg0 = setNames(seq_along(unique(sp0)), nm = unique(sp0))
  out4 = setNames(gg0[as.character(sp0[out3])], nm=NULL)
  out4 = .removeBumps(out4, thr=thr)

  mainSp = apply(table(out4, factor(x, levels=spp)), 1, which.max)

  out = setNames(mainSp[as.character(out4)], nm=NULL)

  return(out)

}

# check for small blocks and suggest merging
.checkBlocks = function(thisBlocks, fac, thr, method) {
  bb = rle(thisBlocks)
  class(bb) = "list"
  bb = data.frame(bb)
  bx = rowsum(bb$lengths, bb$values)
  bx = data.frame(count=bx, blocks=as.numeric(table(bb$values)[rownames(bx)]))
  bx$keep = bx$count > thr*sum(bx$count)
  bx = bx[order(bx$keep, bx$blocks, bx$count), ]
  bx$name = rownames(bx)
  rownames(bx) = NULL
  check = bx
  attr(check, "ok") = TRUE
  # new -> old
  if(all(check$keep)) return(check)
  x = as.matrix(format(unclass(table(fac, thisBlocks))))
  class(x) = "numeric"
  rownames(x) = levels(fac)
  # Method especific --------------------------------------------------------
  x3 = t(t(x)/colSums(x))
  x4 = x3
  old = check[1, "name"]
  for(j in seq_len(nrow(x3))) x4[j, ] = pmin(x3[j,], x3[j, old])
  check$x4 = colSums(x4)[check$name]
  ind = order(check$x4, check$blocks, decreasing = TRUE)
  new = check[ind[2], "name"]
  # end of method
  attr(check, "old") = old
  attr(check, "new") = new
  attr(check, "ok") = FALSE
  return(check)
}

