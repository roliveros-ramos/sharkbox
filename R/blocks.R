

#' Blocks model
#'
#' @param formula
#' @param data
#' @param thr
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
blocks = function(formula, data, thr=0.03, verbose=TRUE, ...) {

  gp = interpret.smc(formula)
  # check only one block term is included
  cl = match.call()
  # frame to sort
  xf = match.call(expand.dots = FALSE)
  xf$formula = gp$sort.formula
  xf[[1]] = quote(stats::model.frame)
  xf$na.action = na.pass # check for NAs manually
  xf = eval(xf, data, parent.frame())
  ind = order(xf[, gp$split.names], xf[, gp$response])
  # frame for blocks
  mf = match.call(expand.dots = FALSE)
  mf$formula = gp$pred.formula
  mf[[1]] = quote(stats::model.frame)
  mf$na.action = na.pass # check for NAs manually
  mf = eval(mf, data, parent.frame())
  # split dataset
  sf <- match.call(expand.dots = FALSE)
  sf$formula = gp$split.formula
  sf[[1]] = quote(stats::model.frame)
  sf$na.action = na.pass # check for NAs manually
  sf = eval(sf, data, parent.frame())
  sf = sf[ind, ]
  # check each split has complete data
  fac = .getFactorLevels(mf, ind)

  dat = split(fac$f, f = sf)
  thisBlocks = fac$levels[unsplit(lapply(dat, .findBlocks, thr=thr), f=sf)]

  check = .checkBlocks(thisBlocks, fac, thr)
  while(!attr(check, "ok")) {
    thisBlocks[thisBlocks==attr(check, "old")] = attr(check, "new")
    if(isTRUE(verbose))
      message(sprintf("Merging '%s' with '%s'", attr(check, "old"), attr(check, "new")))
    check = .checkBlocks(thisBlocks, fac, thr)
  }

  thisBlocks = thisBlocks[order(ind)] # back to original order

  return(thisBlocks)
}



# Internals ---------------------------------------------------------------



.getFactorLevels = function(mf, ind=NULL) {
  fac = factor(apply(mf, 1, paste, collapse="_"))
  lev = levels(fac)
  cod = as.numeric(fac)
  if(!is.null(ind)) cod = cod[ind]
  return(list(f=cod, levels=lev))
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

.checkBlocks = function(thisBlocks, fac, thr) {
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
  x = as.matrix(format(unclass(table(fac$f, thisBlocks))))
  class(x) = "numeric"
  rownames(x) = fac$levels
  x3 = t(t(x)/colSums(x))
  x4 = x3
  old = check[1, "name"]
  for(j in seq_len(nrow(x3))) x4[j, ] = pmin(x3[j,], x3[j, old])
  check$x4 = colSums(x4)[check$name]
  ind = order(check$x4, check$blocks, decreasing = TRUE)
  new = check[ind[2], "name"]
  attr(check, "old") = old
  attr(check, "new") = new
  attr(check, "ok") = FALSE
  return(check)
}

