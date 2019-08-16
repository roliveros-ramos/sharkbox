
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

    .simulator = function(N) {
      out = matrix(nrow=N, ncol=3)
      out[, 1] = .simMC(N, object$groups$TM(N))
      out[, 2] = .simSMC(out[, 1], object$species)
      out[, 3] = .simSMC(out[, 2], object$size)
      colnames(out) = c("group", "species", "size_class")
      class(out) = c("smc_traj", class(out))
      return(out)
    }

    if(length(N)==1) return(.simulator(N))

    out = vector("list", length = length(N))
    for(i in seq_along(N)) out[[i]] = .simulator(N=N[i])

  }

  if(is.null(newdata)) return(simulator)

}



# smc_traj class ----------------------------------------------------------

#' Title
#'
#' @param x
#' @param col
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.smc_traj = function(x, col=NULL, ...) {

  if(is.null(col)) col = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))[-c(6,11)]

  ll = c(0, cumsum(rle(x[,1])$lengths))
  plot.new()
  plot.window(xlim=c(0,N), ylim=c(-0.2,3), xaxs="i", yaxs="i")
  ylim = c(-0.2, -0.2, 0, 0)
  xcol = as.numeric(as.factor(x[,2]))
  points(x[, 3], type="h", col=xcol, pch=19, cex=0.4)
  for(i in seq_len(length(ll)-1)) {
    xx = ll[i + c(0,1,1,0)]
    polygon(x=xx, y=ylim, col=ifelse(i%%2==1, "gray90", "gray60"), border=1)
  }

  return(invisible())

}


