

# Broken bar distribution -------------------------------------------------

#' Broken bar distribution
#'
#' @param n number of observations.
#' @param S number of parts to break the unit bar. See details.
#' @param sequential Boolean. If TRUE, the small pieces is sequentially broken. See details.
#'
#' @return
#' @export
#'
#' @examples
rbrokenbar = function(n, S, sequential=FALSE, shuffle=FALSE) {

  if(S<2) stop("S must be greater or equal than 2.")

  .rbrokenbar1 = function(n, S, shuffle) {
    .brokenbar = function(x) sort(diff(sort(x)), decreasing = TRUE)
    breaks = matrix(0, nrow=n, ncol=S+1)
    breaks[, -c(1,S+1)] = runif(n*S-n)
    breaks[, S+1] = 1
    numbers = t(apply(breaks, 1, .brokenbar))
    if(isTRUE(shuffle)) numbers = t(apply(numbers, 1, sample))
    return(numbers)
  }

  .rbrokenbar2 = function(n, S, shuffle) {
    num = array(dim=c(n, 2, S-1))
    for(i in seq_len(S-1)) num[,,i] = rbrokenbar(n, S=2)
    for(i in seq_len(S-2)) num[,,i+1] = num[,,i+1]*num[,2,i]
    numbers = matrix(nrow=n, ncol=S)
    numbers[, seq_len(S-1)] = num[,1,]
    numbers[, S] = num[,2,S-1]
    if(isTRUE(shuffle)) numbers = t(apply(numbers, 1, sample))
    return(numbers)
  }

  if(S==2) return(.rbrokenbar1(n, S, shuffle))

  if(isTRUE(sequential)) return(.rbrokenbar2(n, S, shuffle))

  return(.rbrokenbar1(n, S, shuffle))

}
