control <- do.call("gam.control", control)

if (nrow(mf) < 2)
  stop("Not enough (non-NA) data to do anything meaningful")
terms <- attr(mf, "terms")
vars <- all.vars1(gp$fake.formula[-2])
inp <- parse(text = paste("list(", paste(vars,
                                         collapse = ","), ")"))
if (!is.list(data) && !is.data.frame(data))
  data <- as.data.frame(data)
dl <- eval(inp, data, parent.frame())
names(dl) <- vars

xsmc = function(formula, data, ...) {
  gp <- interpret.smc(formula)
  # check only one block term is included
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- gp$fake.formula
  mf[[1]] <- quote(stats::model.frame)
  data$block = NA
  mf <- eval(mf, data, parent.frame())
  # split dataset
  # check each split has complete data


  return(mf)
}

x = xsmc(order ~ b(species) + mc(size, by=block), data=dat0)

y = xsmc(order|trip ~ b(species) + mc(size), data=dat0)


  call <- match.call()

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset",
               "weights", "na.action", "etastart",
               "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
