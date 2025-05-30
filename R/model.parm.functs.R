
model.frame.gls <- function(formula, ...) {
  model.frame(formula(formula), nlme::getData(formula),...)
}

model.frame.lme <- function(formula, ...) {
  formula$data
}

model.matrix.aovlist <- function(object, ...) {
    stop(sQuote("predicted.means"), " does not support objects of class ",
         sQuote("aovlist"))
}

model.matrix.gls <- function(object, ...) {
    model.matrix(terms(object), data = nlme::getData(object), ...)
}

model.matrix.lme <- function(object, ...) {
    model.matrix(terms(object), data = model.frame(object), ...)
}

terms.gls <- function(x, ...) {
  terms(model.frame(x),...)
}

terms.lme <- function (x, ...) {
  v <- x$terms
  if (is.null(v)) {
    stop("no terms component")
  }
  return(v)
}

#terms.merMod <- function (object, ...) {
#  v <- terms(object)
#  attr(v, "dataClasses") <- sapply(all.vars(formula(object, fixed.only=TRUE)),function(x) class(model.frame(object)[[x]])[1])
#  return(v)
#}
#
