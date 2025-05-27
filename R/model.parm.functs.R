#' model.frame method for objects of class 'gls'
#'
#' This is a specific \code{model.frame} method for objects of the class
#' \code{gls}. It simply calls \code{nlme::getData}. See \link[nlme]{getData}.
#'
#' @param formula An object of class \code{gls}.
#' @param ... Further arguments passed to \code{model.frame.formula}.
#' @importFrom stats coef model.frame model.matrix terms
#' @importFrom stats vcov
#' @method model.frame gls
#' @export
model.frame.gls <- function(formula, ...) {
  model.frame(formula(formula), nlme::getData(formula),...)
}

#' model.frame method for objects of class 'lme'
#'
#' This is a specific \code{model.frame} method for objects of the class
#' \code{lme}.
#'
#' @param formula An object of class \code{lme}.
#' @param ... Further arguments passed to \code{model.frame.formula}.
#' @importFrom stats coef model.frame model.matrix terms
#' @importFrom stats vcov
#' @method model.frame lme
#' @export
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
