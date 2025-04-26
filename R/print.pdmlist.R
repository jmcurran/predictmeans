#' print method for objects of class 'pdmlist'
#'
#' This is a specific \code{print} method for objects of the class
#' \code{pdmlist}. See \code{\link{predictmeans}}.
#'
#' @param x An object of class \code{pdmlist}.
#' @param ... Additional arguments passed to \code{print}.
#' @method print pdmlist
#' @export
print.pdmlist <- function(x, ...) {
  pos <- grep(
    "predictmeansPlot|predictmeansBKPlot|predictmeansBarPlot|p_valueMatrix",
    names(x)
  )
  x <- x[names(x)[-pos]]
  NextMethod()
}
