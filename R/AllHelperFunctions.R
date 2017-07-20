# accessor methods
#' @describeIn lineageDEDataSet returns the matrix of logcounts.
#' @param x a \code{lineageDEDataSet} object.
#' @export
setMethod(
  f = "logCounts",
  signature = "lineageDEDataSet",
  definition = function(x) x@logCounts
)

#' @describeIn lineageDEDataSet returns the vector of pseudotimes..
#' @export
setMethod(
  f = "t",
  signature = "lineageDEDataSet",
  definition = function(x) x@t
)

#' @describeIn lineageDEDataSet returns the vector of weights.
#' @export
setMethod(
  f = "w",
  signature = "lineageDEDataSet",
  definition = function(x) x@w
)
#'

# accessor methods
#' @describeIn rankingDE returns the ranking dataframe.
#' @param x a \code{rankingDE} object.
#' @export
setMethod(
  f = "ranking.df",
  signature = "rankingDE",
  definition = function(x) x@ranking.df
)

#' @describeIn rankingDE returns list of parameters used to create the ranking.
#' @export
setMethod(
  f = "params",
  signature = "rankingDE",
  definition = function(x) x@params
)
