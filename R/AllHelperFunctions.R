#' @title warpDEDataSet constructor
#' @name warpDEDataSet
#' @description user friendly constructor for warpDEDataSet
#' @param data.frame expression matrix
#' @param matrix pseudotimes values
#' @param matrix weights values
#' @export
warpDEDataSet <- function(counts, t, w){
       cat ("~~~~~ warpDEDataSet: constructor ~~~~~ \n")
       new (Class="warpDEDataSet", counts = log1p(counts), t = t, w = w)
  }

#' @describeIn warpDEDataSet returns the matrix of counts.
#' @param x a \code{warpDEDataSet} object.
#' @export
setMethod(
  f = "counts",
  signature = "warpDEDataSet",
  definition = function(x) x@counts
)

#' @describeIn warpDEDataSet returns the vector of pseudotimes..
#' @export
setMethod(
  f = "t",
  signature = "warpDEDataSet",
  definition = function(x) x@t
)

#' @describeIn warpDEDataSet returns the vector of weights.
#' @export
setMethod(
  f = "w",
  signature = "warpDEDataSet",
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
  definition = function(x) data.frame(x@ranking.df)
)

#' @describeIn rankingDE returns list of parameters used to create the ranking.
#' @export
setMethod(
  f = "params",
  signature = "rankingDE",
  definition = function(x) x@params
)


#' @describeIn rankingDE returns the ranks and criteria of some genes.
#' @export
setMethod(
  f = "g.rank",
  signature = "rankingDE",
  definition = function(x, g.list) x@ranking.df[g.list,]
)
