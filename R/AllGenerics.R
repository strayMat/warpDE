# data object
#' @title Initialize an object of class \code{warpDEDataSet}
#' @name newwarpDEDataSet
#' @docType methods
#'
#' @description Constructs a \code{warpDEDataSet} object. Additional helper
#'   methods for manipulating \code{warpDEDataSet} objects are  also
#'   described below.
#'
#' @param counts matrix. An \code{n} by \code{m} numeric matrix or data
#'   frame giving the filtered, normalized RNA counts. In rows are the genes and in columns are the cells. It should contain the names of the cells and the genes in \code{colnames(counts)} and \code{rownames(counts)}.
#'
#' @param t matrix. An  \code{m} by 2 numeric matrix
#'   denoting each cell's reconstructed ordering position or pseudotime for the two lineages to compare.
#'
#' @param w matrix. An \code{m} by 2 numeric matrix
#'   denoting each cell's probability to belonng to each lineage. The probability should be convex, i.e The sum of the columns shoudl add up to 1.

setGeneric(
  name = "newwarpDEDataSet",
  signature = c('counts','t', 'w'),
  def = function(counts, t, w, ...) {
    standardGeneric("newwarpDEDataSet")
  }
)

# accessor functions
#' @title Returns the counts matrix of the dataset
#'@name counts
#' @param x an object that describes a dataset
#' @return the matrix representing log1p of the data.
#' @export
setGeneric(name = "counts",
           signature = "x",
           def = function(x) standardGeneric("counts"))

#' @title Returns the pseudotimes
#' @name t
#'
#' @description Extract the pseudotimes from a \code{SlingshotDataSet}.
#'
#' @param x an object that describes a dataset.
#' @return the matrix of pseudotimes.
#' @export
setGeneric(name = "t",
           signature = "x",
           def = function(x) standardGeneric("t"))

#' @title Returns the pseudotimes
#' @name w
#'
#' @description Extract the weights (probabilities to beloing to a lineage for a cell) from a \code{SlingshotDataSet}.
#'
#' @param x an object that describes a dataset.
#' @return the matrix of weights.
#' @export
setGeneric(name = "w",
           signature = "x",
           def = function(x) standardGeneric("w"))



# Ranking object
#' @title Initialize an object of class \code{rankingDE}
#' @name newrankingDE
#' @docType methods
#'
#' @description Constructs a \code{rankingDE} object. Additional helper
#'   methods for manipulating \code{rankingDE} objects are  also
#'   described below.
#'
#' @param ranking.df data.frame. An \code{n} by 2 numeric data.frame with the first column "dist" containing the per gene distance computed for the rankinga and the second column "rank" displaying the rank of the gene.
#'
#' @param params list. A list giving details about how the ranking has been obtained: which method, what distance has been used, ...
#'
setGeneric(
  name = "newrankingDE",
  signature = c('ranking.df','params'),
  def = function(ranking.df, params, ...) {
    standardGeneric("newrankingDE")
  }
)

# accessor functions
#' @title Returns the ranking
#' @name ranking.df
#' @param x a \code{rankingDE} object.
#' @return the matrix representing log1p of the data.
#' @export
setGeneric(name = "ranking.df",
           signature = "x",
           def = function(x) standardGeneric("ranking.df"))

#' @title Returns the parameters
#' @name params
#' @param x a \code{rankingDE} object.
#' @return the list of parameters.
#' @export
setGeneric(name = "params",
           signature = "x",
           def = function(x) standardGeneric("params"))


#' @title Returns the rank and distance of some genes
#' @name g.rank
#' @param x a \code{rankingDE} object.
#' @return the ranks and criteria values for the genes of interest.
#' @export
setGeneric(name = "g.rank",
           signature = "x",
           def = function(x, g.list) standardGeneric("g.rank"))
