#' @title Class \code{warpDEDataSet}
#' @aliases warpDEDataSet
#'
#' @description The \code{warpDEDataSet} class holds data relevant for
#'   performing gene differential expression analysis with the \code{warpDE} package, primarily a filtered, normalized count matrix of the data and a vecotr of pseudotimes for the cells as well as a vector of the weigths indicating the probability for each cell to belong to one lineage.
#'   All \code{warpDE} methods can take an object of the class
#'   \code{warpDEDataSet} as input and will output the same.
#'
#' @import methods
#' @export
#'
setClass(
  Class = "warpDEDataSet",
  representation =  representation(
    counts = "matrix",
    t = "matrix",
    w = "matrix"
  ),
  prototype(counts = NULL, t = NULL, w = NULL)
)

setValidity("warpDEDataSet", function(object){
  X <- counts(object)
  n <- nrow(X)
  m <- ncol(X)

  if(!is.numeric(X)) {
    return('counts coordinates must be numeric.')
  }
  if(n==0){
    return('counts has zero rows.')
  }
  if(m==0){
    return('counts has zero columns.')
  }
  if (nrow(t(object))!= m){
    return('pseudotimes should be of the same length as the number of cells.')
  }
  if (nrow(w(object))!= m){
    return('weights should be of the same length as the number of cells.')
  }
  if(ncol(t(object))!= 2){
    return('pseudotimes should be for exactly two lineages and thus have two columns.')
  }
  if(ncol(w(object))!= 2){
    return('weights should be for exactly two lineages and thus have two columns.')
  }
  # something requires row and column names. Automatically input if there is none
  if(is.null(rownames(X))){
    rownames(counts(object)) <- paste('gene',
                                          seq_len(nrow(X)),
                                          sep='-')
  }
  if(is.null(colnames(X))){
    colnames(counts(object)) <- paste('Cell',
                                          seq_len(ncol(X)),
                                          sep='-')
  }
})


#' @title Class \code{rankingDE}
#' @aliases rankingDE
#'
#' @description The \code{rankingDE} class holds da ranking format for the \code{warpDE} package
#'  It is crafted to be part of the workflow for differentially expressed gen inference and the exploration of the genes ranking.

#' @import methods
#' @export
#'
setClass(
  Class = "rankingDE",
  representation =  representation(
    ranking.df = "data.frame",
    params = "list"
  ),
  prototype(ranking.df = NULL, params = NULL)
)
