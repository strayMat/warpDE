
#' @title Wrapper function for all ranking methods
#' @name warpDE_rank
#'
#' @description A wrapper method which takes as input a \code{warpDEDataSet} and various options to return the ranking of the genes following one of three methods available: likelihood ratio test, dynamic time warping distance, or dynamic time warping alignment followed by a likelihood ratio test.
#'
#'
#' @param data a \code{warpDEDataSet} with results to be plotted.
#' @param reg.f character, the regression method to use, either "loess" with the gam package, "ns" which fits a natural cubic spline with gam package or "splines" which fits a smoothing spline (default is "loess").
#' @param span numeric, a smoothing parameter for the regression function if the loess regression (reg.f = 'loess') is used (default is 0.75, see \code{gam::lo} for details about regularization).
#' @param splines.df numeric, a smoothing parameter for the nregression function if natural cubic splines (reg.f = 'ns') are used (default is 4, see \code{splines::s} for details about regularization).
# @param fam character, for the vgam regression, the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").
#' @param ranking_method character, the ranking method to use, either "lrt" for likelihood ratio test, "dtw" for an ordering based solely on dynamic time warping, "warpDE" for a prealable alignment of gene expression with dtw followed by a likelihood ratio test (default is "warpDE").

#' @return returns \code{rankingDE} objects: one for the likelihood ratio test pvalues and one for the aic difference criteria if aic == T.
#'
#' @export

warpDE_rank<- function(data, reg.f = "loess", span = 0.75, splines.df = 4, ranking_method = "warpDE"){
  if (ranking_method == "warpDE"){
    ranking <- likelihood_rank(data, reg.f = reg.f, span = span, s.df = splines.df, dtw = T)$aic
  }
  else if (ranking_method == "dtw"){
    ranking <- dtw_rank(data, reg.f = reg.f, span = span, s.df = splines.df, window.type = "none", window.size = 50)
    }
  else if (ranking_method == "lrt"){
    ranking <- likelihood_rank(data = data, reg.f = reg.f, span = span, s.df = splines.df, dtw = F)$aic
  }
  else {ranking <- "Please enter a valid method: 'warpDE', 'dtw' or 'lrt'."}
  return (ranking)
}
