#' @title Compute pvalues and aic difference criteria
#' @name likelihood_criteria
#'
#' @description Perform a regression with \code{reg_vgam} or \code{reg_loess} function and compute p-values based on a likelihood ratio test as well as a AIC based criteria (computation of the difference of AIC between the null model and the alternative model).
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param reg the regression method to use, either "loess" with the gam package or "splines" which fits a spline with gam package.
#' @param pval logical, if the likelihood criteria is to be computed or not (default is TRUE).
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param fam character, for the vgam regression, the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").

#' @return returns the pvalue and if asked the aic difference for the gene of intereset.
#'
#' @importFrom VGAM logLik
#' @importFrom VGAM AIC
#' @importFrom VGAM df.residual
#' @export

likelihood_criteria <- function(data, gene, reg, pval = T, span = 0.75, s.df = 4, fam = "gaussian"){
  regs <- reg_gam(data, gene, reg, span = span, s.df = s.df)$reg
  reg.null <- regs$null
  reg.alt <- regs$alt
  res <- list()
  res$aic.diff <- AIC(reg.null) - AIC(reg.alt)
  if (pval == T){
    testStatistic <- 2* (logLik(reg.alt)[1] - logLik(reg.null)[1])
    ####### Verify what genes it gives
    res$pval <- pchisq(testStatistic, df = df.residual(reg.null) - df.residual(reg.alt) , lower.tail = F)
  }
  return(res)
}

#' @title Vgam ranks
#' @name likelihood_rank
#'
#' @description compute vgam ranks with likelihood ratio test p-values and aic difference criteria
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param reg the regression method to use, either "loess" with the gam package or "vgam" which fits a spline with loess package.
#' @param pval logical, if the pval criteria is to be computedor not (default is TRUE).
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param fam character, for the vgam regression, the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").
#'
#' @return returns \code{rankingDE} objects: one for the likelihood ratio test pvalues and one for the aic difference criteria if aic == T.
#'
#' @export

likelihood_rank <- function(data,
                            reg.f = "splines",
                            pval = T,
                            span = 0.75,
                            s.df = 4,
                            fam = "gaussian"){
  res <- list()
  criteria <- sapply(rownames(data@counts), function(x) likelihood_criteria(data, x, reg.f, span = span, s.df = s.df, fam = fam))
  aic <- unlist(criteria[1,])
  ranking_aic <- data.frame(aic.diff = aic, rank = length(aic) - rank(aic) + 1)
  ranking_aic <- ranking_aic[order(ranking_aic$rank),]
  if (reg.f =="loess"){smooth.param = span}
  if (reg.f =="splines"){smooth.param = s.df}
  params <- list(method = paste(reg.f, "AIC diff"), smooth.param = smooth.param, fam = fam)
  res$aic <- new("rankingDE", ranking.df = ranking_aic, params = params)
  if (pval ==T){
    pvalues <- unlist(criteria[2,])
    ranking_pval <- data.frame(pval = pvalues, rank = rank(pvalues))
    ranking_pval <- ranking_pval[order(ranking_pval$rank),]
    params$method <- paste(reg.f, "pval")
    res$pval <- new("rankingDE", ranking.df = ranking_pval, params = params)
  }
  return(res)
}
