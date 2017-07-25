#' @title Compute pvalues and aic difference criteria
#' @name likelihood_criteria
#'
#' @description Perform a regression with \code{reg_vgam} or \code{reg_loess} function and compute p-values based on a likelihood ratio test as well as a AIC based criteria (computation of the difference of AIC between the null model and the alternative model).
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param reg the regression method to use, either "loess" with the gam package or "vgam" which fits a spline with loess package.
#' @param aic logical, if the aic criteria is to be computedor not (default is TRUE).
#' @param span numeric, for the loess regression, control the ammount of regularization for the loess regression.
#' @param fam character, for the vgam regression, the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").

#' @return returns the pvalue and if asked the aic difference for the gene of intereset.
#'
#' @importFrom VGAM logLik
#' @importFrom VGAM AIC
#' @importFrom VGAM df.residual
#' @export

likelihood_criteria <- function(data, gene, reg, aic = T, span = 0.5, fam = "gaussian"){
  if (reg == "loess"){
    regs <- reg_loess(data, gene, span = span)$reg
    method = list(reg = reg, span = span)
  }
  else if (reg == "vgam"){
    regs <- reg_vgam(data, gene, fam)$reg
    method <- list(reg = reg, fam = fam)
  }
  reg.null <- regs$null
  reg.alt <- regs$alt
  testStatistic <- 2* (logLik(reg.alt)[1] - logLik(reg.null)[1])
  ############ BIG question here concerning the null model and its degree of freedom################
  pval <- pchisq(testStatistic, df = df.residual(reg.null) - df.residual(reg.alt) , lower.tail = F)
  res <- list(pval = pval)
  if (aic == T){
    res$aic.diff <- AIC(reg.null) - AIC(reg.alt)
  }
  res$method <- method
  return(res)
}

#' @title Vgam ranks
#' @name likelihood_rank
#'
#' @description compute vgam ranks with likelihood ratio test p-values and aic difference criteria
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param reg the regression method to use, either "loess" with the gam package or "vgam" which fits a spline with loess package.
#' @param aic logical, if the aic criteria is to be computedor not (default is TRUE).
#' @param span numeric, for the loess regression, control the ammount of regularization for the loess regression.
#' @param fam character, for the vgam regression, the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").
#'
#' @return returns \code{rankingDE} objects: one for the likelihood ratio test pvalues and one for the aic difference criteria if aic == T.
#'
#' @export

likelihood_rank <- function(data,
                      reg = "loess",
                      aic = T,
                      span = 0.5,
                      fam = "gaussian"
                      ){
  res <- list()
  criteria <- sapply(rownames(data@counts), function(x) likelihood_criteria(data, x, reg))
  pvalues <- unlist(criteria[1,])
  method <- as.list(unlist(criteria[3,1]))
  ranking_pval <- data.frame(pval = pvalues, rank = length(pvalues) - rank(pvalues) + 1)
  ranking_pval <- ranking_pval[order(ranking_pval$rank),]
  method$method <- "pval"
  res$pval <- new("rankingDE", ranking.df = ranking_pval, params = method)
  if (aic ==T){
    aic <-unlist(criteria[2,])
    ranking_aic <- data.frame(aic.diff = aic, rank = length(aic) - rank(aic) + 1)
    ranking_aic <- ranking_aic[order(ranking_aic$rank),]
    method$method <- "AIC difference"
    res$aic <- new("rankingDE", ranking.df = ranking_aic, params = method)
  }
  return(res)
}
