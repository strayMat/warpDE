#' @title Compute vgam pvalues and aic difference criteria
#' @name criteria_vgam
#'
#' @description Perform vgam regression with \code{reg_vgam} function and compute p-values based on a likelihood ratio test as well as a AIC based criteria (computation of the difference of AIC between the null model and the alternative model).
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param fam the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").
#' @param aic logical, if the aic criteria is to be computedor not (default is TRUE).
#'
#' @return returns the pvalue and if asked the aic difference for the gene of intereset.
#'
#' @import VGAM
#' @export

pval_vgam <- function(data, gene, fam = "gaussian", aic = T){
  v_reg <- reg_vgam(data, gene, fam)
  v_reg1 <- v_reg$reg$spl1
  v_reg2 <- v_reg$reg$spl2
  v_reg.d <- v_reg$reg$spl.d
  testStatistic <- -2* (logLik(v_reg.d) - logLik(v_reg2) - logLik(v_reg1))
  pval <- pchisq(testStatistic, df = (length(coef(v_reg2)) + length(coef(v_reg1)))- length(coef(v_reg.d)))
  if (aic == T){
    aic1 <- AIC(v_reg1) + AIC(v_reg2)
    aic.diff <- AIC(v_reg.d) - aic1
    return(list(pval = pval, aic.diff = aic.diff))
  }
  return(pval)
}

#' @title Vgam ranks
#' @name vgam_rank
#'
#' @description compute vgam ranks with likelihood ratio test p-values and aic difference criteria
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param fam the distribution assumption of the residuals; etiher "binomial" or "gaussian" (default is "gaussian").
#' @param aic logical, if the aic criteria is to be computedor not (default is TRUE).
#'
#' @return returns \code{rankingDE} objects: one for the likelihood ratio test pvalues and one for the aic difference criteria if aic == T.
#'
#' @import VGAM
#' @export

vgam_rank <- function(data,
                      fam = "gaussian",
                      aic = T){
  res <- list()
  pval_v <- sapply(rownames(data@counts), function(x) pval_vgam(data, x))
  # only one for the p-values which is disturbing for now
  pvalues_v <- unlist(pval_v[1,])
  v_gaus_ranking_pval <- data.frame(pval = pvalues_v, rank = length(pvalues_v) - rank(pvalues_v) + 1)
  v_gaus_ranking_pval <- v_gaus_ranking_pval[order(v_gaus_ranking_pval$rank),]
  res$pval_vgam <- new("rankingDE", ranking.df = v_gaus_ranking_pval, params = list(method = "vgam likelihood", fam = fam, criteria = "pval"))
  if (aic ==T){
    aic_v <-unlist(pval_v[2,])
    v_gaus_ranking_aic <- data.frame(aic.diff = aic_v, rank = length(aic_v) - rank(aic_v) + 1)
    v_gaus_ranking_aic <- v_gaus_ranking_aic[order(v_gaus_ranking_aic$rank),]
    res$aic_vgam <- new("rankingDE", ranking.df = v_gaus_ranking_aic, params = list(method = "vgam likelihood", fam = fam, criteria = "AIC"))
  }
  return(res)
}



####### Compute alternative and null models and compute their residuals and
#number of parameters and standard error

build_residuals <- function(data, M0 = F){
  n = dim(data$log_counts)[1]
  m = dim(data$log_counts)[2]
  l1_cells <- names(data$t[!is.na(data$t[,1]),1])
  l2_cells <- names(data$t[!is.na(data$t[,2]),2])

  resid1 <- data.frame(matrix(NA, nrow = n, ncol = length(l1_cells)))
  colnames(resid1) <- l1_cells
  rownames(resid1) <- rownames(data$log_counts)
  resid2 <- data.frame(matrix(NA, nrow = n, ncol = length(l2_cells)))
  colnames(resid2) <- l2_cells
  rownames(resid2) <- rownames(data$log_counts)
  residd <- data.frame(matrix(NA, nrow = n, ncol = length(l1_cells) +
                                length(l2_cells)))
  colnames(residd) <- c(l1_cells, paste0(l2_cells,"l2"))
  rownames(residd) <- rownames(data$log_counts)
  nPam <- matrix(NA_integer_, nrow = n, ncol = 3)
  rss <- matrix(NA_integer_, nrow = n, ncol = 3)

  if (M0 == T){
    resid0 <- data.frame(matrix(NA, nrow = n, ncol = m))
    colnames(resid0) <- names(times[,3])
    rownames(resid0) <- rownames(data$log_counts)
    nPam <- matrix(NA_integer_, nrow = n, ncol = 4)
    rss <- matrix(NA_integer_, nrow = n, ncol = 4)
    }

  ## Compute the regression and store the variance matrix as well as the degree opf freedom
  for (g in 1:n){
    gene <- rownames(data$log_counts)[g]
    viz <- loess_regression(gene, df = data, MD = T)
    resid1[g,l1_cells] <- t(array(viz$reg$lo1$y-predict(viz$reg$lo1, viz$reg$lo1$x)))
    resid2[g, l2_cells] <- t(array(viz$reg$lo2$y-predict(viz$reg$lo2, viz$reg$lo2$x)))
    if (M0  == T){
      resid0[g, ] <- t(array(viz$reg$lo0$y-predict(viz$reg$lo0, viz$reg$lo0$x)))
    }
    residd[g,] <- t(array(viz$reg$lod$y-predict(viz$reg$lod, viz$reg$lod$x)))
    nPam[g,] <- sapply(viz$reg, function(x) x$enp)
    rss[g,] <- sapply(viz$reg, function(x) x$s)
  }
  if (M0 == T){
    return(list(resid1 = resid1, resid2 = resid2, resid0 = resid0, residd = residd,
                nb_param = nPam, rss = rss))
  }
  else{
    return(list(resid1 = resid1, resid2 = resid2, residd = residd,
                nb_param = nPam, rss = rss))
  }
}


##### Compute statistics based on the Likelihood estimations and a gaussian model
comp_stats <- function(data, residuals, AIC_pval = T){
  res <- list()
  n <-dim(data$log_counts)[1]
  l1_cells <- names(data$t[!is.na(data$t[,1]),1])
  l2_cells <- names(data$t[!is.na(data$t[,2]),2])
  pvald <- matrix(ncol = 1, nrow = n)
  colnames(pvald) <- "pval_d"
  aicd <- matrix(ncol = 1, nrow = n)
  colnames(aicd) <- "aic_d"
  w1 <- data$w[l1_cells,1]
  w2 <- data$w[l2_cells,2]
  w1sum <-sum(w1)
  w2sum <- sum(w2)
  wconcat <- c(w1,w2)
  n_l1 <- length(l1_cells)
  n_l2 <- length(l2_cells)
  for (g in 1:n){
    gene <- rownames(counts)[g]
    # weighted residuals
    e1_bar <- sum(w1*residuals$resid1[g,]**2)
    e2_bar <- sum(w2*residuals$resid2[g,]**2)
    # - 2loglikelihood for each model, likelihoord ratio test hand made
    L1 <- w1sum*log(2*pi) + w1sum*log(residuals$rss[g,1]**2) + e1_bar/(residuals$rss[g,1]**2)
    L2 <- w2sum*log(2*pi) + w2sum*log(residuals$rss[g,2]**2) + e2_bar/(residuals$rss[g,2]**2)
    # for the whole alternative model
    L_alt  <- L1 + L2
    Ld <- (n_l1 + n_l2)*log(2*pi) + (n_l1 + n_l2)*log(residuals$rss[g,3]**2)
    + sum(wconcat*residuals$residd[g,]**2)/(residuals$rss[g,3]**2)
    # Chi-sauared Pvalues of the likelihood ratio tests :
    pvald[g,] <- 1-pchisq(Ld - L_alt, (residuals$nb_param[g,1] + residuals$nb_param[g,2]) - residuals$nb_param[g,3])
    # AIC handmade:
    aic_alt <- L_alt + 2*(residuals$nb_param[g,1] + residuals$nb_param[g,2])
    aic_d <- Ld + 2*residuals$nb_param[g,3]
    aicd[g,] = aic_d-aic_alt
  }
  if (!is.null(residuals$resid0)){
    pval0 <- matrix(ncol = 1, nrow = n)
    colnames(pvalues) <- "pval_0"
    aic0 <- matrix(ncol = 1, nrow = n)
    colnames(aic) <-"aic_0"
    for (g in 1:n){
      gene <- rownames(counts)[g]
      # weighted residuals
      e1_bar <- sum(w1*residuals$resid1[g,]**2)
      e2_bar <- sum(w2*residuals$resid2[g,]**2)
      # - 2loglikelihood for each model, likelihoord ratio test hand made
      L1 <- w1sum*log(2*pi) + w1sum*log(residuals$rss[g,1]**2) + e1_bar/(residuals$rss[g,1]**2)
      L2 <- w2sum*log(2*pi) + w2sum*log(residuals$rss[g,2]**2) + e2_bar/(residuals$rss[g,2]**2)
      # for the whole alternative model
      L_alt  <- L1 + L2
      # nul models : not the same number of points than the double model and the alternative one...
      L0 <- n*log(2*pi) + n*log(residuals$rss[g,4]**2) + sum(residuals$resid0[g,]**2)/(residuals$rss[g,4]**2)
      # Chi-sauared Pvalues of the likelihood ratio tests :
      pval0[g,] <- 1-pchisq(L0 - L_alt, (residuals$nb_param[g,1] + residuals$nb_param[g,2]) - residuals$nb_param[g,4])
      # AIC handmade:
      aic_alt <- L_alt + 2*(residuals$nb_param[g,1] + residuals$nb_param[g,2])
      aic_0 <- L0 + 2*residuals$nb_param[g,4]
      aic0[g,] = aic_0-aic_alt
    }
    # construct ranking objects
    aic0_ranking <- data.frame(dist = aic0, rank =n+1-rank(aic0))
    rownames(aic0_ranking) <- rownames(counts)
    aic0_ranking <- aic0_ranking[order(aic0_ranking$rank),]

    p0_ranking <- data.frame(dist = pval0, rank = rank(pval0))
    rownames(p0_ranking) <- rownames(counts)
    p0_ranking <- p0_ranking[order(p0_ranking$rank),]
    p0_ranking$rank <- 1:n
    res$pval_0 = p0_ranking
    res$aic_0 = aic0_ranking
  }
  aicd_ranking <- data.frame(dist = aicd, rank = n+1-rank(aicd))
  rownames(aicd_ranking) <- rownames(counts)
  aicd_ranking <- aicd_ranking[order(aicd_ranking$rank),]

  pd_ranking <- data.frame(dist = pvald, rankm = rank(pvald))
  rownames(pd_ranking) <- rownames(counts)
  pd_ranking <- pd_ranking[order(pd_ranking$rank),]
  pd_ranking$rank <- 1:n

  res$pval_d = pd_ranking
  res$aic_d = aicd_ranking

  return(res)
}


F_stat_bs <- function(res, data, g.subset, niter = 1000){
  n <- length(g.subset)
  m <- ncol(data$log_counts)
  l1_cells <- names(data$t[!is.na(data$t[,1]),1])
  l2_cells <- names(data$t[!is.na(data$t[,2]),2])
  w1 <- data$w[l1_cells,1]
  w2 <- data$w[l2_cells,2]
  wconcat <-c(w1, w2)
  tconcat <- c(data$t[l1_cells,1], data$t[l2_cells,2])

  F.stats <- matrix(ncol = n, nrow = niter)
  colnames(F.stats) <- g.subset

  # compute original Fstats
  rss1.o <- sapply(1:n, function(x) sum(w1*(res$resid1[g.subset[x],])**2) + sum(w2*(res$resid2[g.subset[x],])**2))
  rss0.o <- sapply(1:n, function(x) sum(wconcat*(res$residd[g.subset[x],])**2))

  F.stat.o <- (rss0.o-rss1.o)/rss1.o

  # bootstrap
  for (i in 1:niter){
    # sample_indices
    bs.ind1 <- sample(1:length(l1_cells), replace = T)
    bs.ind2 <- sample(1:length(l2_cells), replace = T)
    for (g in 1:n){
      gene <- g.subset[g]
      e1 <- res$resid1[gene, bs.ind1]
      e2 <- res$resid2[gene, bs.ind2]
      ed <- res$residd[gene,]
      ypredd <- c(data$log_counts[gene,l1_cells], data$log_counts[gene,l2_cells]) - ed
      y.bs <- unlist(ypredd + c(e1,e2))
      lo.bs <- loess(y.bs ~ tconcat, weights = wconcat, span = 0.75)
      rss.bs <- sum(wconcat*lo.bs$residuals**2)
      F.stats[i, g] <- (rss.bs - rss1.o[g]) / rss1.o[g]
    }
  }

  p.val <- 1 - sapply(1:n, function(x) ecdf(F.stats[,x])(F.stat.o[x]))
  names(p.val) <- g.subset
  return(list(F.distrib = F.stats, F.stat.o = F.stat.o, p.val = p.val))
}

# # example, need data
# gene <- rownames(counts)[1]
# t <- loess_regression(gene, df, MD = T)
# e<- c(df$log_counts[gene,l1_cells],df$log_counts[gene,l2_cells]) - t$reg$lod$residuals
# ggplot() + geom_point(aes(c(df$t[l1_cells,1],df$t[l2_cells,2]), e), col = "red") +
#   geom_point(aes(c(df$t[l1_cells,1],df$t[l2_cells,2]), t$reg$lod$fitted), col = rgb(0,0,1,alpha = 0.2))

