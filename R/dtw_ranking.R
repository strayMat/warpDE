#' @title Compute a ranking of a list of genes based on the dtw distance between two lineages
#' @name dtw_rank
#'
#' @description
#'
#' @param data a \code{lineageDEDataSet} with genes to be ranked.
#' @param gene character, a gene of interest.
#' @param reg.f a function to perform regression, either "ns" for natural splines, "loess" or "splines" (default is ns).
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param norm  character,("L2" or "L1") the norm to be used for the dtw distance (default is "L2")
#' @param window.size integer, the size of the warping window (default is NULL, no window.size), see details.
#' @param nb.prediction.points integer, if we want to specify the same number of points for each lineage to compute the dtw distance with (default is NULL, we use all points available in each lineage).
#' @param Zscore logical, decide to compute or not the dtw distance with a y-axis normalization of the data (default is FALSE).
#'
#' @return returns a \code{rankingDE} object.
#'
#' @importFrom dtwclust dtw_basic
#' @export
dtw_rank <- function(data,
                     reg.f = "ns",
                     span = 0.75,
                     s.df = 4,
                     norm ="L2",
                     window.size = NULL,
                     equal.size = NULL,
                     Zscore = F){
  logCounts <- log1p(data@counts)
  w <- data@w
  t <- data@t
  n <- nrow(logCounts)
  dtw_dist <- array(NA, dim = n)
  rownames(dtw_dist) <- rownames(logCounts)
  #discard unappropriate cells for new predictions (keep only w >0.5 approximately)
  middle_w1 <- w[w[,1]!=0 & w[,1]!=1,1]
  middle_w2 <- w[w[,2]!=0 & w[,2]!=1,2]
  cells_pred1 <- names(w[w[,1] > (0.5 + sd(middle_w1)),1])
  cells_pred2 <- names(w[w[,2] > (0.5 + sd(middle_w2)),2])
  # eventuellement, subset the cells to have the same population
  if (!is.null(equal.size)){
    cells_pred1 <- sample(cells_pred1, equal.size)
    cells_pred2 <- sample(cells_pred2, equal.size)
  }
  for (g in 1:n){
    gene <- rownames(logCounts)[g]
    y <- logCounts[g, ]

    alt <- reg_gam(data, gene, reg.f = reg.f, span = span, s.df = s.df, null.model = F)$reg$alt
    y1new <- predict(alt, data.frame(x.fit = t[cells_pred1,1], lineage = rep(1, length(cells_pred1))))
    y1new <- y1new[order(t[cells_pred1,1])]
    y2new <- predict(alt, data.frame(x.fit = t[cells_pred2,2], lineage = rep(2, length(cells_pred2))))
    y2new <- y2new[order(t[cells_pred2,2])]

    if (Zscore == T){
      y1new <- zscore(y1new)
      y2new <- zscore(y2new)
    }
    dtw_dist[g] <- dtw_basic(y1new,y2new, window.size = window.size, norm = norm)
  }

  dtw_ranking <- data.frame(dist = dtw_dist )
  rownames(dtw_ranking) <- rownames(logCounts)
  dtw_ranking$rank <- n - rank(dtw_ranking$dist) + 1
  dtw_ranking <- dtw_ranking[order(dtw_ranking$rank),]
  if (reg.f =="loess"){smooth.param = span}
  if (reg.f =="ns"){smooth.param = s.df}
  result <- new("rankingDE", ranking.df = dtw_ranking, params = list(method = paste(reg.f, "dtw"), smooth.param = smooth.param, window.size = window.size, dtw.norm = norm))
  return(result)
}



#' @title Align with dtw two lineages
#' @name dtw_align
#'
#' @description Align two lineage with dtw, fit the aligned data and plot the results
#'
#' @param data a \code{lineageDEDataSet} with genes to be ranked.
#' @param gene character, a gene of interest.
#' @param reg.f a function to perform regression, either "ns" for natural splines, "loess" or "splines" (default is ns).
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param norm  character,("L2" or "L1") the norm to be used for the dtw distance (default is "L2").
#' @param window.type charcater, the window type for the computation of dtw, (see \code{dtw::dtw} for precisions; default is "none").
#' @param window.size integer, the size of the warping window (default is NULL, no window.size), see details.
#' @param align.show logical, if we want to plot the alignment.
#' @param legend.show logical, if the legend is wanted (default is FALSE).
#'
#' @return returns a \itemize{\item \code{pl} a plotof the aligned lineages. \item \code{reg} the alternative and null model for likelihood comaprison.}
#'
#' @importFrom dtw dtw
#' @export
dtw_align <- function(data,
                      gene,
                      reg.f = "ns",
                      span = 0.75,
                      s.df = 4,
                      norm ="L2",
                      window.type = "none",
                      window.size = 50,
                      align.show = F,
                      legend.show = F){
  logCounts <- log1p(data@counts)
  w <- data@w
  t <- data@t
  n <- nrow(logCounts)
#discard unappropriate cells for new predictions (keep only w >0.5 approximately)
  middle_w1 <- w[w[,1]!=0 & w[,1]!=1,1]
  middle_w2 <- w[w[,2]!=0 & w[,2]!=1,2]
  middle_w = c(middle_w1, middle_w2)
  cells_pred1 <- names(w[w[,1] > (0.5 + sd(middle_w)),1])
  cells_pred2 <- names(w[w[,2] > (0.5 + sd(middle_w)),2])
  cells_shared1 <- names(w[w[,1] <= (0.5 + sd(middle_w)) & w[,1] != 0,1])
  cells_shared2 <- names(w[w[,2] <= (0.5 + sd(middle_w)) & w[,2] != 0,2])
  ## gene level
  y <- logCounts[gene, ]
  alt <- reg_gam(data, gene, reg.f = reg.f, span = span, s.df = s.df, null.model = F)$reg$alt
  y1new <- predict(alt, data.frame(x.fit = t[cells_pred1,1], lineage = rep(1, length(cells_pred1))))
  y1new <- y1new[order(t[cells_pred1,1])]
  y2new <- predict(alt, data.frame(x.fit = t[cells_pred2,2], lineage = rep(2, length(cells_pred2))))
  y2new <- y2new[order(t[cells_pred2,2])]
  window.size <- floor((length(y1new) +  length(y2new))/10)
  dtw_dist <- dtw(y1new,y2new, window.type = window.type, window.size = window.size)
  align1 <- dtw_dist$index1
  align2 <- dtw_dist$index2
  rm.1 <- duplicated(align1)
  rm.2 <- duplicated(align2)
  y_w1 <- y[cells_pred1][align1]
  y_w2 <- y[cells_pred2][align2]
  w_w1 <- w[cells_pred1,1][align1]
  w_w2 <- w[cells_pred2,2][align2]
  # renornalize the unshared part to be the same proporation compared to the shared part as in primary data
  t_warp <- (1:length(align1) + max(t[cells_shared1,1],t[cells_shared2,2]))*(max(t[cells_pred1,1],t[cells_pred2,2])-max(t[cells_shared1,1],t[cells_shared2,2]))/length(align1)
  t_w1 <- t_warp
  t_w2 <- t_warp
  # compute a mean time for the replicated cells
  for (e in y_w1){
    t_w1[y_w1==e] <- mean(t_w1[y_w1==e])
  }
  for (e in y_w2){
    t_w2[y_w2==e] <- mean(t_w2[y_w2==e])
  }
  # remove the duplicated cells
  y_w1 <- y_w1[!rm.1]
  y_w2 <- y_w2[!rm.2]
  t_w1 <- t_w1[!rm.1]
  t_w2 <- t_w2[!rm.2]
  w_w1 <- w_w1[!rm.1]
  w_w2 <- w_w2[!rm.2]

  # add the shared part
  y_warp1 <- c(y[cells_shared1], y_w1)
  y_warp2 <- c(y[cells_shared2], y_w2)
  w_warp1 <- c(w[cells_shared1,1], w_w1)
  w_warp2 <- c(w[cells_shared2,2], w_w2)
  t_warp1 <- c(t[cells_shared1,1], t_w1)
  t_warp2 <- c(t[cells_shared2,2], t_w2)

  reg.df <- data.frame(y.fit = c(y_warp1, y_warp2), x.fit = c(t_warp1, t_warp2), w.fit = c(w_warp1, w_warp2), lineage = c(rep(1, length(y_warp1)), rep(2, length(y_warp2))))

  if (reg.f == "loess"){
    alt <- gam(y.fit ~ lineage + lo(x.fit, span = span) + lo(x.fit, span = span):lineage , weights = w.fit, data = reg.df)
  }
  else if (reg.f == "ns"){
    alt <- gam(y.fit ~ lineage + ns(x.fit, s.df) + ns(x.fit, s.df):lineage , weights = w.fit, data = reg.df)
  }
  else if (reg.f == "splines"){
    alt <- gam(y.fit ~ lineage + s(x.fit, 4) + s(x.fit, 4):lineage , weights = w.fit, data = reg.df)
  }
  regs <- list(alt = alt)
  if (reg.f == "loess"){
    null.m <- gam(y.fit ~ lo(x.fit, span = span), weights = w.fit, data = reg.df)
  }
  else if (reg.f == "ns"){
    null.m <- gam(y.fit ~ ns(x.fit, df = s.df), weights = w.fit, data = reg.df)
  }
  else if (reg.f == "splines"){
    null.m <- gam(y.fit ~ s(x.fit, df = 4), weights = w.fit, data = reg.df)
  }
  regs$null <- null.m

  if (align.show == T){

    pl <- ggplot() +
      geom_point(aes(t_warp2, y_warp2), col = "#377EB8", alpha = w_warp2, shape = 1) +
      geom_point(aes(t_warp1, y_warp1), col = "#E41A1C", alpha = w_warp1, shape = 1) +
      ggtitle(gene, paste(reg.f, "regression warped")) + coord_cartesian(ylim=c(-1.5,10))+ xlab("times") + ylab("logcounts")

    t1new <- seq(0, max(t_warp1), length.out = length(t_warp1))
    t2new <- seq(0,max(t_warp2), length.out = length(t_warp2))
    y_pred.alt1 <- predict(alt, data.frame(x.fit = t1new, lineage = rep(1, length(t1new))))
    y_pred.alt2 <- predict(alt, data.frame(x.fit = t2new, lineage = rep(2, length(t2new))))
    plalt <- pl + geom_line(aes(t1new, y_pred.alt1), col = "#E41A1C") + geom_line(aes(t2new, y_pred.alt2), col = "#377EB8")
    y_pred.null <- predict(null.m, data.frame(x.fit = t1new))
    plalt <- plalt + geom_line(aes(t1new, y_pred.null, colour = "null model"), linetype = 2, na.rm = T)
    if (legend.show == F){
      plalt <- plalt + theme(legend.position = "none")
    }
    return(list(pl = plalt, reg = regs))
  }
  return(regs)
}




















# Bootstraps over cells : for confidence interval and sd
cells_bootstrap <- function(df, ranking.subset, nb_bootstrap = 1000){
  sub.size <- nrow(ranking.subset)
  prediction_length <- 100
  bs_cells_dtw <- matrix(nrow = nb_bootstrap, ncol = sub.size)
  genes <- rownames(ranking.subset)
  colnames(bs_cells_dtw) <- genes

  for (i in 1:nb_bootstrap){
    # resample cells
    bs.sample <- sample(1:ncol(df$log_counts), replace =  T)
    bs.df <- list(log_counts = df$log_counts[,bs.sample], w = df$w[bs.sample,], t = df$t[bs.sample,])
    t1new <- seq(0, max(bs.df$t[,1]), length.out = prediction_length)
    t2new <- seq(0, max(bs.df$t[,2]), length.out = prediction_length)
    bs.dtw <- array(sub.size)
    for (g in 1:length(genes)){
      bs.lo1 <- loess(bs.df$log_counts[genes[g],] ~ bs.df$t[,1], weights = bs.df$w[,1])
      bs.lo2 <- loess(bs.df$log_counts[genes[g],] ~ bs.df$t[,2], weights = bs.df$w[,2])
      bs.dtw[g] <- dtw_basic(predict(bs.lo1, t1new), predict(bs.lo2, t2new), window.size = 20, norm = "L2")
      # # #plot curve
      # loess_regression(genes[g], df = df)$pl +
      #  geom_line(aes(t1new, predict(bs.lo1, t1new))) +
      #  geom_line(aes(t2new, predict(bs.lo2, t2new)))
    }
    bs_cells_dtw[i,] <- bs.dtw
  }
  return(list(bs.distribs = bs_cells_dtw, original_ranking = ranking.subset))
}
