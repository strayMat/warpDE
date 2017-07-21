#' @title Compute a ranking of a list of genes based on the dtw distance between two lineages
#' @name dtw_rank
#'
#' @description
#'
#' @param data a \code{lineageDEDataSet} with genes to be ranked.
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
    lo1 <- loess(y ~ t[,1], weights = w[,1], span = 0.75)
    lo2 <- loess(y ~ t[,2], weights = w[,2], span = 0.75)
    y1new <- predict(lo1, t[cells_pred1,1])
    y1new <- y1new[order(t[cells_pred1,1])]
    y2new <- predict(lo2, t[cells_pred2,2])
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
  result <- new("rankingDE", ranking.df = dtw_ranking, params = list(method = "dtw_basic", window.size = window.size, dtw.norm = norm))
  return(result)
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


### ne marche pas mais essayer de trouver comment faire pour afficher sur un meme axe les eloginements
dtw_pplot <- function(gene,data, equal.size = NULL, window.size = NULL) {
  middle_w1 <- data$w[data$w[,1]!=0 & data$w[,1]!=1,1]
  middle_w2 <- data$w[data$w[,2]!=0 & data$w[,2]!=1,2]
  cells_pred1 <- names(data$w[df$w[,1] > (0.5 + sd(middle_w1)),1])
  cells_pred2 <- names(data$w[df$w[,2] > (0.5 + sd(middle_w2)),2])
  # eventuellement, subset the cells to have the same population
  if (!is.null(equal.size)){
    cells_pred1 <- sample(cells_pred1, equal.size)
    cells_pred2 <- sample(cells_pred2, equal.size)
  }
  y <- data$log_counts[gene, ]
  lo1 <- loess(y ~ data$t[,1], weights = data$w[,1], span = 0.75)
  lo2 <- loess(y ~ data$t[,2], weights = data$w[,2], span = 0.75)
  y1new <- predict(lo1, data$t[cells_pred1,1])
  y2new <- predict(lo2, data$t[cells_pred2,2])
  dtw_dist <- dtw_basic(y1new, y2new, backtrack = T, window.size = window.size)
  xstr1 <- data$t[cells_pred1,1][dtw_dist$index1]
  xstr2 <- data$t[cells_pred2,2][dtw_dist$index2]
  pl <- ggplot() + geom_line(aes(xstr, y1new[dtw_dist$index1]), col = "red") +
    geom_line(aes(xstr1, y2new[dtw_dist$index2]), col = "blue")
  plot(pl)
  return(dtw)
}
