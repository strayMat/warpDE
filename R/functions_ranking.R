#library(VGAM)
# Compare ranking : two functions to give an idea of the similarity between two ranking.
# Another good measure of similarity is the kendall's tau which counts the number of permutations.
# these function takes as input an output of the format of the dtw_ranking function:
# ranking = dataframe(dist = distance values for the ranking, rank = ranks of the elements )
# rownames(ranking) are gene names

compare_ranking <- function(ranking1, ranking2, nmax = min(dim(ranking1)[1], dim(ranking2)[1])){
  best1 <- ranking1[order(ranking1[,2]),]
  best2 <- ranking2[order(ranking2[,2]),]
  ind <- round(seq(10,nmax, length.out = 20),0)
  inter <- matrix(NA,nrow = 1, ncol = length(ind))
  for (i in 1:length(ind)){
    inter[1,i] <- round(length(intersect(rownames(best1)[1:ind[i]], rownames(best2)[1:ind[i]]))/ind[i],2)
  }
  colnames(inter) <- ind
  return(inter)
}

###### More visual
plot_2rankings <- function(ranking1, ranking2, nmax = min(dim(ranking1)[1], dim(ranking2)[1])){
  xlab <- deparse(substitute(ranking1))
  ylab <- deparse(substitute(ranking2))
  ranking1 <- ranking1[ranking1[,2]<nmax,]
  ranking2 <- ranking2[ranking2[,2]<nmax,]
  commons <- intersect(rownames(ranking1), rownames(ranking2))
  sc_plot <-ggplot() + geom_point(aes(ranking1[commons, 2], ranking2[commons, 2]), col = rgb(0,0,0,0.4))+
    xlab(xlab)  + ylab(ylab) + geom_line(aes(x = 1:max(ranking1$rank, ranking2$rank), y = 1:max(ranking1$rank, ranking2$rank)), col = "red")
  return(sc_plot)
}

###### kendall tau adapted to my format of rankings
ken_tau <- function(ranking1,ranking2, nmax = min(dim(ranking1)[1], dim(ranking2)[1])){
  ranking1 <- ranking1[ranking1[,2]<nmax,]
  ranking2 <- ranking2[ranking2[,2]<nmax,]
  commons <- intersect(rownames(ranking1), rownames(ranking2))
  return(kendall.tau(ranking1[commons,2], ranking2[commons,2]))
}

####### Elbow curve for quick vizualization of the rankings
elbow_curve <- function(ranking, xmax = nrow(ranking)){
  x <- seq(1,xmax)
  y <- as.numeric(ranking$dist[order(-ranking$dist)])[1:xmax]

  logscale <- grid.arrange(ggplot() + geom_point(aes(x = x, y = y)) + xlab("rank of gene") + ylab("dtw distance") + expand_limits(x = 0, y =0) + ggtitle("dtw distances against ranking"), ggplot() + geom_point(aes(x = x, y =y)) + xlab("rank of gene") + ylab("dtw distance") + expand_limits(x = 0, y =0) + scale_y_log10()+ ggtitle("log10(dtw) distances against ranking"))

  # "Best" subset of the data to fit an exponential
  xfit <- 50
  y <- as.numeric(ranking$dist[order(-ranking$dist)])[1:xmax]
  d <- ggplot() + geom_point(aes(x = x, y = y)) + xlab("rank of gene") + ylab("dtw distance") + expand_limits(x = 0, y =0)
  tan <- lm(y[1:xfit] ~ x[1:xfit])
  s <-summary(tan)
  #threshold
  thres <- round(as.numeric(-tan$coefficients[1]/tan$coefficients[2]))
  thres_plot <- d + geom_abline(intercept = tan$coefficients[1], slope = tan$coefficients[2], color = "red")+
    annotate("text", label = paste("threshold:", thres), x  = thres + xmax/13, y = - 0.005, size = 4)+
    geom_point(aes(x = thres, y = 0), color = "red")
  return(list(pl = thres_plot, thres = thres, logpl = logscale))
}
