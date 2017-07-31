#' @title Compare two rankings by counting the shared elements along the rankings.
#' @name quantiles_shared
#'
#' @description Compares two \code{rankingDE} by computing the shared elements in different cumulative quantiles of the rankings.
#'
#' @param ranking1 a first \code{rankingDE} object.
#' @param ranking2 a second \code{rankingDE} object.
#' @param quantiles quantiles of the distribution that we wat to compare (default is is 20-quantiles).
#' @param nmax integer, design the last elemtns where to stop the comparison (default is the minimum length of the two rankings).

#'
#' @return returns a dataframe displaying the number of common elements between the two ranking in the cumulative quantiles.
#'
#' @export

quantiles_shared <- function(ranking1, ranking2, quantiles = c(0.01,seq(0.05,1, length.out = 20)), nmax = min(dim(ranking.df(ranking1))[1], dim(ranking.df(ranking2))[1])){
  method1 <- params(ranking1)$method
  method2 <- params(ranking2)$method
  if (!is.null(params(ranking1)$reg)){
    method1 <- paste(method1, (params(ranking1)$reg))
  }
  if (!is.null(params(ranking2)$reg)){
    method2 <- paste(method2, (params(ranking2)$reg))
  }
  r1 <- ranking.df(ranking1)
  r2 <- ranking.df(ranking2)

  ind <- round(quantiles*nmax,0)
  inter <- matrix(NA,nrow = 1, ncol = length(ind))
  for (i in 1:length(ind)){
    inter[1,i] <- round(length(intersect(rownames(r1)[1:ind[i]], rownames(r2)[1:ind[i]]))/ind[i],2)
  }
  colnames(inter) <- ind
  rownames(inter) <- paste(method1,"vs", method2)
  return(inter)
}

#' @title Compare two rankings by scatter plotting them
#' @name plot_rankings
#'
#' @description Compares two \code{rankingDE} by plotting the ranks of their elements against each other.
#'
#' @param ranking1 a first \code{rankingDE} object.
#' @param ranking2 a second \code{rankingDE} object.
#' @param nmax integer, design the last elemtns where to stop the comparison (default is the minimum length of the two rankings).
#'
#' @return returns a plot comparing two plots rank by rank.
#'
#' @import ggplot2
#' @export
plot_rankings <- function(ranking1, ranking2, nmax = min(dim(ranking.df(ranking1))[1], dim(ranking.df(ranking2))[1])){
  xlab <- params(ranking1)$method
  ylab <- params(ranking2)$method
  if (!is.null(params(ranking1)$reg)){
    xlab <- paste(xlab, (params(ranking1)$reg))
  }
  if (!is.null(params(ranking2)$reg)){
    ylab <- paste(ylab, (params(ranking2)$reg))
  }
  r1 <- ranking.df(ranking1)
  r2 <- ranking.df(ranking2)
  r1 <- r1[r1[,2]<nmax,]
  r2 <- r2[r2[,2]<nmax,]
  commons <- intersect(rownames(r1), rownames(r2))
  sc_plot <-ggplot() + geom_point(aes(r1[commons, 2], r2[commons, 2]), col = rgb(0,0,0,0.4))+
    xlab(xlab)  + ylab(ylab) + geom_line(aes(x = 1:max(r1$rank, r2$rank), y = 1:max(r1$rank, r2$rank)), col = "red")
  return(sc_plot)
}


#' @title Compute kendall's tau between two rankings
#' @name kendall
#'
#' @description Compares two \code{rankingDE} by computing their kendall's tau from package \code{VGAM}.
#'
#' @param ranking1 a first \code{rankingDE} object.
#' @param ranking2 a second \code{rankingDE} object.
#' @param nmax integer, design the last elemtns where to stop the comparison (default is the minimum length of the two rankings).
#'
#'  @details Kendall's tau is basically the number of permutation between the ranks of the elements of two rankings : tau = ((number of concordant pairs) - (number of discordant pairs))/N
#'
#' @return returns Kendall's tau.
#'
#' @importFrom VGAM kendall.tau
#' @export
kendall <- function(ranking1,ranking2, nmax = min(dim(ranking.df(ranking1))[1], dim(ranking.df(ranking2))[1])){
  r1 <- ranking.df(ranking1)
  r2 <- ranking.df(ranking2)
  r1 <- r1[r1[,2]<nmax,]
  r2 <- r2[r2[,2]<nmax,]
  commons <- intersect(rownames(r1), rownames(r2))
  return(kendall.tau(r1[commons,2], r2[commons,2]))
}


#' @title Compare two rankings with different means
#' @name rankings_compare
#'
#' @description Compares two \code{rankingDE} by plotting them against each othe, computing the number of shared elements in the cumulative quantiles of the rankings and computing their kendall's tau.
#'
#' @param ranking1 a first \code{rankingDE} object.
#' @param ranking2 a second \code{rankingDE} object.

#' @return returns the three measures of similarity between the rankings.
#'
#' @export

rankings_compare <- function(ranking1,
                             ranking2,
                             nmax = min(dim(ranking.df(ranking1))[1], dim(ranking.df(ranking2))[1])){
  pl <- plot_rankings(ranking1, ranking2, nmax = nmax)
  shared.q <- quantiles_shared(ranking1, ranking2, nmax = nmax )
  kendalls.tau <- kendall(ranking1, ranking2, nmax = nmax)
  return(list(pl = pl, shared.q = shared.q, kendalls.tau = kendalls.tau))
}


#' @title Elbow curve for a ranking
#' @name elbow_curve
#'
#' @description Plot an elbow curve of the distribution of distances from a \code{rankingDE}. It gives an idea of the shape of the distribution and can help select a threshold for the differntially expressed genes.
#'
#' @param ranking a \code{rankingDE} object.
#' @param xmax integer, at which rank we cut off the curve. (default is length(ranking))
#' @param xfit integer, cutoff for the plot of the tangent at the origin (default is 0.005*length(ranking)).
#'
#' @return returns \itemize{\item{a plot of the distances against the ranks of the elements in the ranking.}\item{rank threshold}}
#'
#' @export
elbow_curve <- function(ranking, xmax = nrow(ranking), xfit = 0.005* nrow(ranking) ){
  ranking <- ranking.df(ranking)
  x <- seq(1,xmax)
  y <- as.numeric(ranking[order(-ranking[,1]),1])[1:xmax]

  xfit <- 50
  y <- as.numeric(ranking[order(-ranking[,1]),1])[1:xmax]
  d <- ggplot() + geom_point(aes(x = x, y = y)) + xlab("rank of gene") + ylab("dtw distance") + expand_limits(x = 0, y =0)
  tan <- lm(y[1:xfit] ~ x[1:xfit])
  s <-summary(tan)
  #threshold
  thres <- round(as.numeric(-tan$coefficients[1]/tan$coefficients[2]))
  thres_plot <- d + geom_abline(intercept = tan$coefficients[1], slope = tan$coefficients[2], color = "red")+
    annotate("text", label = paste("threshold:", thres), x  = thres + xmax/13, y = - 0.005, size = 4)+
    geom_point(aes(x = thres, y = 0), color = "red")
  return(list(pl = thres_plot, thres = thres))
}


#' @title Heatmap visualization of multiple kendall comparisons
#' @name kendall.heatmap
#'
#' @description Compute pairwise kendall taus between the elements of a list from a \code{rankingDE}.
#'
#' @param ranking list, a list of \code{rankingDE} objects.
#' @param labels character, a charcater vector for renaming the rankings to compare (default will use the  parameters of the rankings).
#' @param ... other parameters for the function \code{corrplot}.
#'
#' @return returns a heatmap of kendall's tau
#'
#' @importFrom VGAM kendall.tau
#' @importFrom corrplot corrplot
#' @export
#'
kendall.heatmap <- function(rankings, labels = NULL, addCoef.col = "black", cl.pos = "n", ...){
  lr <- length(rankings)
  tau_matrix <- matrix(nrow = lr, ncol = lr)
  for (r1 in 1:lr){
    for (r2 in 1:lr){
      tau_matrix[r1,r2] <- kendall(rankings[[r1]], rankings[[r2]])
    }
  }
  rnames <- sapply(rankings, function(x) paste(x@params$method, x@params$reg))
  if (!is.null(labels)){
   rnames <- labels
  }
  colnames(tau_matrix) <- rnames
  rownames(tau_matrix) <- rnames

  corrplot(tau_matrix, addCoef.col = addCoef.col, cl.pos = cl.pos, tl.col = "black", ...)
}
