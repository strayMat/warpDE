
## Individual gene plots
## plot_gene
#' @title Plot one gene raw data and its regressions
#' @name plot_gene
#'
#' @description Tools for visualizing gene signals for one given gene.
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param regression logical, if the regression is to be computed and plotted or not (defualt is TRUE).
#' @param MD logical, if the null model is to be computed and plotted or not (default is TRUE).
#' @param span numeric, a span parameter for the loess regression (default is 0.75, same as the \code{loess} function).
#' @param npred logical, if the unshared part of the data is to be plotted or not(default is FALSE).
#' @param sd.show logical, if the plot of standard deviation is wanted (default is FALSE).
#'
#' @return returns \itemize{\item{\code{pl},the visualization of the data}
#' \item{\code{reg}, the regression objects for both lineages and the null model.}}
#'
#' @import ggplot2
#' @importFrom msir loess.sd
#' @export
plot_gene <- function(data,
                      gene,
                      regression = T,
                      MD = T,
                      span = 0.75,
                      npred = F,
                      sd.show = F,
                      legend.show = F
                      ){
          t = data@t
          w = data@w
          prediction_length <- 100
          #time to predict the new data
          t1new <- seq(0, max(t[,1], na.rm = T), length.out = prediction_length)
          t2new <- seq(0,max(t[,2], na.rm = T), length.out = prediction_length)
          y <- data@logCounts[gene,]

          pl <- ggplot() +
            geom_point(aes(t[,2], y), col = "blue", shape = 2, alpha = w[,2], na.rm =T) +
            geom_point(aes(t[,1], y), col = "red", shape = 1, alpha = w[,1], na.rm =T) +
            ggtitle(gene, "(loess)") + coord_cartesian(ylim=c(-1.5,10))+ xlab("times") + ylab("logcounts")

          if (regression == F){
            return(pl)
          }
          else{
            #alternative model
            lo1 <- loess(y ~ t[,1], weights = w[,1], span = span)
            lo2 <- loess(y ~ t[,2], weights = w[,2], span = span)
            # prediciton at new points and curve plotting
            y_pred1 <- predict(lo1, t1new)
            y_pred2 <- predict(lo2, t2new)
            plalt <- pl +
              geom_line(aes(t1new, y_pred1), color = "red", linetype = 1, na.rm =T) +
              geom_line(aes(t2new, y_pred2), color = "blue", linetype = 1, na.rm =T)

            reg <- list(lo1 = lo1, lo2 = lo2)

             # null models
            if (MD == T){
              lodouble <- loess(c(y, y) ~ c(t[,1], t[,2]), weights = c(w[,1], w[,2]), span = span, na.action = na.exclude)
              plalt <- plalt + geom_line(aes(t1new, predict(lodouble, t1new), colour = "null model"), linetype = 2, na.rm = T)
              reg$lod <- lodouble
            }
            if (sd.show == T){
              ## confidence intervals:
              se1 <- loess.sd(x = t[,1], y = y, weights = w[,1], span = span, nsigma = 1, na.action = na.exclude)
              se2 <- loess.sd(x = t[,2], y = y, weights = w[,2], span = span, nsigma = 1, na.action = na.exclude)

              plalt <- plalt + geom_ribbon(aes(x = se1$x, ymin = se1$lower, ymax = se1$upper), alpha = 0.2, fill = "red") +
                geom_ribbon(aes(x = se2$x, ymin = se2$lower, ymax = se2$upper), alpha = 0.2, fill = "blue")
            }
            if (npred == T){
              #discard unappropriate cells for new predictions (keep only w >0.5 approximately)
              middle_w1 <- w[w[,1]!=0 & w[,1]!=1,1]
              middle_w2 <- w[w[,2]!=0 & w[,2]!=1,2]
              cells_pred1 <- names(w[w[,1] > (0.5 + sd(middle_w1)),1])
              cells_pred2 <- names(w[w[,2] > (0.5 + sd(middle_w2)),2])
              plalt <- plalt + geom_point(aes(t[cells_pred1,1], predict(reg$lo1, t[cells_pred1,1])), alpha = 0.4)+
                geom_point(aes(t[cells_pred2,2], predict(reg$lo2, t[cells_pred2,2])), alpha = 0.4)
            }
            if (legend.show == F){
              plalt <- plalt + theme(legend.position = "none")
            }
          }
          return(res <- list(pl = plalt, reg = reg))
}


## multi gene plots fior ranking objects
#' @title Plot several genes expression patterns
#' @name plot_gene
#'
#' @description Tools for visualizing gene signals for sveral ranked genes and display their ranking informations.
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param ranking a \code{rankingDE} object, rankings rows of the genes of interest in the ranking dataframe.
#' @param subset.genes character vector, the names of the genes of interest.
#' @param MD logical, if the plot of null model is wanted (default is FALSE).
#' @param grid.size 2 by 2 vector, for the number of rows and the number of columns that we want for the plot (the size of grid must #' be greater than the number of genes of interest), default is NULL fo a squared grid.
#' @return a visualization of the genes of interest.
#'
#' @import ggplot2
#' @import cowplot
#' @export

plot_multigenes <- function(data,
                            ranking,
                            subset.genes,
                            MD = F,
                            grid.size = NULL
                            ){
  if (is.null(grid.size)){
    c <- ceiling(sqrt(length(subset.genes)))
    grid.size <- c(c,c)
  }
  method = ""
  if (grepl("dtw", ranking@params$method)){method = "dtw"}
  if (grepl("likelihood", ranking@params$method)){method = "lkl"}
  subset.genes <- data.frame(ranking@ranking.df)[subset.genes,]
    graphs <- lapply(rownames(subset.genes), function(x) plot_gene(gene = x, data = data, MD = MD)$pl +
      labs(subtitle = paste0(method,".dist: ",round(subset.genes[x,1],1), " | ",method,".rank:", subset.genes[x,2])))
    return(plot_grid(plotlist = graphs, ncol = grid.size[2], nrow = grid.size[1]))
}

