## Individual gene plots with loess
## reg_gam
#' @title Plot one gene raw data and its regressions
#' @name reg_gam
#'
#' @description Tools for visualizing gene signals for one given gene with the gam framework.
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param reg.f a function to perform regression, either "ns" for natural splines, "loess" or "s" (default is "loess").
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param regression logical, if the loess regression is to be computed and plotted or not (default is TRUE).
#' @param null.model logical, if the null model is to be computed and plotted or not (default is TRUE).
#' @param npred logical, if the unshared part of the data is to be plotted or not(default is FALSE).
#' @param sd.show logical, if the plot of standard deviation is wanted (default is FALSE).
#' @param legend.show logical, if the legend is wanted (default is FALSE).
#'
#' @return returns \itemize{\item{\code{pl},the visualization of the data}
#' \item{\code{regs}, the regression objects for both lineages and the null model.}}
#'
#' @import ggplot2
#' @importFrom splines ns
#' @importFrom gam gam
#' @import gam
#' @importFrom gam lo
#' @export
reg_gam <- function(data,
                    gene,
                    reg.f = "loess",
                    span = 0.75,
                    s.df = 4,
                    regression = T,
                    null.model = T,
                    npred = F,
                    sd.show = F,
                    legend.show = F){
  d <- s.df
  t = data@t
  w = data@w
  y = log1p(data@counts[gene,])
  # discard NAS from analysis
  w1 <- w[w[,1]!=0,1]
  w2 <- w[w[,2]!=0,2]
  l1_cells <- names(w1)
  l2_cells <- names(w2)
  y1 <- y[l1_cells]
  y2 <- y[l2_cells]
  t1 <- t[l1_cells,1]
  t2 <- t[l2_cells,2]
  reg.df1 <- data.frame(y.fit = y1, x.fit = t1, w.fit = w1, lineage = rep(1, length(y1)))
  reg.df2 <- data.frame(y.fit = y2, x.fit = t2, w.fit = w2, lineage = rep(2, length(y2)))
  #time to predict the new data
  t1new <- seq(0, max(t1), length.out = length(l1_cells))
  t2new <- seq(0,max(t2), length.out = length(l2_cells))
  pl <- ggplot() +
    geom_point(aes(t2, y2), col = "#377EB8", alpha = w2, shape = 1) +
    geom_point(aes(t1, y1), col = "#E41A1C", alpha = w1, shape = 1) +
    ggtitle(gene, paste(reg.f, "regression")) + coord_cartesian(ylim=c(-1.5,10))+ xlab("times") + ylab("logcounts")
  if (regression == F){
    return(pl)
  }
  else{
    reg.df.d <- rbind(reg.df1, reg.df2)
    if (reg.f == "loess"){
      alt <- gam(y.fit ~ lineage + lo(x.fit, span = span) + lo(x.fit, span = span):lineage , weights = w.fit, data = reg.df.d)
    }
    else if (reg.f == "ns"){
      alt <- gam(y.fit ~ lineage + ns(x.fit, df = s.df) + ns(x.fit, df = s.df):lineage , weights = w.fit, data = reg.df.d)
    }
    if (reg.f == "s"){
      alt <- gam(y.fit ~ lineage + s(x.fit, df = d) + s(x.fit, df = d):lineage , weights = w.fit, data = reg.df.d)
    }
    y_pred.alt1 <- predict(alt, data.frame(x.fit = t1new, lineage = rep(1, length(t1new))))
    y_pred.alt2 <- predict(alt, data.frame(x.fit = t2new, lineage = rep(2, length(t2new))))
    plalt <- pl + geom_line(aes(t1new, y_pred.alt1), col = "#E41A1C") + geom_line(aes(t2new, y_pred.alt2), col = "#377EB8")
    regs <- list(alt = alt)
    if (null.model == T){
      # null model
      if (reg.f == "loess"){
        null.m <- gam(y.fit ~ lo(x.fit, span = span), weights = w.fit, data = reg.df.d)
      }
      else if (reg.f == "ns"){
        null.m <- gam(y.fit ~ ns(x.fit, df = s.df), weights = w.fit, data = reg.df.d)
      }
      if (reg.f == "s"){
        null.m <- gam(y.fit ~  s(x.fit, df = d), weights = w.fit, data = reg.df.d)
      }
      y_pred.null <- predict(null.m, data.frame(x.fit = t1new))
      plalt <- plalt + geom_line(aes(t1new, y_pred.null, colour = "null model"), linetype = 2, na.rm = T)
      regs$null = null.m
    }
  }
  if (sd.show == T){
      #TO ADD
      # ## confidence intervals:
      # se1 <- loess.sd(x = t[,1], y = y, weights = w[,1], nsigma = 1, na.action = na.exclude)
      # se2 <- loess.sd(x = t[,2], y = y, weights = w[,2], nsigma = 1, na.action = na.exclude)
      #
      # plalt <- plalt + geom_ribbon(aes(x = se1$x, ymin = se1$lower, ymax = se1$upper), alpha = 0.2, fill = "#E41A1C") +
      #   geom_ribbon(aes(x = se2$x, ymin = se2$lower, ymax = se2$upper), alpha = 0.2, fill = "#377EB8")
  }
  if (npred == T){
      #discard unappropriate cells for new predictions (keep only w >0.5 approximately)
      middle_w1 <- w[w[,1]!=0 & w[,1]!=1,1]
      middle_w2 <- w[w[,2]!=0 & w[,2]!=1,2]
      cells_pred1 <- names(w[w[,1] > (0.5 + sd(middle_w1)),1])
      cells_pred2 <- names(w[w[,2] > (0.5 + sd(middle_w2)),2])
      plalt <- plalt + geom_point(aes(t[cells_pred1,1], predict(regs$alt, data.frame(x.fit = t[cells_pred1,1], lineage = rep(1, length(cells_pred1))))), alpha = 0.4, col = "black", size = 1)+
        geom_point(aes(t[cells_pred2,2], predict(regs$alt, data.frame(x.fit = t[cells_pred2,2], lineage = rep(2, length(cells_pred2))))), alpha = 0.4, col = "black", size = 1)
  }
  if (legend.show == F){
      plalt <- plalt + theme(legend.position = "none")
  }
  res <- list(pl = plalt, reg = regs)
  return(res)
}

#' @title Plot several genes expression patterns
#' @name plot_multigenes
#'
#' @description Tools for visualizing gene signals for sveral ranked genes and display their ranking informations.
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param ranking a \code{rankingDE} object, rankings rows of the genes of interest in the ranking dataframe.
#' @param subset.genes character vector, let the user specifies the names of the genes of interest (default is NULL).
#' @param quick.view character, "first" (default), "last" : specify if we want to see the 8 first or end genes of the ranking.
#' @param reg.f a function to perform regression, either "ns" for natural splines, "loess" or "splines" (default is "loess").
#' @param span numeric, a smoothing parameter for the regression function (default is 0.75, see \code{gam::lo} for details).
#' @param s.df numeric, a smoothing parameter for the nsplines regregression (default is 4, see \code{splines::s} for details about regularization).
#' @param null.model logical, if the plot of null model is wanted (default is FALSE).
#' @param grid.size 2 by 2 vector, for the number of rows and the number of columns that we want for the plot (the size of grid must #' be greater than the number of genes of interest), default is NULL fo a squared grid.
#' @return a visualization of the genes of interest.
#'
#' @import ggplot2
#' @import cowplot
#' @export


plot_multigenes <- function(data,
                            ranking,
                            subset.genes = NULL,
                            quick.view = "first",
                            reg.f = "loess",
                            span = 0.75,
                            s.df = 4,
                            null.model = F,
                            grid.size = NULL
){
  if (is.null(subset.genes)){
    if (quick.view == "first"){
      subset.genes <- rownames(ranking@ranking.df)[1:8]
    }
    else if (quick.view == "last"){
      subset.genes <- tail(rownames(ranking@ranking.df),8)
    }
    if (is.null(grid.size)){
      grid.size <- c(2,4)
    }
  }
  if (is.null(grid.size)){
    c <- ceiling(sqrt(length(subset.genes)))
    grid.size <- c(c,c)
  }
  method = ranking@params$method

  subset.genes <- data.frame(ranking@ranking.df)[subset.genes,]
  graphs <- lapply(rownames(subset.genes), function(x) reg_gam(gene = x, data = data, reg.f = reg.f, span = span, s.df =s.df, null.model = null.model)$pl + labs(subtitle = paste0(method,".dist: ",round(subset.genes[x,1],1), " | ",method,".rank:", subset.genes[x,2])))
  return(plot_grid(plotlist = graphs, ncol = grid.size[2], nrow = grid.size[1]))
}

#' @title Differential Expression analysis according with time
#' @name timeDE
#'
#' @description Tools for visualizing gene signals in time for one lineage analysis.
#'
#' @param data a \code{lineageDEDataSet} with results to be plotted.
#' @param gene character, a gene of interest.
#' @param lineage integer, if working with multi lineages data, it sspecifies the lineage of interest).

#'
#' @return returns \itemize{\item{\code{pl},the visualization of the data}
#' \item{\code{reg}, the regression objects for the lineage of interest and the null model (flat line).}}
#'
#' @import ggplot2
#'
#' @import gam
#' @export
timeDE <- function(data,
                   gene,
                   lineage,
                   span = 0.5){
  t = data@t
  w = data@w
  y = log1p(data@counts[gene,])
  # discard NAS from analysis
  w <- w[w[,lineage]!=0,lineage]
  l_cells <- names(w)
  y <- y[l_cells]
  t <- t[l_cells,lineage]
  reg.df <- data.frame(y.fit = y, x.fit = t, w.fit = w)
  #time to predict the new data
  tnew <- seq(0, max(t), length.out = length(l_cells))
  pl <- ggplot() +
    geom_point(aes(t, y), shape = 1, alpha = w) + ggtitle(gene, "single lineage DE") +
    coord_cartesian(ylim=c(-1.5,10))+ xlab("times") + ylab("logcounts")
  alt <- gam(y.fit ~ lo(x.fit, span = span), weights = w.fit, data = reg.df)
  null <- mean(y)
  reg <- list(alt = alt, null = null)
  y_pred.alt <- predict(alt, data.frame(x.fit = tnew))
  plalt <- pl + geom_line(aes(tnew, y_pred.alt), col = "#E41A1C") + geom_hline(yintercept = null, col = "#377EB8")
  return(list(pl = plalt, reg = reg, pval = summary(reg$alt)[4][[1]][1, 5]))
}



#' ## Individual gene plots with vgam
#' #' @title Plot one gene raw data and its vgam regressions
#' #' @name reg_vgam
#' #'
#' #' @description Tools for visualizing gene signals for one given gene and compute vgam regression on the count data.
#' #'
#' #' @param data a \code{lineageDEDataSet} with results to be plotted.
#' #' @param gene character, a gene of interest.
#' #' @param model which disrtibution assumption is made on the residuals; either "negbinomial" or  gausian (default is gaussian).
#' #' @param regression logical, if the loess regression is to be computed and plotted or not (default is TRUE).
#' #' @param npred logical, if the unshared part of the data is to be plotted or not(default is FALSE).
#' #' @param null.model logical, if the null model is to be computed and plotted or not (default is TRUE).
#' #' @param legend.show logical, if the legend is to be shown or not (default is FALSE).
#' #'
#' #' @return returns \itemize{\item{\code{pl},the visualization of the data}
#' #' \item{\code{reg}, the vgam regression objects for both lineages and the null model.}}
#' #'
#' #' @import ggplot2
#' #' @import VGAM
#' #' @export
#'
#' reg_vgam <- function(data,
#'                      gene,
#'                      model = "gaussian",
#'                      regression = T,
#'                      npred = F,
#'                      null.model = T,
#'                      legend.show = F){
#'   # choice of the model
#'   if (model %in% "negbinomial"){
#'     y <- round(data@counts[gene, ],0)
#'     y[y<0] <- 0
#'     fam <- negbinomial(imu = mean(y), isize = 1)
#'   }
#'   else{
#'     fam = gaussianff()
#'     y <- log1p(data@counts[gene,])
#'   }
#'   w <- data@w
#'   t <- data@t
#'
#'   # discard NAS from analysis
#'   w1 <- w[w[,1]!=0,1]
#'   w2 <- w[w[,2]!=0,2]
#'   l1_cells <- names(w1)
#'   l2_cells <- names(w2)
#'   y1 <- y[l1_cells]
#'   y2 <- y[l2_cells]
#'   t1 <- t[l1_cells,1]
#'   t2 <- t[l2_cells,2]
#'   reg.df1 <- data.frame(y.fit = y1, x.fit = t1, w.fit = w1, lineage = rep(1, length(y1)))
#'   reg.df2 <- data.frame(y.fit = y2, x.fit = t2, w.fit = w2, lineage = rep(2, length(y2)))
#'   #time to predict the new data
#'   t1new <- seq(0, max(t[,1], na.rm = T), length.out = length(l1_cells))
#'   t2new <- seq(0,max(t[,2], na.rm = T), length.out = length(l2_cells))
#'
#'   pl <- ggplot() +
#'     geom_point(aes(t2, y2), col = "#377EB8", alpha = w2, shape = 1) +
#'     geom_point(aes(t1, y1), col = "#E41A1C", alpha = w1, shape = 1) +
#'     ggtitle(gene, "splines regression") + coord_cartesian(ylim=c(-1.5,10))+ xlab("times") + ylab("logcounts")
#'   if (regression == F){
#'     return(pl)
#'   }
#'   else{
#'     reg.df.d <- rbind(reg.df1, reg.df2)
#'     spl.alt <- vglm(y.fit ~ lineage + sm.s(x.fit, df = 3) + sm.s(x.fit, df = 3):lineage, family = fam, weights = w.fit, data = reg.df.d)
#'     reg <- list(alt = spl.alt)
#'     if (null.model == T){
#'       # null model
#'       spl.d <- vglm(y.fit ~ sm.s(x.fit, df = 3), weights = w.fit, family = fam, data = reg.df.d)
#'       reg$null <- spl.d
#'     }
#'     # prediciton at new points and curve plotting
#'     if (model %in% "negbinomial"){
#'       y_pred.alt1 <- predict(spl.alt, data.frame(x.fit = t1new, lineage = rep(1, length(t1new))))[,1]
#'       y_pred.alt2 <- predict(spl.alt, data.frame(x.fit = t2new, lineage = rep(2, length(t2new))))[,1]
#'       if (null.model ==T){
#'         y_pred.d <- predict(spl.d, data.frame(x.fit = t1new))[,1]
#'       }
#'     }
#'     else{
#'       y_pred.alt1 <- predict(spl.alt, data.frame(x.fit = t1new, lineage = rep(1, length(t1new))))
#'       y_pred.alt2 <- predict(spl.alt, data.frame(x.fit = t2new, lineage = rep(2, length(t2new))))
#'       if (null.model == T) {
#'         y_pred.d <- predict(spl.d, data.frame(x.fit = t1new))
#'       }
#'     }
#'     plalt <- pl + geom_line(aes(t1new, y_pred.alt1), col = "#E41A1C") + geom_line(aes(t2new, y_pred.alt2), col = "#377EB8" )
#'     if (null.model == T){
#'       # null model
#'       plalt <- plalt + geom_line(aes(t1new, y_pred.d, col = "null model"), linetype = 2, na.rm =T)
#'     }
#'     if (npred == T){
#'       #discard shared cells for new predictions (keep only w >0.5 approximately)
#'       middle_w1 <- w[w[,1]!=0 & w[,1]!=1,1]
#'       middle_w2 <- w[w[,2]!=0 & w[,2]!=1,2]
#'       cells_pred1 <- names(w[w[,1] > (0.5 + sd(middle_w1)),1])
#'       cells_pred2 <- names(w[w[,2] > (0.5 + sd(middle_w2)),2])
#'       plalt <- plalt + geom_point(aes(t[cells_pred1,1], predict(reg$alt, data.frame(x.fit = t[cells_pred1,1], lineage = rep(1, length(cells_pred1))))), alpha = 0.4, col = "#E41A1C", size = 1)+
#'         geom_point(aes(t[cells_pred2,2], predict(reg$alt, data.frame(x.fit = t[cells_pred2,2], lineage = rep(2, length(cells_pred2))))), alpha = 0.4, col = "#377EB8", size = 1)
#'     }
#'     if (legend.show == F) {
#'       plalt <- plalt + theme(legend.position = "none")
#'     }
#'   }
#'   res <- list(pl = plalt, reg = reg)
#'   return(res)
#' }
