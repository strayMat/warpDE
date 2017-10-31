#' @title Interactive plot for a ranking
#' @name warpDE_app
#'
#' @description Given a ranking of genes of a scRNA-Seq experiment, open an interactive interface of the best n genes given the ranking (where n is interactively chosen by the user).
#' @param ranking a ranking of the gens of interest.
#' @param data a \code{warpDEDataSet} with genes as rows and cells as columns.
#'
#' @return returns nothing.
#'
#' @import shiny
#' @import plotly
# @importFrom plotly plot_ly
# @importFrom plotly add_markers
# @importFrom plotly layout
# @importFrom plotly plotlyOutput
# @importFrom plotly renderPlotly
#' @export

warpDE_app <- function(data, ranking){
  genes <- rownames(ranking@ranking.df)
  nb_genes <- 30
  palette3d <- c("#C78C0A",  "#D1B00A", "#BAB700","#0CC700", "#75D100")

  ui <- fluidPage(
    titlePanel("Interactive 3D plots for gene expression"),
    sidebarPanel(
      selectInput("gene.choice", label = "gene choice", choices = genes[1:nb_genes])
    ),
    plotlyOutput("plot")
  )

  server <- function(input, output) {
    # compute the 3d low dimensionnal representation of the cells.
    pca <- prcomp(t(log1p(df@counts)))
    output$plot <- renderPlotly({
      g <- input$gene.choice
      gene.expression <- log1p(data@counts[g,])
      g.df <- data.frame(pc1 = pca$x[,1], pc2 = pca$x[,2],pc3 = pca$x[,3], col = gene.expression)

      plot_ly(g.df, x = ~pc1,  y = ~pc2, z = ~pc3,
              mode = "marker",
              marker = list(size = 4, opacity = 1), color = ~col, colors = c('#474647', '#1CD600')) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'PC 1'),
                            yaxis = list(title = 'PC 2'),
                            zaxis = list(title = 'PC 3')),
               annotations = list(
                 x = 1.13,
                 y = 1.05,
                 text = paste(g,'log_expression'),
                 xref = 'paper',
                 yref = 'paper',
                 showarrow = FALSE
               ))
    })
  }

  shinyApp(ui, server)

}
