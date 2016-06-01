#' @title Internal - Generatate marker expression scores describing phenotypes
#' 
#' @description 
#' This function is used internally to generate a melted numeric matrix of discrete expression scores for each marker of each cluster.
#' 
#' @details 
#' NA values are removed 
#' 
#' @param SPADEResults a SPADEResults object
#' @param num a numeric value specifying the number of markers expression categories
#'  
#' @return a numeric matrix of expression scores
#' 
#' @import gtools plyr
#' @importFrom plyr is.discrete
computePhenoTable <- function(SPADEResults, num = 5){

    data        <- SPADEResults@marker.expressions

    data        <- na.omit(data)# NA values are removed, generate a warning ?
    data.melted <- reshape2::melt(data, id.vars = c("sample", "cluster"))
    
    colnames(data.melted) <- c("sample", "cluster", "marker", "value")
    data.melted$marker    <- as.vector(data.melted$marker)
    means                 <- plyr::ddply(data.melted, c("cluster", "marker"), function(df){mean(df$value, na.rm = TRUE)}) #NA values are removed
    colnames(means)       <- c("cluster", "marker", "value")

    for(i in 1:nrow(means)){
        
        cluster <- means[i, "cluster"]
        value   <- means[i, "value"]

        min     <- SPADEResults@bounds[1, means[i, "marker"]]
        max     <- SPADEResults@bounds[2, means[i, "marker"]]

        seq     <- seq(from = min, to = max, length.out = num)
        means[i, "value"] <- which.min(abs(value - seq))
        
    }
    
    means <- means[gtools::mixedsort(colnames(means))]

    return(means)
    
}


#' @title Internal - Create a list of elements allowing to build a heatmap   
#'
#' @description 
#' This function is used internally to build the element needed for an heatmap
#' 
#' @param matrix a numeric matrix containing the markers expression categories
#' @param dendrogram.type a character specifying the look of dendrograms ("rectangle" or "triangle", "rectangle" by default)
#' @param num a numeric value specifying the number of markers expression categories
#' @param clustering.markers a character vector of clustering markers
#' @return a list of 3 plots (top dendrogram, right dendrogram, heatmap)
#'
#' @import ggplot2 reshape2 grDevices
ggheatmap <- function(matrix, dendrogram.type = "rectangle", num = 5, clustering.markers = NULL) {#TO ADD, dists = c("euclidian","euclidian")

    row.hc <- hclust(dist(matrix), "ward.D")
    col.hc <- hclust(dist(t(matrix)), "ward.D")
    
    row.dendro <- ggdendro::dendro_data(as.dendrogram(row.hc), type = dendrogram.type)
    col.dendro <- ggdendro::dendro_data(as.dendrogram(col.hc), type = dendrogram.type)
    
    col.plot <- g_dendro(col.dendro, col=TRUE)
    row.plot <- g_dendro(row.dendro, row=TRUE)
    
    col.ord <- match(col.dendro$labels$label, colnames(matrix))
    row.ord <- match(row.dendro$labels$label, rownames(matrix))
    
    mat.ordered  <- matrix[row.ord,col.ord]

    data.frame           <- as.data.frame(mat.ordered)
    data.frame$markers   <- rownames(mat.ordered)
    data.frame$markers   <- with(data.frame, factor(markers, levels = markers, ordered = TRUE))
    melted.data.frame    <- reshape2::melt(data.frame, id.vars = "markers")
        
    colfunc <- grDevices::colorRampPalette(c("#FFFFFF", "#ECE822", "#F9A22B", "#EE302D", "#A32D33"))#white -> yellow -> orange -> red -> brown

    melted.data.frame$value <- as.factor(melted.data.frame$value)
    
    centre.plot <- ggplot2::ggplot(melted.data.frame, ggplot2::aes_string(x = "variable", y = "markers")) + 
                   ggplot2::geom_tile(ggplot2::aes_string(fill = "value"), colour = "black") +
                   ggplot2::scale_fill_manual(values = colfunc(num), guide = ggplot2::guide_legend(title          = "relative.expression",
                                                                                                   direction      = "horizontal",
                                                                                                   ncol           = 5,
                                                                                                   byrow          = TRUE,
                                                                                                   label.theme    = ggplot2::element_text(size = 10, angle = 0), 
                                                                                                   label.position = "bottom",
                                                                                                   label.hjust    = 0.5,
                                                                                                   title.position = "top")) + 
                   ggplot2::theme(legend.text          = ggplot2::element_text(size = 4),
                                  panel.background     = ggplot2::element_rect("white"),
                                  axis.text.x          = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                                  legend.position      = c(ifelse(num >= 5, 0.6, 1 - (num * 0.1)), 0.5),
                                  legend.background    = ggplot2::element_blank())

    if (!is.null(clustering.markers)){
        clustering.markers <- is.element(data.frame$markers, clustering.markers)
        bold.markers       <- ifelse(clustering.markers,"bold","plain")
        centre.plot        <- centre.plot + ggplot2::theme(axis.text.y = ggplot2::element_text(face = bold.markers))#bold.markers
    }

    ret <- list(col = col.plot, row = row.plot, centre = centre.plot)

    return(ret)
    
}

#' @title Internal - Build a dendrograms plot
#'
#' @description 
#' This function is used internally to generate a 'ggplot' dendrogram.
#'
#' @details 
#' It is to note that 'row' and 'col' are mutuality excluded (both cannot be both TRUE) with priority to row.
#' 
#' @param dist a numeric matrix containing distances between objects
#' @param row a logical value specifying if the horizontal dendrogram must be computed
#' @param col a logical value specifying if the vertical dendrogram must be computed
#' 
#' @return a 'ggplot' dendrogram object
#'
#' @import ggplot2 ggdendro grid
g_dendro <- function(dist, row=!col, col=!row) {

    p <- ggplot2::ggplot() +
         ggplot2::geom_segment(data = ggdendro::segment(dist),
                               ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
         ggplot2::labs(x = NULL, y = NULL) +
         ggdendro::theme_dendro() +
         ggplot2::theme(axis.line        = ggplot2::element_blank(),
                        axis.text.x      = ggplot2::element_blank(),
                        axis.text.y      = ggplot2::element_blank(),
                        axis.ticks       = ggplot2::element_blank(),
                        axis.title.x     = ggplot2::element_blank(),
                        axis.title.y     = ggplot2::element_blank(),
                        legend.position  = "none",
                        panel.background = ggplot2::element_blank(),
                        panel.border     = ggplot2::element_blank(),
                        panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        plot.background  = ggplot2::element_blank(),
                        plot.margin      = grid::unit(c(0,0,0,0), "cm"),
                        panel.margin     = grid::unit(c(0,0,0,0), "cm"))
    if(row) {
        p <- p + ggplot2::scale_x_continuous(expand = c(0.005,0.005)) +
                 ggplot2::coord_flip()                
    } 
    else {
        p <- p +
                ggplot2::scale_x_continuous(expand = c(0.005,0.005))
    }
    return(p)
}

#' @title Internal - Extraction of ggplot legend
#'
#' @description 
#' This function is used internally to extract the legend from a 'ggplot' object.
#'
#' @param gplot a 'ggplot' plot
#' 
#' @return a 'ggplot' legend object
#'
#' @import ggplot2 gtable
g_legend <- function(gplot){
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(gplot))
    tmp <- gtable::gtable_filter(tmp, "guide-box")
    
    return(tmp$grobs[[TRUE]])
}

#' @title Internal - Extraction of ggplot axes
#'
#' @description 
#' This function is used internally to extract axes from a 'ggplot' objet.
#'
#' @details 
#' It is to note that 'x' and 'y' are mutuality excluded (both cannot be both TRUE) with priority to 'x'.
#'
#' @param gplot a 'ggplot' plot
#' @param x.axis a logical value specifying if the x-axis must be extract
#' @param y.axis a logical value specifying if the y-axis must be extract
#' 
#' @return a 'ggplot' axis object
#'
#' @import ggplot2 gtable
g_axis <- function(gplot, x.axis =! y.axis, y.axis =! x.axis ){
    if (x.axis){
        name <- "axis-b"
    }else{
        name <- "axis-l"
    }
    built <- ggplot2::ggplot_build(gplot)
    tmp   <- ggplot2::ggplot_gtable(built)
    tmp   <- gtable::gtable_filter(tmp, name)
    return(tmp$grobs[[TRUE]])
}

#' @title Internal - Generate an heatmap by assembling elements
#'
#' @description 
#' This function is used internally to displays the heatmap elements build by 'ggheatmap()'
#'
#' @param list the list of ggplot object provided by ggheatmap
#' @param col.width size of horizontal dendrogram
#' @param row.width size of vertical dendrogram
#' 
#' @return a ggplot2 axis
#'
#' @import ggplot2 gridExtra grid
ggheatmap.plot <- function(list, col.width=0.15, row.width=0.15) {

    layout <- rbind(c(2, 1, NA),
                    c(5, 3, 4),
                    c(NA, 6, NA))

    legend <- g_legend(list$centre)
    x.axis <- g_axis(list$centre, x.axis = TRUE)        
    y.axis <- g_axis(list$centre, y.axis = TRUE)
    
    center.without_legend = list$centre + ggplot2::theme(axis.line        = ggplot2::element_blank(),
                                                         axis.text.x      = ggplot2::element_blank(),
                                                         axis.text.y      = ggplot2::element_blank(),
                                                         axis.ticks       = ggplot2::element_blank(),
                                                         axis.title.x     = ggplot2::element_blank(),
                                                         axis.title.y     = ggplot2::element_blank(),
                                                         legend.position  = "none",
                                                         panel.background = ggplot2::element_blank(),
                                                         panel.border     = ggplot2::element_blank(),
                                                         panel.grid.major = ggplot2::element_blank(),
                                                         panel.grid.minor = ggplot2::element_blank(),
                                                         plot.background  = ggplot2::element_blank(),
                                                         plot.margin      = grid::unit(c(0,0,0,0), "cm"),
                                                         panel.margin     = grid::unit(c(0,0,0,0), "cm"))
                                             
    ret <- gridExtra::arrangeGrob(list$col, #1 on the layout
                                  legend, #2 on the layout
                                  center.without_legend, #3 on the layout
                                  list$row, #4 on the layout
                                  y.axis, #5 on the layout
                                  x.axis, #6 on the layout
                                  layout_matrix = layout,
                                  widths        = grid::unit(c(col.width, 1-(2*col.width), col.width), "null"),
                                  heights       = grid::unit(c(row.width, 1-(2*row.width), row.width), "null"),
                                  top           = "Heatmap Viewer")

    return(ret)
}

#' @title Internal - Generate a circle representation
#'
#' @description 
#' This function is used internally to generate a packed circles representation
#'
#' @param circles a 2 column dataframe the clusters to be displayed and theirs sizes
#' @param class a numeric specifyng the class number to be displayed
#' @param color a character specifying the color of the packed circles representation
#' @param npoint a numeric specifying the levels of details of polygones
#' @param limits a numeric specifying the size of the coordinate system centered on (0,0)
#'
#' @return a ggplot2 object
#'
#' @import ggplot2 ggrepel gridExtra grid packcircles
buildCircles <- function(circles,
                         color   = "grey80",
                         class   = NA,
                         npoint  = 100,
                         limits  = 30000,
                         maxiter = 100) {

    xyr <- data.frame(x = runif(nrow(circles), 0, 1),
                      y = runif(nrow(circles), 0, 1),
                      r = circles$size)

    res  <- packcircles::circleLayout(xyr, xlim = c( - limits, limits), ylim = c( - limits, limits), maxiter = 1000, wrap = FALSE)
    data <- packcircles::circlePlotData(layout = res$layout, npoints = npoint)
    text <- cbind(res$layout, cluster = circles$cluster)
    
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(paste0("Class ", class)) +
            ggplot2::geom_polygon(ggplot2::aes_string(x = "x", y = "y", group = "id"),
                                  fill = color,
                                  color = "grey80",
                                  alpha = 0.2) +
            ggrepel::geom_text_repel(data = text, ggplot2::aes_string(x = "x", y = "y", label = "cluster"), size = 3,
                                     box.padding = grid::unit(0.35, "lines"),
                                     point.padding = grid::unit(0.3, "lines")) +
            ggplot2::coord_equal(xlim = c( - limits, limits), ylim = c( - limits, limits)) +
            ggplot2::theme(axis.line = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           legend.position = "none",
                           panel.background = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background = ggplot2::element_blank(),
                           plot.margin = grid::unit(c(0, 0, 0, 0), "cm"),
                           panel.margin = grid::unit(c(0, 0, 0, 0), "cm"))

    return(plot)

}

#' @title Internal - Generate an legend for circles representation
#'
#' @description 
#' This function is used internally to generate the legend of a packed circles representation
#'
#' @param circles a 3 colmuns data frame with the x, y coordinate of points and their raduis
#' @param npoint a numeric specifying the levels of details of polygones
#' @param limits a numeric specifying the size of the coordinate system centered on (0,0)
#' 
#' @return a ggplot2 object
#'
#' @import ggplot2 ggrepel gridExtra grid packcircles
buildCirclesLegend <- function(circles = data.frame(x = c(-29500, -19000, -8000, 3000, 20000),
                                                    y = c(20000, 20000, 20000, 20000, 20000),
                                                    r = c(500, 1000, 2000, 5000, 10000)),
                               npoint  = 100,
                               limits  = 30000) {

    data <- packcircles::circlePlotData(layout = circles, npoints = npoint)
    text <- circles
    colnames(text) <- c("x", "y", "cluster")

    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle("number.of.cells") +
            ggplot2::geom_polygon(ggplot2::aes_string(x = "x", y = "y", group = "id"),
                                  fill  = "white",
                                  color = "grey80",
                                  alpha = 0.2) +
            ggplot2::geom_text(data = text, ggplot2::aes_string(x = "x", y = "y-20000", label = "cluster"), size = 3) +
            ggplot2::coord_equal(xlim = c( - limits, limits), ylim = c( - limits, limits)) +
            ggplot2::theme(axis.line        = ggplot2::element_blank(),
                           axis.text.x      = ggplot2::element_blank(),
                           axis.text.y      = ggplot2::element_blank(),
                           axis.ticks       = ggplot2::element_blank(),
                           axis.title.x     = ggplot2::element_blank(),
                           axis.title.y     = ggplot2::element_blank(),
                           legend.position  = "none",
                           panel.background = ggplot2::element_blank(),
                           panel.border     = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background  = ggplot2::element_blank(),
                           plot.margin      = grid::unit(c(0, 0, 0, 0), "cm"),
                           panel.margin     = grid::unit(c(0, 0, 0, 0), "cm"))

    return(plot)
                       
}