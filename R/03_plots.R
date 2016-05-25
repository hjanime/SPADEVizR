#' @title Visualization of abundant clusters
#'
#' @description 
#' Generates a scatter plot representation showing for each cluster its mean abundance and associated p-value.
#' 
#' @details 
#' By default, only significant abundant clusters are labeled. Labels for all clusters can be displayed by setting the 'show.all_labels' parameter to TRUE. 
#' 
#' @param AC an object of class AC (object returned by the 'computeAC()' function)
#' @param show.cluster_sizes a logical specifying if dot sizes are proportional to cell counts
#' @param show.all_labels a logical specifying if all cluster labels must be shown (or just significant clusters)
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggrepel grid
#' 
#' @export
abundantClustersViewer <- function(AC,
                                   show.cluster_sizes = TRUE,
                                   show.all_labels    = FALSE,
                                   show.on_device     = TRUE) {

    AC@result <- cbind (AC@result, cluster.size = AC@cluster.size)
    data.text <- AC@result

    if (!show.all_labels){
        data.text <- subset(AC@result, AC@result$significant)
    }

    plot <-  ggplot2::ggplot(data = AC@result) +
             ggplot2::ggtitle(paste0("Cell abundances of clusters (", format(sum(AC@cluster.size), big.mark=" "), " cells)", sep = "")) +
             ggplot2::geom_hline(yintercept = AC@th.mean,
                                 linetype   = "dashed",
                                 alpha      = 0.3,
                                 color      = "red",
                                 size       = 1) +  
             ggplot2::geom_vline(xintercept = -log10(AC@th.pvalue),
                                 linetype   = "dashed",
                                 alpha      = 0.3,
                                 color      = "red",
                                 size       = 1)

    if (show.cluster_sizes){
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significant", size = "cluster.size"),
                                           shape = 21,
                                           colour = "black",
                                           stroke = 1) +
                       ggplot2::scale_size(name = "number.of.cells", guide = ggplot2::guide_legend(order = 2))
    }else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significant"),
                                           shape = 21,
                                           colour = "black",
                                           stroke = 1)
    }                  
    
    x.max    <- ceiling(max(-log10(AC@th.pvalue), -log10(AC@result$pvalue)))
    x.breaks <- c(seq(0, x.max, by = 1), round(-log10(AC@th.pvalue), 2))
    
    y.max    <- ceiling(max(AC@th.mean, AC@result$mean))
    y.breaks <- seq(0, y.max, by = 1)

    plot <- plot +  ggrepel::geom_text_repel(data          = data.text, 
                                             ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", label = "cluster"),
                                             size          = 3,
                                             box.padding   = grid::unit(0.35, "lines"),
                                             point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey", "red"), guide = ggplot2::guide_legend(order = 1)) +
                    ggplot2::scale_x_continuous(limits = c(0, x.max), minor_breaks = NULL, breaks = x.breaks) + 
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab("-log10(p-value)") +
                    ggplot2::ylab(ifelse(AC@use.percentages, "mean (% of cells)", "mean (# of cells)")) +                  
                    ggplot2::theme_bw()
    
    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)

}

#' @title Visualization of differentially enriched clusters
#'
#' @description 
#' Generates a Volcano plot representation showing for each cluster
#' 
#' @details 
#' By default, only significant differentially enriched clusters are labeled. Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param DEC an object of class 'DEC' (object returned by the 'computeDEC()' function)
#' @param fc.log2 a logical specifying if fold-change or log2(fold-change) is use 
#' @param show.cluster.sizes a logical specifying if dot sizes are proportional to cell counts
#' @param show.all_labels a logical specifying if all cluster labels must be show or just significant cluster
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#'  
#' @import ggplot2 grid
#' 
#' @export
volcanoViewer <- function(DEC                = NULL,
                          fc.log2            = TRUE,
                          show.cluster.sizes = TRUE,
                          show.all_labels    = FALSE,
                          show.on_device     = TRUE) {
    
    th.fc <- DEC@th.fc
    
    if(fc.log2){
        th.fc <- log2(th.fc)
        temp <- DEC@result$fold.change
        for(i in 1:nrow(DEC@result)){
            DEC@result$fold.change[i] <- ifelse (temp[i] > 0, log2(DEC@result$fold.change[i]), -log2(abs(DEC@result$fold.change[i])))
        }
    }
    
    DEC@result <- cbind (DEC@result, cluster.size = DEC@cluster.size)
    
    data.text <- DEC@result
    if (!show.all_labels){
        data.text <- subset(DEC@result, DEC@result$significant)
    }
    x.min    <- floor(min(DEC@result$fold.change))
    x.max    <- ceiling(max(DEC@result$fold.change))
    x.max    <- max(x.max, abs(x.min))
    x.breaks <- c(round(c(-th.fc, th.fc), 2), seq(-x.max, x.max, by = 1))

    y.max    <- ceiling(max(-log10(DEC@th.pvalue), -log10(DEC@result$pvalue)))
    y.breaks <- c(seq(0, y.max, by = 1), round(-log10(DEC@th.pvalue), 2))

    title.details <- ifelse(DEC@use.percentages, "using % of cells", "using # of cells")
        
    plot <- ggplot2::ggplot(data = DEC@result, ggplot2::aes_string(x = "fold.change", y = "-log10(pvalue)")) +
            ggplot2::ggtitle(paste0("Volcano plot showing differentially enriched clusters ", title.details, " (", format(sum(DEC@cluster.size), big.mark=" "), " cells)", sep = "")) +
            ggplot2::geom_vline(xintercept = c(th.fc, -th.fc),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1) +  
            ggplot2::geom_hline(yintercept = -log10(DEC@th.pvalue),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1)
    if (show.cluster.sizes){
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significant", size = "cluster.size"),
                                           shape = 21,
                                           colour = "black",
                                           stroke = 1) +
                       ggplot2::scale_size(name = "number.of.cells")
    }else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significant"),
                                           shape = 21,
                                           colour = "black", 
                                           stroke = 1)
    }
                 
    plot <- plot + ggrepel::geom_text_repel(data          = data.text,
                                            ggplot2::aes_string(label = "cluster"),
                                            size          = 3,
                                            box.padding   = grid::unit(0.35, "lines"),
                                            point.padding = grid::unit(0.3, "lines")) +
                   ggplot2::scale_fill_manual(values = c("grey", "red")) +
                   ggplot2::scale_x_continuous(limits = c(-x.max, x.max), minor_breaks = NULL, breaks = x.breaks) +
                   ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                   ggplot2::xlab(paste0(ifelse(fc.log2,"log2(fold.change)","fold.change"),"\ncond2 < enriched > cond1")) +
                   ggplot2::ylab("-log10(p-value)") +
                   ggplot2::theme_bw()
    
    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}

#' @title Visualization of correlated clusters
#'
#' @description 
#' Generate a scatter plot representation showing for each cluster
#' 
#' @details 
#' By default, only significant correlated clusters are labeled. Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param CC an object of class 'CC' (object returned by the 'computeCC()' function)
#' @param show.cluster.sizes a logical specifying if dot sizes are proportional to cell counts
#' @param show.all_labels a logical specifying if all cluster label must be show or just significant cluster
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggrepel grid
#' 
#' @export
correlatedClustersViewer <- function(CC,
                                     show.cluster.sizes = TRUE,
                                     show.all_labels = FALSE,
                                     show.on_device = TRUE) {
    
    CC@result <- cbind (CC@result, cluster.size = CC@cluster.size)
    
    data.text <- CC@result
    if (!show.all_labels){
        data.text <- subset(CC@result, CC@result$significant)
    }
    title.details <- ifelse(CC@use.percentages, "using % of cells", "using # of cells")
    plot <- ggplot2::ggplot(data = CC@result) +
            ggplot2::ggtitle(paste0("Correlation of clusters kinetics ", title.details, " (", format(sum(CC@cluster.size), big.mark=" "), " cells)", sep = "")) +
            ggplot2::geom_hline(yintercept = -log10(CC@th.pvalue),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1) +  
            ggplot2::geom_vline(xintercept = c(-CC@th.correlation, CC@th.correlation),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1)
    if(show.cluster.sizes){
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significant", size = "cluster.size"),
                                           shape = 21,
                                           colour = "black",
                                           stroke = 1) +
                       ggplot2::scale_size(name = "number.of.cells")
    }else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significant"),
                                           shape = 21,
                                           colour = "black",
                                           stroke = 1)
    }

    x.breaks <- c( - CC@th.correlation, CC@th.correlation, seq(-1, 1, by = 0.1))

    y.max    <- ceiling(max(-log10(CC@th.pvalue), -log10(CC@result$pvalue)))
    y.breaks <- c(seq(0, y.max, by = 1), round(-log10(CC@th.pvalue), 2))
    
    plot <- plot +  ggrepel::geom_text_repel(data          = data.text, 
                                             ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", label = "cluster"),
                                             size          = 3,
                                             box.padding   = grid::unit(0.35, "lines"),
                                             point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey", "red")) +
                    ggplot2::scale_x_continuous(minor_breaks = NULL, limits = c(-1, 1), breaks = x.breaks) +
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab(paste(CC@method, "coeficient of correlation")) +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::theme_bw()

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)

}

#' @title classificationViewer
#'
#' @description 
#' Generate a graph representation of classified clusters
#' 
#' @details 
#' Clusters of the same classe are shown using a circular graph. Circular graphs are sorted by the number of cluster in each class. 
#'
#' @param CCR an object of class 'CCR' (object returned by the 'classifyclusteringResults()' function)
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggnetwork network grDevices gridExtra
#' @importFrom network %s% set.vertex.attribute add.edges get.edge.attribute add.vertices %c% list.vertex.attributes set.edge.attribute is.directed set.vertex.attribute get.vertex.attribute delete.edges is.bipartite get.edges list.edge.attributes delete.vertices
#' @export
#' 
classificationViewer <- function(CCR,
                                 show.on_device = TRUE) {

    classes        <- CCR@classes
    classes        <- na.omit(classes)
    sorted.classes <- names(sort(table(classes$class), decreasing = TRUE))

    colours <- grDevices::rainbow(n = length(sorted.classes), alpha = 0.5)

    plots   <- list()

    for (i in sorted.classes){

        same.class <- classes[classes$class == i,]
        index <- length(plots) + 1
        if (nrow(same.class) >= 2) {

            x <- c()
            y <- c()

            previous <- NA
            for (j in 1:nrow(same.class)) {

                x <- c(x, j)
                y <- c(y, previous)
                previous <- j

            }
            y[1] <- previous

            #            temp <- base::tempfile()
            #            sink(temp) 
            graph <- network::network.initialize(nrow(same.class), directed = FALSE)
            network::set.vertex.attribute(x = graph, attrname = "cluster", value = as.character(same.class$cluster))
            network::add.edges(x = graph, tail = x, head = y)

            graph <- ggnetwork::ggnetwork(graph, layout = "circle")
            #            sink()
            #            
            #            if (file.info(temp)$size > 0){
            #                cat(file = temp)
            #            }
            #            file.remove(temp) 

            
            plots[[index]] <- ggplot2::ggplot(data = graph, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
                              ggnetwork::geom_edges(linetype = "twodash", color = "grey90", size = 1, curvature = 0.1) +
                              ggnetwork::geom_nodes(size = 6, fill = "grey", color = colours[as.numeric(i)], shape = 21, stroke = 3) +
                              ggnetwork::geom_nodetext(ggplot2::aes_string(label = "cluster"), size = 2) +
                              ggnetwork::theme_blank()

        } else {

            plots[[index]] <- ggplot2::ggplot(data = data.frame(x = 0, y = 0), ggplot2::aes_string(x = "x", y = "y")) +
                              ggplot2::geom_point(size = 6, fill = "grey", color = colours[as.numeric(i)], shape = 21, stroke = 3) +
                              ggplot2::geom_text(x = 0, y = 0, label = as.character(same.class$cluster), size = 2) +
                              ggnetwork::theme_blank()
            
        }                       
        
    }

    if(length(plots) == 0){
        plots <- list(grid::rectGrob())
    }

    plot <- gridExtra::arrangeGrob(grobs = plots, top = paste0("ClassificationViewer - (based on ", CCR@type, " using ", CCR@method, " method)"))

    if (show.on_device) {
        grid::grid.newpage()
        grid::grid.draw(plot)
    }

    invisible(plot)

}

#' @title Graphical representation for some SPADEVizR objects
#'
#' @description 
#' This function generates a graphical representation for 'AC', 'DEC', 'CC', 'PhenoProfiles' and 'EnrichmentProfiles' objects.
#'
#' @param x a 'AC', 'DEC', 'CC', 'PhenoProfiles' or 'EnrichmentProfiles' object
#' @param y a supplementary parameter transmited respectively to 'abundantClustersViewer()', 'volcanoViewer()' or 'correlatedClustersViewer()' functions
#' @param ... some supplementaries parameters transmited respectively to \code{\link[SPADEVizR]{abundantClustersViewer}}, \code{\link[SPADEVizR]{volcanoViewer}} or \code{\link[SPADEVizR]{correlatedClustersViewer}} functions
#' 
#' @return a 'ggplot' object
#'  
#' @name plot
#' @rdname plot-methods
setGeneric("plot", function(x, y=NULL, ...){ standardGeneric("plot") })


#' @rdname plot-methods
#' @export
setMethod("plot", c("DEC", "missing"),
        function(x, ...){
            return(volcanoViewer(x, ...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("AC", "missing"),
        function(x, y, ...){
            return(abundantClustersViewer(x, ...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("CC", "missing"),
        function(x, y, ...){
            return(correlatedClustersViewer(x, ...))
        }
)


#' @rdname plot-methods
#' @export
setMethod("plot", c("CCR", "missing"),
        function(x, ...){
            return(classificationViewer(x, ...))
        }
)