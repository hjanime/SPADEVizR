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
#' 
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggrepel
#' 
#' @export
abundantClustersViewer <- function(AC,
                                   show.cluster_sizes = TRUE,
                                   show.all_labels    = FALSE){

    AC@result <- cbind (AC@result, cluster.size = AC@cluster.size)
    
    data.text <- AC@result
    if (!show.all_labels){
        data.text <- subset(AC@result, AC@result$significance)
    }

    plot <-  ggplot2::ggplot(data = AC@result) +
             ggplot2::ggtitle(paste0("Cells abundance of clusters(", format(sum(AC@cluster.size), big.mark=" "), " cells)", sep = "")) +
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
    if (show.cluster_sizes > 0){
       plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significance", size = "cluster.size"), shape = 21, colour = "black", stroke = 1)
    }
    else {
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significance"), shape = 21, colour = "black", stroke = 1)
    }                  
    
    x.max <- ceiling(max(-log10(AC@th.pvalue),-log10(AC@result$pvalue)))
    x.breaks <- seq(0, x.max, by = 1)
    
    y.max <- ceiling(max(AC@th.mean,AC@result$mean))
    y.breaks <- seq(0, y.max, by = 1)

    plot <- plot +  ggrepel::geom_text_repel(data = data.text, 
                                             ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", label = "cluster"),
                                             size          = 3,
                                             box.padding   = grid::unit(0.35, "lines"),
                                             point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey","red")) +
                    ggplot2::scale_x_continuous(limits = c(0, x.max), minor_breaks = NULL, breaks = x.breaks) + 
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab("-log10(p-value)") +
                    ggplot2::ylab(ifelse(AC@use.percentages,"mean (% of cells)","mean (# of cells)")) +                  
                    ggplot2::theme_bw()
    
    return(plot)
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
#'
#' @return a 'ggplot' object
#'  
#' @import ggplot2
#' 
#' @export
volcanoViewer <- function(DEC                = NULL,
                          fc.log2            = TRUE,
                          show.cluster.sizes = TRUE,
                          show.all_labels    = FALSE){
    
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
        data.text <- subset(DEC@result, DEC@result$significance)
    }
    min           <- floor(min(DEC@result$fold.change))
    max           <- ceiling(max(DEC@result$fold.change))
    max           <- max(max,abs(min))
    x.breaks      <- c(round(c(-th.fc, th.fc),2), seq(-max, max, by = 1))
    
    title.details <- ifelse(DEC@use.percentages,"using % of cells","(using # of cells")
        
    plot <- ggplot2::ggplot(data = DEC@result, ggplot2::aes_string(x = "fold.change", y = "-log10(pvalue)")) +
            ggplot2::ggtitle(paste0("Volcano plot showing differentially enriched clusters ", title.details, " (", format(sum(DEC@cluster.size), big.mark=" "), " cells)", sep = "")) +
            ggplot2::geom_vline(xintercept = c(th.fc,-th.fc),
                                       linetype   = "dashed",
                                       alpha      = 0.3,
                                       color      = "red",
                                       size       = 1) +  
            ggplot2::geom_hline(yintercept = -log10(DEC@th.pvalue),
                                       linetype   = "dashed",
                                       alpha      = 0.3,
                                       color      = "red",
                                       size       = 1)
    if (show.cluster.sizes > 0){
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significance", size = "cluster.size"), shape = 21, colour = "black", stroke = 1)
    }
    else {
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significance"), shape = 21, colour = "black", stroke = 1)
    }
                 
    plot <- plot + ggrepel::geom_text_repel(data     = data.text,
                                            ggplot2::aes_string(label = "cluster"),
                                            size          = 3,
                                            box.padding   = grid::unit(0.35, "lines"),
                                            point.padding = grid::unit(0.3, "lines")) +
                   ggplot2::scale_fill_manual(values = c("grey","red")) +
                   ggplot2::scale_x_continuous(limits = c(-max, max) ,minor_breaks = NULL, breaks = x.breaks) +
                   ggplot2::scale_y_continuous(minor_breaks = NULL, breaks = round(-log10(c(DEC@th.pvalue,1,0.1,0.01,0.001)),2)) +
                   ggplot2::xlab(paste0(ifelse(fc.log2,"log2(fold.change)","fold.change"),"\ncond2 < enriched > cond1")) +
                   ggplot2::ylab("-log10(p-value)") +
                   ggplot2::theme_bw()
    
    return(plot)
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
#' 
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggrepel
#' 
#' @export
correlatedClustersViewer <- function(CC,
                                     show.cluster.sizes = TRUE,
                                     show.all_labels    = FALSE){
    
    CC@result <- cbind (CC@result, cluster.size = CC@cluster.size)
    
    data.text <- CC@result
    if (!show.all_labels){
        data.text <- subset(CC@result, CC@result$significance)
    }
    title.details <- ifelse(CC@use.percentages,"using % of cells","using # of cells")
    plot <- ggplot2::ggplot(data = CC@result) +
            ggplot2::ggtitle(paste0("Correlation of clusters kinetics ", title.details," (", format(sum(CC@cluster.size), big.mark=" "), " cells)", sep = "")) +
            ggplot2::geom_hline(yintercept = -log10(CC@th.pvalue),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1) +  
            ggplot2::geom_vline(xintercept = c(-CC@th.correlation,CC@th.correlation),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1)
    if(show.cluster.sizes > 0){
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significance", size = "cluster.size"), shape = 21, colour = "black", stroke = 1)
    }
    else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significance"), shape = 21, colour = "black", stroke = 1)
    }
    
    y.max <- ceiling(max(-log10(CC@th.pvalue),-log10(CC@result$pvalue)))
    y.breaks <- seq(0, y.max, by = 1)
    
    plot <- plot +  ggrepel::geom_text_repel(data = data.text, 
                    ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", label = "cluster"),
                                        size          = 3,
                                        box.padding   = grid::unit(0.35, "lines"),
                                        point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey","red")) +
                    ggplot2::scale_x_continuous(minor_breaks = NULL, limits = c(-1,1),
                                                breaks = c(-CC@th.correlation,CC@th.correlation,seq(-1, 1, by = 0.1))) +
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab(paste(CC@method,"coeficient of correlation")) +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::theme_bw()
    return(plot)
}



#' @title profilesViewer
#'
#' @description 
#' Generate a graph representation of PhenoProfiles classes
#' 
#' @details 
#' xxx
#'
#' @param profile.object a PhenoProfiles object or an EnrichmentProfiles object
#' 
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggnetwork network
#' 
#' @export
#' 
profilesViewer <- function (profile.object){

    classes <- profile.object@classes    
    classes <- na.omit(classes)
    all.sorted.classes <- names(sort(table(classes$class), decreasing = TRUE))

    plots <- list()
    
    for (i in all.sorted.classes){

        same.class <- classes[classes$class == i,]
        
        if (nrow(same.class) >= 2){
            
            x      <- c()
            y      <- c()
            
            previous <- NA
            for (j in 1:nrow(same.class)){
                
                x        <- c(x,j)
                y        <- c(y,previous)
                previous <- j
                
            }
            y[1] <- previous 
            
            graph <- network::network.initialize(nrow(same.class), directed = FALSE)
            
            network::set.vertex.attribute(x = graph, attrname = "cluster", value = as.character(same.class$cluster))

            graph <- network::add.edges(graph,x,y)
            graph <- ggnetwork::ggnetwork(graph, layout = "circle")

            index <- length(plots) + 1
            plots[[index]] <- ggplot2::ggplot(data = graph, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
                              ggnetwork::geom_edges(linetype = "twodash", color = "grey90", size = 1, curvature = 0.1) +
                              ggnetwork::geom_nodes(size = 6, fill = "grey", color = "black", shape = 21, stroke = 3) +
                              ggnetwork::geom_nodetext(ggplot2::aes_string(label = "cluster"), size = 2) +
                              ggnetwork::theme_blank()
          }
    }

    ret <- gridExtra::grid.arrange(grobs = plots, top = paste0("Profiles Viewer (",names(profile.object)," using ",profile.object@method," method)"))
    
    return(ret)

}

#' @title Graphical representation for some SPADEVizR objects
#'
#' @description 
#' This function generates a graphical representation for 'AC', 'DEC', 'CC', 'PhenoProfiles' and 'EnrichmentProfiles' objects.
#'
#' @param x a 'AC', 'DEC', 'CC', 'PhenoProfiles' or 'EnrichmentProfiles' object
#' 
#' @return a 'ggplot' object
#'  
#' @name plot
#' @rdname plot-methods
setGeneric("plot", function(x,y=NULL,...){ standardGeneric("plot") })


#' @rdname plot-methods
#' @export
setMethod("plot",c("DEC","missing"),
        function(x,...){
            return(volcanoViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("AC","missing"),
        function(x,y,...){
            return(abundantClustersViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("CC","missing"),
        function(x,y,...){
            return(correlatedClustersViewer(x,...))
        }
)


#' @rdname plot-methods
#' @export
setMethod("plot",c("PhenoProfiles","missing"),
        function(x,...){
            return(profilesViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("EnrichmentProfiles","missing"),
        function(x,...){
            return(profilesViewer(x,...))
        }
)