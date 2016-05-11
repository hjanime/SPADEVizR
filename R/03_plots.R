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
#' @param show.cluster.sizes a logical specifying if dot sizes are proportional to cell counts
#' 
#' @return a 'ggplot' object
#' 
#' @import ggplot2 ggnetwork grDevices network 
#' 
#' @export
#' 
profilesViewer <- function (profile.object,
                            show.cluster.sizes = TRUE){

    classes <- profile.object@classes    
    classes <- na.omit(classes)
    all.sorted.classes <- names(sort(table(classes$class), decreasing = TRUE))

    plots <- list()
    
    for (i in all.sorted.classes){
        print(i)
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

            network::set.vertex.attribute(x = graph, attrname = "cluster", value = same.class$cluster)
            
            profile.object@cluster.size[same.class$cluster]
            
            if (show.cluster.sizes){
                print(profile.object@cluster.size[same.class$cluster])
                network::set.vertex.attribute(x = graph, attrname = "size", value = profile.object@cluster.size[same.class$cluster])
            }

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

#' @title Visualization of cluster phenotypes
#' 
#' @description 
#' Generates a parallel coordinate plot representation showing for each cluster the marker median expressions.
#' 
#' @details 
#' If 'Results' is an object of class 'SPADEResults' object then the quantile ranges will be displayed using a ribbon.
#' 
#' The 'samples' parameter required a named vector providing the correspondence between a sample name (in rowname) and the logical value TRUE to show this sample or FALSE otherwise.
#' 
#' The 'show.mean' parameter allows to visualize three kinds of information:
#' \itemize{
#' \item "none" value will show marker median expressions for each selected samples;
#' \item "only" value will show only the mean of median maker expressions for all selected samples (displayed as black dashed line);
#' \item "both" value will show marker median expressions for each selected samples together with the mean of median maker expressions for all selected samples.
#' }
#' 
#' @param Results a SPADEResults or Result object
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param samples a named vector specifying the samples to used
#' @param markers a character vector specifying the markers to be displayed 
#' @param show.mean a character specifying if marker means expression should be displayed, possible value are among : "none", "only" or "both"
#' 
#' @return a list of ggplot objects
#'
#' @import reshape2 ggplot2
#' 
#' @export
clusterViewer <- function(Results,
                          clusters      = NULL,
                          samples       = NULL,
                          markers       = NULL,
                          show.mean     = "both"){

    if(show.mean != "none" && show.mean != "both" && show.mean != "only"){
        stop("Error : show.mean must be one of those : 'none' 'both' 'only' ")
    }

    if(is.null(samples)){
        data        <- Results@marker.expressions
        cells.count <- Results@cells.count
    }else{
        data        <- subset(Results@marker.expressions, sample %in% names(samples[ samples == TRUE]), drop = FALSE)
        cells.count <- Results@cells.count[,c(names(samples[ samples == TRUE] )), drop = FALSE]
    }

    data <- na.omit(data)
    
    if(!is.null(clusters)){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters.select <- data[,"cluster"] %in% clusters
        data            <- data[clusters.select,]
    }else{
        clusters <- unique(data[,"cluster"])
    }
    
    if(!is.null(markers)){
        data.keys <- data[,c("sample","cluster")]
        data <- cbind(data.keys,data[,markers])
    }
    
    if (names(Results) == "SPADEResults"){
        markers <- colnames(data[, grep ("cluster|sample",colnames(data), invert = TRUE)])
        clustering.markers <- is.element(markers,Results@marker.names[Results@marker.clustering])
        bold.markers <- ifelse(clustering.markers,"bold","plain")
    }
    
    data <- reshape2::melt(data, id = c("sample","cluster"),stringsAsFactors=FALSE)
    colnames(data) <- c("sample","cluster","marker","value")
    
    if (names(Results) == "SPADEResults"){
        for(i in 1:nrow(data)){
            data[i,"lower.bound"] <- Results@quantiles[1,as.character(data[i,"marker"])]
            data[i,"upper.bound"] <- Results@quantiles[2,as.character(data[i,"marker"])]
        }
    }
    
    plots <- list()
    
    for(current.cluster in clusters){
        
        data.temp  <- data[data["cluster"] == current.cluster,]
        if (names(Results) == "SPADEResults"){
            max.value <- max(c(data.temp$value,data.temp$upper.bound))
            min.value <- min(c(data.temp$value,data.temp$lower.bound))
        }else{
            max.value <- max(c(data.temp$value))
            min.value <- min(c(data.temp$value))
        }
        
        max.value <- max.value + 0.1 * max.value
        min.value <- min.value + 0.1 * min.value
         
        i <- length(plots) + 1
        
        cells.number <- sum(cells.count[current.cluster,])
                
        plots[[i]] <- ggplot2::ggplot(data = data.temp) +
                      ggplot2::ggtitle(paste("cluster ",current.cluster," - Cluster Viewer (", format(cells.number, big.mark=" "), " cells)" ,sep = ""))                      
        
        if(show.mean == "both" || show.mean == "none"){
            plots[[i]] <- plots[[i]] + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "value", group = "sample", color = "sample"), size = 1)                               
        }     
        if(show.mean == "only" || show.mean == "both"){
            wide <- reshape2::dcast(data.temp,sample ~ marker)
            means <- apply(wide[,2:ncol(wide)],2,mean)
            df.means <- data.frame(marker = names(means), means = means)
            df.means$show.mean <- rep("show.mean",nrow(df.means))
            
            plots[[i]] <- plots[[i]] + ggplot2::geom_line(data = df.means, ggplot2::aes_string(x = "marker", y = "means", group="show.mean"), linetype = "dashed", size = 1)
        }
        if(show.mean == "only"){
            plots[[i]] <- plots[[i]] + ggplot2::theme(legend.position  = "none")
        }
        
        if(names(Results) == "SPADEResults"){
            plots[[i]] <- plots[[i]] + ggplot2::geom_ribbon(ggplot2::aes_string(x = "as.numeric(marker)", ymin = "lower.bound", ymax = "upper.bound"),alpha = 0.1, fill = "grey20")
        }
        
        plots[[i]] <- plots[[i]] + ggplot2::scale_x_discrete(limits = markers) +
                                   ggplot2::scale_y_continuous(limits = c(min.value,max.value),breaks = round(seq(0, max.value , by = 1),0)) +
                                   ggplot2::theme_bw()         
  
        if (names(Results) == "SPADEResults"){
            plots[[i]] <- plots[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0, vjust = 1, face = bold.markers),
                                                      legend.text = ggplot2::element_text(size = 6))                 
        }else{
            plots[[i]] <- plots[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0, vjust = 1),
                                                      legend.text = ggplot2::element_text(size = 6))
        }
                                   
        
    }
    
    return(plots)    
}

#' @title Visualization of all clusters phenotypes as an heatmap
#' 
#' @description 
#' Generates an heatmap representation showing for all clusters the marker median expressions.
#' 
#' @details 
#' For each marker, median expressions are discretized in severals categories corresponding to the heat intensities. 
#' This number of categories is provided using 'num' parameter.
#' 
#' By default the 'pheno.table' parameter is null and will be automatically compute using the 'computePhenoTable' function.
#' 
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param num a numeric value specifying the number of markers expression categories to use
#'
#' @return a list of ggplot objects
#'
#' @import reshape2
#' 
#' @export
phenoViewer <- function(SPADEResults,
                        num         = 5){
    #print(SPADEResults)
    if (names(SPADEResults) == "Results"){
       stop("Error : phenoViewer required a SPADEResults object")
    }                
                    
    pheno.table <- computePhenoTable(SPADEResults, num)
        
    pheno.table <- reshape2::dcast(pheno.table, cluster ~ marker)
    pheno.table <- pheno.table[,2:ncol(pheno.table)]
    pheno.table <- t(pheno.table)
    pheno.table <- as.matrix(pheno.table)

    plot.elements <- ggheatmap(pheno.table, num = num, clustering.markers = SPADEResults@marker.names[SPADEResults@marker.clustering])

    heatmap <- ggheatmap.plot(plot.elements)
    
    return(heatmap)
    
}

#' @title Visualization of combined SPADE trees
#' 
#' @description 
#' Generates a tree representation showing combined SPADE trees. 
#' 
#' 
#' @details 
#' The size of tree nodes are related to the number of cells in each cluster. 
#' If the 'stat.object' parameter is provided node outlines are colored according to clusters signifiance.
#' If the 'marker' parameter is provided, the nodes are colored according to mean expression for the selected marker using selected samples.
#' 
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param samples a named vector providing the correspondence between samples name (in rowname) and the logical value TRUE to use these samples (all samples by default)
#' @param marker a character specifying the marker name to display
#' @param stat.object an AC, DEC or CC object to highligth identified significant clusters in the SPADE tree
#'
#' @return a list of ggplot objects
#'
#' @import reshape2 ggplot2
#' 
#' @export
treeViewer <- function(SPADEResults,
                       samples       = NULL,
                       stat.object   = NULL,
                       marker        = NULL){
    
    if (names(SPADEResults) == "Results"){
        stop("Error : treeViewer required a SPADEResults object")
    }
                   
    data   <- SPADEResults@cells.count
    
    if(!is.null(samples)){ 
        data   <- data[names(samples[ samples == TRUE ]), drop = FALSE]
    }else{
        data   <- data

    }
        
    vertex.size <- apply(data,1,sum)

    pos.vertex  <- data.frame(id   = as.character(1:nrow(SPADEResults@graph.layout)),
                              x    = SPADEResults@graph.layout[,1],
                              y    = SPADEResults@graph.layout[,2],
                              size = vertex.size)
                      
    if (!is.null(marker)){
        if(!is.null(samples)){ 
              expr   <- subset(SPADEResults@marker.expressions, sample %in% names(samples[ samples == TRUE]), drop = FALSE)
        }else{
              expr   <- SPADEResults@marker.expressions
        }
                          
       expr <- expr[,c("cluster", "sample", marker)]
       expr <- reshape2::dcast(expr, cluster ~ sample, value.var = marker)
       rownames(expr) <- expr$cluster
       expr <- expr[, -1]
       mean.expr <- apply(expr, 1, mean, na.rm = TRUE)
       pos.vertex <- cbind(pos.vertex, marker = mean.expr)
       colnames(pos.vertex)[5] <- marker
   }   
   print(head(pos.vertex))                   
      
    edges    <- igraph::get.edgelist(SPADEResults@graph,names = FALSE)
    pos.edge <- data.frame(x    = SPADEResults@graph.layout[edges[,1],1],
                           xend = SPADEResults@graph.layout[edges[,2],1],
                           y    = SPADEResults@graph.layout[edges[,1],2],
                           yend = SPADEResults@graph.layout[edges[,2],2])

    cells.number <- sum(colSums(data))
    
    plot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste("Tree Viewer (", format( cells.number, big.mark = " "), " cells)", sep = "")) +
            ggplot2::geom_segment(data = pos.edge, ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"))
    
    if (!is.null(stat.object)){
        stat.object.name <- names(stat.object)
        pos.vertex[, stat.object.name] <- stat.object@result$significance
        if(!is.null(marker)){
            plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker, colour = stat.object.name), stroke = 3, shape = 21)  
        }else{
            plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size", colour = stat.object.name), fill = "grey80", stroke = 3, shape = 21)  
        }
    }else{
        if(!is.null(marker)){
            plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker), stroke = 3, shape = 21)
        }else{
            plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size"), fill = "grey80", stroke = 3, shape = 21)
        }
    }
    plot <- plot + ggplot2::scale_fill_gradient(low = "#ECE822", high = "#EE302D") +#low = yellow, high = red
                   ggplot2::scale_color_manual(values = c("black","blue")) +
                   ggplot2::scale_size_area(max_size = 15) +
                   ggrepel::geom_label_repel(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", label = "id"),size = 4, color = "black", box.padding = grid::unit(0.1, "lines"), point.padding = grid::unit(0.1, "lines")) +
                   ggplot2::coord_fixed() +
                   ggplot2::theme(panel.background = ggplot2::element_blank(),
                                  panel.border     = ggplot2::element_blank(),
                                  axis.text.x      = ggplot2::element_blank(),
                                  axis.text.y      = ggplot2::element_blank(),
                                  axis.title.x     = ggplot2::element_blank(),
                                  axis.title.y     = ggplot2::element_blank(),
                                  panel.grid.minor = ggplot2::element_blank(),
                                  panel.grid.major = ggplot2::element_blank(),
                                  axis.ticks       = ggplot2::element_blank(),
                                  legend.position  = "right")
    return(plot)
}

#' @title Visualization of cluster sizes
#'
#' @description 
#' Generate a two dimensional vizualisation showing the number of cells (sum of selected samples) of each cluster.
#' 
#' @details 
#' xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param samples a named vector providing the correspondence between samples name (in row-names) and the logical value TRUE to use these samples (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param min.cells a numeric specifying the minimum number of cell (sum of all selected samples) to display a cluster
#' @param sort a logical specifying if clusters will be to be sorted (descending) based on the sum of all selected samples for each cluster
#' @param show.samples a logical specifying if the number of cells for all selected samples will be displayed
#' 
#' @return a 'ggplot' object
#' 
#' @import ggplot2 reshape2
#' 
#' @export
countViewer <- function(Results,
        samples      = NULL,
        clusters     = NULL,
        min.cells    = 0,
        sort         = TRUE,
        show.samples = TRUE){
    
    if(is.null(samples)){ 
        data <- Results@cells.count
    }else{
        data  <- Results@cells.count[, names(samples[ samples == TRUE]), drop = FALSE]
    }
    
    #print(head(data))
    
    if(!is.null(clusters)){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        data <- data[clusters, , drop = FALSE]
    }
    
    cells.number <- sum(colSums(data))
    
    data  <- cbind(data, "sum.of.samples" = apply(data, 1, sum))
    data  <- data[data$sum.of.samples > min.cells, ]
    
    data <- cbind("cluster" = rownames(data), data)
    
    if (sort){
        data <- transform(data, "cluster" = reorder(cluster, -sum.of.samples))
    }else{
        data$cluster <- as.factor(data$cluster)
    }
    
    data.melted <- reshape2::melt(data, id = c("cluster"))  
    colnames(data.melted) <- c("cluster", "sample", "value")
    
    data.melted$total <- ifelse(data.melted[,"sample"] == "sum.of.samples","sum of selected samples","")
    
    plot <- ggplot2::ggplot(data = data.melted) +
            ggplot2::ggtitle(paste("Count Viewer (", format(cells.number, big.mark=" "), " cells)", sep = "")) +
            ggplot2::geom_point(data = subset(data.melted, sample == "sum.of.samples"),
                    ggplot2::aes_string(x = "cluster", y = "value", size = "value", shape = "total"),
                    fill = "grey40") +
            ggplot2::scale_shape_manual(values = 21)        
    if(show.samples){                
        plot <- plot + ggplot2::geom_jitter(data   = subset(data.melted, sample != "sum.of.samples"),
                height = 0,
                width  = 0.5,
                ggplot2::aes_string(x = "cluster", y = "value", size = "value", fill = "sample"),
                shape  = 21,
                alpha  = 0.4)
    }
    plot <- plot + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,1.1*max(data.melted$value))) +
            ggplot2::ylab("# of cells") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0, vjust = 1))
    
    return(plot)
}

#' @title Visualization of cluster enrichment profiles kinetics
#' 
#' @description 
#' Generates a kinetics plot representation showing for each cluster its enrichment profiles at each time-point of each individual.
#' 
#' @details 
#' xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' 
#' @return a 'ggplot' object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
kineticsViewer <- function(Results,
                           assignments,
                           clusters        = NULL,
                           use.percentages = TRUE){
    
    if(missing(assignments) || is.null(assignments)){
        stop("Error : the 'assignments' parameter is required")   
    }
                       
    data <- Results@cells.count
    cells.count <- data[,rownames(assignments), drop = FALSE]
    
    if(use.percentages){
        data.percent <- prop.table(as.matrix(data), 2) * 100
        print(data.percent)
        data         <- data.frame(data.percent)#row.names = rownames(data),
        legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
    
    if(is.null(clusters)){
        clusters <- rownames(data)
        message("All clusters will be compute")
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data     <- data[clusters,]
        message("These clusters will be compute :")
        message(paste0(clusters,collapse="\t"))
    }else{
        stop("Error : Unknown cluster ")
    }
        
    data <- cbind(cluster = rownames(data), data)
    data.melted <- reshape2::melt(data, id = "cluster")
    
    colnames(data.melted) <- c ("cluster", "sample", "value")

    data.melted$individuals <- assignments[data.melted$sample,'individuals']
    data.melted$timepoints  <- assignments[data.melted$sample,'timepoints']

    plots <- list()
    for(current.cluster in clusters){
        #print(clusters)
        data.temp  <- data.melted[data.melted$cluster == current.cluster,]
        
        max.value  <- max(data.temp$value)
        max.value  <- max.value + max.value*0.1 + 1
        
        i <- length(plots) + 1
        #print(current.cluster)
        #print(cells.count)
        cells.number <- sum(cells.count[current.cluster,])

        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "as.factor(timepoints)", y = "value", group = "individuals", color = "individuals")) +
                      ggplot2::ggtitle(paste("cluster ",current.cluster," - Kinetics Viewer (", format(cells.number, big.mark = " "), " cells)",sep = "")) +
                      ggplot2::geom_line() +
                      #ggalt:::geom_xspline(size = 1, na.rm = TRUE) + #add smooth curve with spline function
                      ggplot2::geom_point(na.rm = TRUE) +
                      ggplot2::scale_x_discrete(expand = c(0,0.05))
       if(use.percentages){    
           plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value),breaks = round(seq(0, max.value)), minor_breaks = NULL)
       }else{
           plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value))
       }

       plots[[i]] <- plots[[i]] + ggplot2::ylab(legendy) +
                                  ggplot2::xlab("timepoints") +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0),
                                                 legend.text = ggplot2::element_text(size = 6))            
    }
    return(plots)
}

#' @title Visualization of cluster enrichment profiles conditions
#'
#' @description
#' Generate a boxplot representation displaying the cell abundance for each cluster. Clusters are gathered by given biological conditions.
#' 
#' @details
#' Cells clusters are colored based on theirs associated biological samples. 
#' 
#' @param Results a SPADEResults or Results object
#' @param conditions conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample : NA to exclude a sample from tests
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param use.percentages a logical specifying if the visualization must be performed on percentage
#' @param show.legend a logical specifying if the legend must be displayed
#' @param show.violin a logical specifying if the count distribution must be displayed
#'
#' @return a 'ggplot' object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
boxplotViewer <- function(Results,
                          conditions,
                          clusters        = NULL,
                          use.percentages = TRUE,
                          show.legend     = FALSE,
                          show.violin     = TRUE){
                      
    data <- Results@cells.count
    
    cells.count <- data[,names(conditions[!is.na(conditions)]), drop = FALSE]
    
    if(use.percentages){
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent)
        legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
    
    if(is.null(clusters)){
        clusters <- rownames(data)
        message("All clusters will be compute")
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data     <- data[clusters,]
        message("These clusters will be compute:")
        message(paste0(clusters,collapse="\t"))
    }else{
        stop("Error : Unknown cluster ")
    }
          
    data <- cbind(cluster = rownames(data), data)
    data.melted <- reshape2::melt(data, id = "cluster")
                      
    colnames(data.melted) <- c ("cluster", "sample", "value")
                    
    data.melted$cond <- conditions[data.melted$sample]
    data.melted <- data.melted[!is.na(data.melted$cond),]
    #print(data.melted)
    plots <- list()
    
    for(current.cluster in clusters){
        data.temp  <- data.melted[data.melted$cluster == current.cluster,]
        #print(data.temp$value)
        max.value <- max(data.temp$value)
        max.value <- max.value + 0.1*max.value + 1
                
        data.temp$cond <- as.factor(data.temp$cond)
        
        i <- length(plots) + 1
        
        cells.number <- sum(cells.count[current.cluster,])

        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "cond", y = "value")) +
                ggplot2::ggtitle(paste("cluster ",current.cluster," - Boxplot Viewer (", format(cells.number, big.mark=" "), " cells)",sep = "")) +
                ggplot2::geom_boxplot() +
                ggplot2::geom_jitter(ggplot2::aes_string(color = "sample"), width = 0.2, show.legend = show.legend)
        
        if(show.violin){
            plots[[i]] <- plots[[i]] + ggplot2::geom_violin(alpha = 0.05, fill = "red", colour = "red")
        }
        
        if(use.percentages){    
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value),breaks = round(seq(0, max.value)), minor_breaks = NULL)
        }else{
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value))
        }
        
        plots[[i]] <- plots[[i]] + ggplot2::ylab(legendy) +
                                   ggplot2::xlab("biological conditions") +
                                   ggplot2::theme_bw() +
                                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0),
                                                  legend.text = ggplot2::element_text(size = 6))
                                   
    }  
    return(plots)
}

#' @title Visualization of marker co-expressions
#'
#' @description 
#' Generate a distogram representation showing the marker co-expressions.
#' 
#' @details 
#' A pearson correlation matrix is calculated between selected markers using selected clusters and samples.
#' High positive correlated markers are shown by a green tile at their perpendicular intersection. 
#' In the same way, absence of correlation are shown by black tiles and negative correlation by red tiles.
#'
#' @param Results a SPADEResults or Results object
#' @param clusters a character vector containing the clusters names to be use (by default all clusters will be used)
#' @param samples a character vector specifying the samples to use 
#' @param markers a character vector specifying the markers to be displayed 
#' 
#' @return a list of 'ggplot' objects
#' 
#' @import reshape2 ggplot2
#' 
#' @export
distogramViewer <- function(Results,
                            clusters = NULL,
                            samples  = NULL,
                            markers  = NULL){

    if(is.null(samples)){
        data        <- Results@marker.expressions
        cells.count <- Results@cells.count
    }else{
        data        <- subset(Results@marker.expressions, sample %in% names(samples[ samples == TRUE]), drop = FALSE)
        cells.count <- Results@cells.count[, c(names(samples[ samples == TRUE] )), drop = FALSE]
    }
    
    #print(head(data))
    
    if(!is.null(clusters)){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters.select <- data[,"cluster"] %in% clusters
        data            <- data[clusters.select,]
        cells.count     <- cells.count[clusters,]
    }
    
    #print(head(data))
    
    data <- data[ , -c(1, 2) ]
    if(!is.null(markers)){
        data <- cbind(data[,markers])
    }
       
    data <- na.omit(data)# NA values are removed
    
    cormat <- round(cor(data), 2)
    dist <- as.dist(1 - cormat)
    hc <- hclust(dist)
    cormat <-cormat[hc$order, hc$order]
    
    cormat[upper.tri(cormat, diag = TRUE)] <- NA
    
    markers <- colnames(cormat)
    
    dimnames(cormat) <- NULL
    
    melted.cormat <- reshape2::melt(cormat)#TODO maybe use a data.frame
    
    bold.markers <- "plain"
    if (names(Results) == "SPADEResults"){
        clustering.markers <- is.element(markers, Results@marker.names[Results@marker.clustering])
        bold.markers <- ifelse(clustering.markers, "bold", "plain")
    }
    
    plot <- ggplot2::ggplot(data = melted.cormat, ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")) + 
            ggplot2::ggtitle("Distogram of marker phenotypes correlations") +
            ggplot2::geom_tile(color = "white", ) +
            ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black", 
                                          midpoint = 0, limit = c(-1,1), na.value = 'white',
                                          name = "Pearson\nCorrelation") +
            ggplot2::annotate(geom     = "text",
                              x        = 1:length(markers),
                              y        = 1:length(markers),
                              color    = "blue",
                              angle    = -45,
                              size     = 4,
                              label    = markers,
                              hjust    = 1,
                              fontface = bold.markers) +
            ggplot2::coord_fixed() +
            ggplot2::theme(axis.line        = ggplot2::element_blank(),
                           axis.text.x      = ggplot2::element_blank(),
                           axis.text.y      = ggplot2::element_blank(),
                           axis.ticks       = ggplot2::element_blank(),
                           axis.title.x     = ggplot2::element_blank(),
                           axis.title.y     = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(),
                           panel.border     = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background  = ggplot2::element_blank(),
                           legend.justification = c(1, 0),
                           legend.position = c(0.4, 0.7),
                           legend.direction = "horizontal") +
            ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth       = 7,
                            barheight      = 1,
                            title.position = "top",
                            title.hjust    = 0.5))
    
    return(plot)
}

#' @title Visualization of cluster abundance dynamics
#'
#' @description 
#' Generate a streamgraph representation showing the dynamic evolution of the number of cells in clusters across samples.
#' 
#' @details 
#' xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param order a named vector providing the correspondence between a sample name (in rownames) and an integer ordering samples (NA to exclude a sample)
#' @param clusters a character vector containing the clusters names to be vizualised (by default all clusters will be displayed)
#' @param use.relative a logical specifying if the visualization should be performed on relative abundance
#'
#' @return a 'ggplot' object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
streamgraphViewer <- function(Results,
                              order          = NULL,
                              clusters       = NULL,
                              use.relative   = FALSE){
    
    data <- Results@cells.count
       
    if(!is.null(order)){
        data    <- data[,names(order[!is.na(order)]), drop = FALSE]
    }
    
    if(is.null(clusters)){
        clusters <- rownames(data)
        message("All clusters will be compute")
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        data     <- data[clusters,]
        message("These clusters will be compute :")
        message(paste0(clusters,collapse="\t"))
    }else{
        stop("Error : Unknown cluster ")
    }
    
    cells.number <- sum(colSums(data))
    
    if (use.relative){
        data         <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data)
    }

    data <- cbind(cluster = rownames(data), data)
    melted.data           <- reshape2::melt(data, id = "cluster")
    colnames(melted.data) <- c("cluster","sample","value")

    melted.data$cluster <- factor(melted.data$cluster, levels = gtools::mixedsort(unique(melted.data$cluster)))
    melted.data <- melted.data[order(melted.data$sample, melted.data$cluster, decreasing = TRUE),]    
    
    dt.melted.data <- data.table::setDT(melted.data)
    dt.melted.data[, ymax := cumsum(value) - (sum(value)/2), by = sample]
    dt.melted.data[, ymin := ymax - value, by = sample]
    dt.melted.data[, label := format(round(ymax - min(ymin), 2), big.mark=" "), by = sample]
    dt.melted.data[, ybase := min(ymin), by = sample]
    
    print(dt.melted.data)
    
    title <- paste("Streamgraph with ",ifelse(use.relative,"relative","absolute")," abundance (", format(cells.number, big.mark=" "), " cells)", sep = "")
    plot = ggplot2::ggplot(data = dt.melted.data) +
           ggplot2::ggtitle(title) +
           ggplot2::geom_ribbon(ggplot2::aes_string(x = "sample", ymin = "ymin", ymax = "ymax", group = "cluster", fill = "cluster")) +
           ggplot2::geom_point(ggplot2::aes_string(x = "sample", y = "ymax", group = "cluster"), shape = 45) +
           ggplot2::geom_point(ggplot2::aes_string(x = "sample", y = "ybase", group = "cluster"), shape = 45) +
           ggplot2::geom_text(ggplot2::aes_string(x = "sample", y = "ymax", label = "label"), angle = 360, hjust = 1.1, size = 3) +
           ggplot2::geom_text(ggplot2::aes_string(x = "sample", y = "ybase"), label = "0", angle = 360, hjust = 1.1, size = 3) +
           ggplot2::theme_bw() +
           ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                          axis.text.x = ggplot2::element_text(angle = 290, hjust = 0),
                          axis.line        = ggplot2::element_blank(),
                          axis.text.y      = ggplot2::element_blank(),
                          axis.ticks       = ggplot2::element_blank(),
                          axis.title.x     = ggplot2::element_blank(),
                          axis.title.y     = ggplot2::element_blank(),
                          panel.background = ggplot2::element_blank(),
                          panel.border     = ggplot2::element_blank(),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank())


    return(plot)
          
}

#' @title Visualization of SPADE cluster similarities using MDS
#'
#' @description 
#' Generate a Multidimensional Scaling (MDS) representation showing the similarities between SPADE results based on theirs phenotypes.
#' 
#' @details 
#' talk about space xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param space a character specifying the space ("clusters" or "samples", cluster by default)
#' @param dist.method a character string containing the name of the distance measure to use
#' 
#' @return a list of ggplot objects
#' 
#' @import MASS ggplot2 ggrepel grDevices
#' 
#' @export
MDSViewer <- function(Results,
                      use.percentages = TRUE,
                      assignments,
                      clusters        = NULL,
                      space           = "clusters",
                      dist.method     = "euclidean"){
    
    data <- Results@cells.count

    cells.count <- data
    
    if(use.percentages){
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent)
        legendy      <- "% of cells relative to parent"
    }else{
        legendy      <- "# of cells"
    }
    
    if(is.null(clusters)){
        clusters <- rownames(data)
        message("All clusters will be compute")
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data     <- data[clusters,]
        message("These clusters will be compute :")
        message(paste0(clusters,collapse="\t"))
    }else{
        stop("Error : Unknown cluster ")
    }
    
    data <- cbind(cluster = rownames(data), data)               
    
    if(space == "samples"){
        if(missing(assignments) || is.null(assignments)){
            stop("Error : the 'assignments' parameter is required when the 'space' parameter is \"sample\"")   
        }
        data <- t(data[,colnames(data) != "cluster"])
    }else if(space != "clusters"){
        stop("Error in \"space\" parameter")
    }
    
    dist <- dist(data, method = dist.method)
    message("MDS computation")
    
    fit    <- MASS::isoMDS(dist, k = 2)
    stress <- fit$stress
    fit    <- fit$point
    message("done")
    
    x       <- fit[,1]
    y       <- fit[,2]
    
    data_i = data.frame(x = x,y = y)
    
    min.lim <- min(min(x),min(y))*1.1
    max.lim <- max(max(x),max(y))*1.1
    lim     <- max(abs(min.lim),abs(max.lim))
    min.lim <- -lim
    max.lim <- lim
    
    
    if(space == "samples"){
        
        data_i   <- cbind (data_i, assignments)
        samples <- rownames(data)
        
        data_i$individuals <- as.factor(data_i$individuals)
        data_i$timepoints  <- as.factor(data_i$timepoints)
        
        cells.number <- sum(colSums(cells.count[,rownames(assignments)]))
        
        data.table_i <- data.table::data.table(data_i, key = "individuals")
        hulls <- data.table_i[, .SD[grDevices::chull(x, y)], by = "individuals"]

        plot <- ggplot2::ggplot(data = data_i)  +
                ggplot2::ggtitle(paste("MDS with samples (", format(cells.number, big.mark = " "), " cells)", sep = "")) +
                ggplot2::geom_hline(yintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_vline(xintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_polygon(data = hulls, ggplot2::aes_string(x = "x", y = "y", group = "individuals", fill = "individuals"), colour = "black", alpha = 0.3) +
                ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", colour = "individuals", shape = "timepoints"), size = 4) +
                ggplot2::xlim(min.lim, max.lim) +
                ggplot2::ylim(min.lim, max.lim) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                        panel.border     = ggplot2::element_rect(fill = NA),
                        axis.text.x      = ggplot2::element_blank(),
                        axis.text.y      = ggplot2::element_blank(),
                        axis.title.x     = ggplot2::element_blank(),
                        axis.title.y     = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.grid.major = ggplot2::element_blank(),
                        axis.ticks       = ggplot2::element_blank(),
                        legend.position  = "right") +
                ggplot2::annotate(geom  = "text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1, label = paste0("Kruskal Stress : ",round(stress,2)))
        
    }else{
        data_i <- cbind (data_i, cluster = data[,"cluster"])
        data_i$cluster <- as.factor(data_i$cluster)
        
        cells.number <- sum(colSums(cells.count))
        plot <- ggplot2::ggplot(data = data_i)  +
                ggplot2::ggtitle(paste("MDS with clusters (", format(cells.number, big.mark=" "), " cells)", sep = "")) +
                ggplot2::geom_hline(yintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_vline(xintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", color = "cluster"), size = 2) +
                ggrepel::geom_text_repel(ggplot2::aes_string(x = "x", y = "y", label = "cluster", color = "cluster"), size = 5) +
                ggplot2::xlim(min.lim, max.lim) +
                ggplot2::ylim(min.lim, max.lim) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                        panel.border     = ggplot2::element_rect(fill = NA),
                        axis.text.x      = ggplot2::element_blank(),
                        axis.text.y      = ggplot2::element_blank(),
                        axis.title.x     = ggplot2::element_blank(),
                        axis.title.y     = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.grid.major = ggplot2::element_blank(),
                        axis.ticks       = ggplot2::element_blank(),
                        legend.position  = "none") +
                ggplot2::annotate(geom  = "text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1, label = paste0("Kruskal Stress : ",round(stress,2)))
    }
    
    return(plot)
}




#' @title biplotViewer
#'
#' @description 
#' Generates a biplot representation with two markers
#'
#' @details 
#' In such representation, each dot corresponds to a cell profile and dot are ploted in a 2-dimentional space corresponding to the marker expressions. 
#' resample.ratio nb cell not the same than displayed dots 
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param x.marker a character indicating the marker name of the first dimension
#' @param y.marker a character indicating the marker name of the second dimension
#' @param samples xxx
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be used)
#' @param default.min a numeric value indicating the lower bound of the biplot representation 
#' @param sample.merge xxx
#' @param resample.ratio a numric ratio (between 0 and 1) specifying the resample ratio to show less dots 
#' 
#' @return if return.gg is TRUE, the function returns a list of ggplot objects
#'
#' @import grDevices ggplot2
#' 
#' @export
biplotViewer <- function(SPADEResults,
                         x.marker,
                         y.marker,
                         samples        = NULL, 
                         clusters       = NULL,
                         sample.merge   = FALSE,
                         resample.ratio = NULL){# take care about log10(x+1) transformation
    if (names(SPADEResults) == "Results"){
        stop("Error : biplotViewer required a SPADEResults object")
    }  
    if (is.null(SPADEResults@flowset)){
        stop("Error : The 'flowset' slot of the 'SPADEResults' object is not loaded, use the function load.flowSet before using the 'biplotViewer' function")
    }
        
    flowset <- SPADEResults@flowset
    
    x.data <- c()
    y.data <- c()
    facet  <- c()
    flowset.samples <- flowCore::sampleNames(flowset)

    plots <- list()
    
    if(is.null(samples)){
        samples        <- rep(TRUE, length(SPADEResults@sample.names)) 
        names(samples) <- SPADEResults@sample.names
    }
    
    for (sample in flowset.samples){
        if (samples[sample]){
            flowframe <- flowset[[which(sample == flowset.samples),]]
            exprs <- flowframe@exprs
            if (!is.null(clusters)){
                exprs <- subset(exprs, exprs[,"cluster"] %in% clusters)
            }
            x.data <- c(x.data,exprs[,x.marker])
            y.data <- c(y.data,exprs[,y.marker])
            
            if (!sample.merge){
                facet  <- c(facet,rep(sample,length(exprs[,x.marker])))
            }
        }    
    }
    
    if(is.null(samples)){
        cells.count <- SPADEResults@cells.count
    }else{
        cells.count <- SPADEResults@cells.count[,c(names(samples[ samples == TRUE] )), drop = FALSE]
    }
    cells.number.by.sample <- colSums(cells.count)

    if (sample.merge){
       data <- data.frame(x = x.data, y = y.data)
    }else{
       data <- data.frame(x = x.data, y = y.data, facet = paste0(facet, " (",format(cells.number.by.sample[facet], big.mark=" ")," cells)"))
    }
    
    if (!is.null(resample.ratio)){
        if (resample.ratio > 0 && resample.ratio < 1){
            data <- data[sample(nrow(data), round((nrow(data)*resample.ratio))),]
        } else {
            stop("resample.ratio must be > 0 and < 1 or null")
        }
    }
   

    x.max              <- max(data["x"])*1.1
    y.max              <- max(data["y"])*1.1
    
    colramp          <- colorRampPalette(c("yellow","red"))
    data$cols        <- grDevices::densCols(data$x,data$y,colramp = colramp)
       
    cells.number <- sum(cells.number.by.sample)
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(paste0(" BiplotViewer (",format(cells.number, big.mark = " "), " cells)", sep = "")) + #paste0("Sample used : ",paste0(names(samples[samples == TRUE]), collapse = ", "))
            ggplot2::geom_point(ggplot2::aes_string(x = "x",y = "y",colour = "cols"),size = 0.25) +
            ggplot2::stat_density2d(ggplot2::aes_string(x = "x",y = "y"), size = 0.2, colour = "blue", linetype = "dashed") +
            ggplot2::scale_color_identity() +
            ggplot2::xlab(x.marker) +
            ggplot2::ylab(y.marker) +
            ggplot2::coord_cartesian(xlim = c(-1, x.max), ylim = c(-1, y.max)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "black",linetype = "dotted"))
    if (!sample.merge){
        plot <- plot + ggplot2::facet_wrap(~ facet, scales = "free")
    }
    
    return(plot)
    
}