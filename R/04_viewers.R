#' @title Visualization of cluster sizes
#'
#' @description 
#' Generate a two dimensional vizualisation showing the number of cells (sum of selected samples) of each cluster.
#' 
#' @param Results a SPADEResults or Results object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param min.cells a numeric specifying the minimum number of cell (sum of all selected samples) to display a cluster
#' @param sort a logical specifying if clusters will be to be sorted (descending) based on the sum of all selected samples for each cluster
#' @param show.samples a logical specifying if the number of cells for all selected samples will be displayed
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @import ggplot2 reshape2
#' 
#' @export
countViewer <- function(Results,
                        samples        = NULL,
                        clusters       = NULL,
                        min.cells      = 0,
                        sort           = TRUE,
                        show.samples   = TRUE,
                        show.on_device = TRUE) {
    
    if(is.null(samples)){ 
        data <- Results@cells.count
    }else{
        data  <- Results@cells.count[, samples, drop = FALSE]
    }
    
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
    
    data.melted           <- reshape2::melt(data, id = c("cluster"))  
    colnames(data.melted) <- c("cluster", "sample", "value")
    data.melted$total     <- ifelse(data.melted[, "sample"] == "sum.of.samples", "sum of selected samples","")
    
    plot <- ggplot2::ggplot(data = data.melted) +
            ggplot2::ggtitle(paste("Count Viewer (", format(cells.number, big.mark=" "), " cells)", sep = ""))
    if(show.samples){                
        plot <- plot + ggplot2::geom_jitter(data   = subset(data.melted, sample != "sum.of.samples"),
                                            ggplot2::aes_string(x = "cluster", y = "value", size = "value", fill = "sample"),
                                            height = 0,
                                            width  = 0.5,
                                            shape  = 21,
                                            alpha  = 0.4)
    }

    step <- 1000

    cluster.maxsize <- max(data.melted$value)
    dot_size.breaks <- seq(0, round(ceiling(cluster.maxsize / step) * step), by = step)

    plot <- plot + ggplot2::geom_point(data = subset(data.melted, sample == "sum.of.samples"),
                                       ggplot2::aes_string(x = "cluster", y = "value", size = "value", shape = "total"),
                                       fill = "grey40") +
                   ggplot2::scale_size(name = "number.of.cells", breaks = dot_size.breaks, range = c(0, 10)) +
                   ggplot2::scale_shape_manual(values = 21) +
                   ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1 * max(data.melted$value))) +
                   ggplot2::ylab("# of cells") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
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
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param marker a character specifying the marker name to display
#' @param highlight an AC, DAC or CC object to highligth identified significant clusters in the SPADE tree
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a list of 'ggplot' objects
#'
#' @import reshape2 ggplot2 grid igraph
#' @importFrom igraph permute normalize tree parent
#' @export
treeViewer <- function(SPADEResults,
                       samples        = NULL,
                       highlight      = NULL,
                       marker         = NULL,
                       show.on_device = TRUE) {
    
    if (names(SPADEResults) == "Results"){
        stop("Error : treeViewer required a SPADEResults object")
    }
    
    data   <- SPADEResults@cells.count
    
    if(!is.null(samples)){ 
        data   <- data[ , samples, drop = FALSE]
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
            expr   <- subset(SPADEResults@marker.expressions, sample %in% samples, drop = FALSE)
        }else{
            expr   <- SPADEResults@marker.expressions
        }

        print(head(expr))

        expr           <- expr[,c("cluster", "sample", marker)]
        expr <- reshape2::dcast(expr, cluster ~ sample, value.var = marker)
        
        rownames(expr) <- expr$cluster
        expr           <- expr[, -1]
        mean.expr      <- apply(expr, 1, mean, na.rm = TRUE)

        pos.vertex     <- cbind(pos.vertex, marker = mean.expr)
        colnames(pos.vertex)[ncol(pos.vertex)] <- marker

        max.mean.expr <- ceiling(max(mean.expr, na.rm = TRUE))
        seq.mean.expr <- seq(from = -1, to = max.mean.expr, by = 1)

    }   
    
    edges    <- igraph::get.edgelist(SPADEResults@graph,names = FALSE)
    pos.edge <- data.frame(x    = SPADEResults@graph.layout[edges[, 1], 1],
                           xend = SPADEResults@graph.layout[edges[, 2], 1],
                           y    = SPADEResults@graph.layout[edges[, 1], 2],
                           yend = SPADEResults@graph.layout[edges[, 2], 2])
    
    cells.number <- sum(colSums(data))
    
    plot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste("Tree Viewer (", format( cells.number, big.mark = " "), " cells)", sep = "")) +
            ggplot2::geom_segment(data = pos.edge, ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"))
    
    if (!is.null(highlight)) {
        highlight.name <- names(highlight)
        pos.vertex[, highlight.name] <- highlight@result$significant
        if(!is.null(marker)){
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker, colour = highlight.name),
                                               stroke = 2.5,
                                               shape = 21) +
                           ggplot2::scale_fill_gradient(low = "#ECE822", high = "#EE302D", limits = c(-1, max.mean.expr), breaks = seq.mean.expr)
        }else{
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", colour = highlight.name),
                                               fill   = "grey80",
                                               stroke = 2.5,
                                               shape  = 21)
        }
    }else{
        if(!is.null(marker)){
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker),
                                               stroke = 2.5,
                                               shape = 21) +
                           ggplot2::scale_fill_gradient(low = "#ECE822", high = "#EE302D", limits = c(-1, max.mean.expr), breaks = seq.mean.expr) #low = yellow, high = red
        }else{
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size"),
                                               fill   = "grey80",
                                               stroke = 2.5,
                                               shape  = 21)
        }
    }
    plot <- plot + ggplot2::scale_color_manual(values = c("black", "blue")) +
                   ggplot2::scale_size_area(max_size = 15) +
                   ggrepel::geom_label_repel(data          = pos.vertex, 
                                             ggplot2::aes_string(x = "x", y = "y", label = "id"),
                                             size          = 4,
                                             color         = "black",
                                             box.padding   = grid::unit(0.1, "lines"),
                                             point.padding = grid::unit(0.1, "lines")) +
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
                                  legend.position = "right")

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
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
#' If the 'Results' parameter is a 'SPADEResults' object, markers used by SPADE to clustered cell populations are shown in bold.
#'
#' @param Results a SPADEResults or Results object
#' @param num a numeric value specifying the number of markers expression categories to use
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a list of 'ggplot' objects
#'
#' @import reshape2
#' 
#' @export
heatmapViewer <- function(Results,
                          num = 5,
                          show.on_device = TRUE) {
    
    pheno.table <- computePhenoTable(Results, num)
    
    pheno.table <- reshape2::dcast(pheno.table, cluster ~ marker)
    cluster     <- pheno.table$cluster
    pheno.table <- pheno.table[, 2:ncol(pheno.table)]
    pheno.table <- t(pheno.table)
    pheno.table <- as.matrix(pheno.table)
    
    colnames(pheno.table) <- cluster
    print(pheno.table)

    if(names(Results) == "SPADEResults"){
        clustering.markers <- Results@marker.names[Results@marker.clustering]
    }else{
        clustering.markers <- NULL
    }

    plot.elements <- ggheatmap(pheno.table, num = num, clustering.markers = clustering.markers)

    plot <- ggheatmap.plot(plot.elements)

    if (show.on_device) {
        grid::grid.draw(plot)
    }

    invisible(plot)

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
#' @param conditions conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param use.percentages a logical specifying if the visualization must be performed on percentage
#' @param show.legend a logical specifying if the legend must be displayed
#' @param show.violin a logical specifying if the count distribution must be displayed
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#' @param verbose a logical specifying if the details of computation must be printed
#'
#' @return a 'ggplot' object
#' 
#' @import reshape2 ggplot2 gtools
#' 
#' @export
boxplotViewer <- function(Results,
                          conditions,
                          clusters        = NULL,
                          use.percentages = TRUE,
                          show.legend     = FALSE,
                          show.violin     = TRUE,
                          show.on_device  = TRUE,
                          verbose         = FALSE) {

    data <- Results@cells.count[, names(conditions[!is.na(conditions)]), drop = FALSE]
    cells.count <- data

    if(use.percentages){
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent)
        legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
    
    if(is.null(clusters)){
        clusters <- rownames(data)
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The 'clusters' parameter must be a character vector")
        }
        clusters    <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data        <- data[clusters,]
    }else{
        stop("Error : Unknown cluster ")
    }
    
    data <- cbind(cluster = rownames(data), data)
    data.melted <- reshape2::melt(data, id = "cluster")
    
    colnames(data.melted) <- c ("cluster", "sample", "value")
    
    data.melted$cond <- conditions[data.melted$sample]
    data.melted      <- data.melted[!is.na(data.melted$cond),]
    data.melted$cond <- factor(data.melted$cond, levels = gtools::mixedsort(unique(data.melted$cond)))

    plots            <- list()
    
    for (current.cluster in clusters) {
        if (verbose) {
            message(paste0("\tCluster ", current.cluster, " on ", length(clusters)))
        }
        data.temp <- data.melted[data.melted$cluster == current.cluster,]

        max.value <- max(data.temp$value, na.rm = TRUE)
        max.value <- max.value + 0.1*max.value + 1
        
        data.temp$cond <- as.factor(data.temp$cond)
        
        i <- length(plots) + 1
        
        cells.number <- sum(cells.count[current.cluster,])
        
        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "cond", y = "value")) +
                      ggplot2::ggtitle(paste("Boxplot Viewer - cluster ", current.cluster, " (", format(cells.number, big.mark=" "), " cells) ", sep = " ")) +
                      ggplot2::geom_boxplot() +
                      ggplot2::geom_jitter(ggplot2::aes_string(color = "sample"), width = 0.2, show.legend = show.legend)
        
        if(show.violin){
            plots[[i]] <- plots[[i]] + ggplot2::geom_violin(alpha = 0.05, fill = "red", colour = "red")
        }
        
        if(use.percentages){    
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value), breaks = round(seq(0, max.value)), minor_breaks = NULL)
        }else{
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value))
        }
        
        plots[[i]] <- plots[[i]] + ggplot2::ylab(legendy) +
                                   ggplot2::xlab("biological conditions") +
                                   ggplot2::theme_bw() +
                                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
                                                  legend.text = ggplot2::element_text(size = 6))
        
    }

    if (show.on_device) {
        gridExtra::grid.arrange(grobs = plots)
    }

    invisible(plots)
}

#' @title Visualization of cluster enrichment profiles kinetics
#' 
#' @description 
#' Generates a kinetics plot representation showing for each cluster its enrichment profiles at each time-point of each individual.
#' 
#' @details 
#' Time-points are sorted in a way that strings with embedded numbers are in the correct order
#' 
#' @param Results a SPADEResults or Results object
#' @param assignments a 2 column data.frame with the sample names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#' @param verbose a logical specifying if the details of computation must be printed
#'
#' @return a 'ggplot' object
#' 
#' @import reshape2 ggplot2 gtools
#' 
#' @export
kineticsViewer <- function(Results,
                           assignments,
                           clusters        = NULL,
                           use.percentages = TRUE,
                           show.on_device  = TRUE,
                           verbose         = FALSE) {
    
    if(missing(assignments) || is.null(assignments)){
        stop("Error : the 'assignments' parameter is required")   
    }
    
    data <- Results@cells.count[, rownames(assignments), drop = FALSE]
    cells.count <- data
    
    if(use.percentages){
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent)#row.names = rownames(data),
        legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
    
    if(is.null(clusters)){
        clusters <- rownames(data)
    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data     <- data[clusters,]
    }else{
        stop("Error : Unknown cluster ")
    }
    
    data <- cbind(cluster = rownames(data), data)
    data.melted <- reshape2::melt(data, id = "cluster")
    
    colnames(data.melted) <- c ("cluster", "sample", "value")
    
    data.melted$individuals <- assignments[data.melted$sample,'individuals']
    data.melted$timepoints  <- assignments[data.melted$sample,'timepoints']
    data.melted$timepoints  <- factor(data.melted$timepoints, levels = gtools::mixedsort(unique(data.melted$timepoints)))

    plots <- list()
    for (current.cluster in clusters) {
        if (verbose) {
            message(paste0("\tCluster ", current.cluster, " on ", length(clusters)))
        }
        data.temp <- data.melted[data.melted$cluster == current.cluster,]

        max.value <- max(data.temp$value)
        max.value <- max.value + max.value * 0.1 + 1

        i <- length(plots) + 1

        cells.number <- sum(cells.count[current.cluster,])

        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "as.factor(timepoints)", y = "value", group = "individuals", color = "individuals")) +
                      ggplot2::ggtitle(paste("Kinetics Viewer - cluster ", current.cluster, " (", format(cells.number, big.mark = " "), " cells) ", sep = " ")) +
                      ggplot2::geom_line() +
                      ggplot2::geom_point(na.rm = TRUE) +
                      ggplot2::scale_x_discrete(expand = c(0, 0.05))
        if (use.percentages) {
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value), breaks = round(seq(0, max.value)), minor_breaks = NULL)
        } else {
            plots[[i]] <- plots[[i]] + ggplot2::scale_y_continuous(limits = c(0, max.value))
        }

        plots[[i]] <- plots[[i]] + ggplot2::ylab(legendy) +
                ggplot2::xlab("timepoints") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
                               legend.text = ggplot2::element_text(size = 6))
    }

    if (show.on_device) {
        gridExtra::grid.arrange(grobs = plots)
    }

    invisible(plots)
}

globalVariables(c("cluster", "sum.of.samples", "ymax", "ymin", "value", "ybase"))

#' @title Visualization of cluster abundance dynamics
#'
#' @description
#' Generate a streamgraph representation showing the dynamic evolution of the number of cells in clusters across samples.
#' The 'clusters' parameter is required.
#'
#' @details
#' The order of samples in the 'samples' vector corespond to the order where the sample will be displayed
#'
#' @param Results a SPADEResults or Results object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be vizualised
#' @param use.relative a logical specifying if the visualization should be performed on relative abundance
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 	
#' @import data.table reshape2 ggplot2
#' 
#' @export
streamgraphViewer <- function(Results,
                              samples        = NULL,
                              clusters       = NULL,
                              use.relative   = FALSE,
                              show.on_device = TRUE) {
    
    data <- Results@cells.count
    
    if (!is.null(samples)) {
        data <- data[, samples, drop = FALSE]
    }
    
    if(is.null(clusters)){
        stop("Error streamgraphViewer: 'clusters' parameter is required")

    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters <- unique(clusters)
        data     <- data[clusters,]

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
    colnames(melted.data) <- c("cluster", "sample", "value")

    melted.data$cluster <- factor(melted.data$cluster, levels = clusters)
    melted.data <- melted.data[order(melted.data$cluster, decreasing = TRUE),]
    melted.data <- melted.data[order(melted.data$sample),]

    dt.melted.data <- data.table::setDT(melted.data)
    dt.melted.data[, ymax := cumsum(value) - (sum(value) / 2), by = sample]
    dt.melted.data[, ymin := ymax - value, by = sample]
    dt.melted.data[, label := format(round(ymax - min(ymin), 2), big.mark = " "), by = sample]
    dt.melted.data[, ybase := min(ymin), by = sample]
    
    title <- paste("Streamgraph with ",ifelse(use.relative,"relative","absolute")," abundance (", format(cells.number, big.mark=" "), " cells)", sep = "")
    plot  <- ggplot2::ggplot(data = dt.melted.data) +
             ggplot2::ggtitle(title) +
             ggplot2::geom_ribbon(ggplot2::aes_string(x = "sample", ymin = "ymin", ymax = "ymax", group = "cluster", fill = "cluster"), color = "grey40", size = 0.1) +
             ggplot2::geom_point(ggplot2::aes_string(x = "sample", y = "ymax", group = "cluster"), shape = 45) +
             ggplot2::geom_point(ggplot2::aes_string(x = "sample", y = "ybase", group = "cluster"), shape = 45) +
             ggplot2::geom_text(ggplot2::aes_string(x = "sample", y = "ymax", label = "label"), check_overlap = TRUE, angle = 360, hjust = 1.1, size = 3) +
             ggplot2::geom_text(ggplot2::aes_string(x = "sample", y = "ybase"), label = "0", check_overlap = TRUE, angle = 360, hjust = 1.1, size = 3) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.text      = ggplot2::element_text(size = 6),
                            axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1),
                            axis.line        = ggplot2::element_blank(),
                            axis.text.y      = ggplot2::element_blank(),
                            axis.ticks       = ggplot2::element_blank(),
                            axis.title.x     = ggplot2::element_blank(),
                            axis.title.y     = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(),
                            panel.border     = ggplot2::element_blank(),
                            panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
    
}

#' @title Visualization of cluster phenotypes
#' 
#' @description 
#' Generates a parallel coordinate plot representation showing for each cluster the marker median expressions.
#' 
#' @details 
#' The ranges of value between marker bounds (using the 'bounds' slot) will be displayed using a grey ribbon.
#' 
#' The 'show.mean' parameter allows to visualize three kinds of information:
#' \itemize{
#' \item "none" value will show marker median expressions for each selected samples;
#' \item "only" value will show only the mean of median maker expressions for all selected samples (displayed as black dashed line);
#' \item "both" value will show marker median expressions for each selected samples together with the mean of median maker expressions for all selected samples.
#' }
#' 
#' If the 'Results' parameter is a 'SPADEResults' object, markers used by SPADE to clustered cell populations are shown in bold.
#'
#' @param Results a SPADEResults or Result object
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param markers a character vector specifying the markers to be displayed 
#' @param show.mean a character specifying if marker means expression should be displayed, possible value are among : "none", "only" or "both"
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#' @param verbose a logical specifying if the details of computation must be printed
#'
#' @return a list of 'ggplot' objects
#'
#' @import reshape2 ggplot2
#' 
#' @export
phenoViewer <- function(Results,
                        clusters       = NULL,
                        samples        = NULL,
                        markers        = NULL,
                        show.mean      = "both",
                        show.on_device = TRUE,
                        verbose        = FALSE) {
    
    if(show.mean != "none" && show.mean != "both" && show.mean != "only"){
        stop("Error : show.mean must be one of those : 'none', 'both' or 'only' ")
    }
    
    if(is.null(samples)){
        data        <- Results@marker.expressions
        cells.count <- Results@cells.count
    }else{
        data        <- subset(Results@marker.expressions, sample %in% samples, drop = FALSE)
        cells.count <- Results@cells.count[, samples, drop = FALSE]
    }
    
    data <- na.omit(data)# NA values are removed, generate a warning ?
    
    if(!is.null(clusters)){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters.select <- data[,"cluster"] %in% clusters
        data            <- data[clusters.select,]
        cells.count     <- cells.count[clusters,]
    }else{
        clusters <- unique(data[,"cluster"])
    }
    
    if(!is.null(markers)){
        data.keys <- data[,c("sample","cluster")]
        data      <- cbind(data.keys,data[,markers])
    }
    
    if (names(Results) == "SPADEResults"){
        markers            <- colnames(data[, grep ("cluster|sample", colnames(data), invert = TRUE)])
        clustering.markers <- is.element(markers, Results@marker.names[Results@marker.clustering])
        bold.markers       <- ifelse(clustering.markers, "bold", "plain")
    }
    
    data           <- reshape2::melt(data, id = c("sample", "cluster"), stringsAsFactors = FALSE)
    colnames(data) <- c("sample", "cluster", "marker", "value")
    
    for(i in 1:nrow(data)){
        data[i, "lower.bound"] <- Results@bounds[1, as.character(data[i, "marker"])]
        data[i, "upper.bound"] <- Results@bounds[2, as.character(data[i, "marker"])]
    }
    
    plots <- list()
    
    for(current.cluster in clusters){
        if (verbose) {
            message(paste0("\tCluster ", current.cluster, " on ", length(clusters)))
        }
        data.temp  <- data[data["cluster"] == current.cluster,]
        max.value <- max(c(data.temp$value, data.temp$upper.bound))
        min.value <- min(c(data.temp$value, data.temp$lower.bound))
                
        max.value <- max.value + 0.1 * max.value
        min.value <- min.value + 0.1 * min.value
        
        i <- length(plots) + 1
        
        cells.number <- sum(cells.count[current.cluster,])
        
        plots[[i]] <- ggplot2::ggplot(data = data.temp) +
                      ggplot2::ggtitle(paste("Pheno Viewer - cluster ", current.cluster, " (", format(cells.number, big.mark = " "), " cells)", sep = ""))
        
        if(show.mean == "both" || show.mean == "none"){
            plots[[i]] <- plots[[i]] + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "value", group = "sample", color = "sample"),
                                                          size = 1)
        }     
        if(show.mean == "only" || show.mean == "both"){
            wide               <- reshape2::dcast(data.temp,sample ~ marker)
            means              <- apply(wide[, 2:ncol(wide)], 2, mean)
            df.means           <- data.frame(marker = names(means), means = means)
            df.means$show.mean <- rep("show.mean", nrow(df.means))
            
            plots[[i]] <- plots[[i]] + ggplot2::geom_line(data     = df.means,
                                                          ggplot2::aes_string(x = "marker", y = "means", group = "show.mean"),
                                                          linetype = "dashed",
                                                          size     = 1)
        }
        if(show.mean == "only"){
            plots[[i]] <- plots[[i]] + ggplot2::theme(legend.position = "none")
        }
        
        plots[[i]] <- plots[[i]] + ggplot2::geom_ribbon(ggplot2::aes_string(x = "as.numeric(marker)", ymin = "lower.bound", ymax = "upper.bound"),
                                                        alpha = 0.1,
                                                        fill  = "grey20")
        
        
        plots[[i]] <- plots[[i]] + ggplot2::scale_x_discrete(limits = markers) +
                                   ggplot2::scale_y_continuous(limits = c(min.value, max.value), breaks = round(seq(0, max.value, by = 1), 0)) +
                                   ggplot2::theme_bw()         
        
        if (names(Results) == "SPADEResults"){
            plots[[i]] <- plots[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 1, face = bold.markers))                 
        }else{
            plots[[i]] <- plots[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 1))
        }
        
        plots[[i]] <- plots[[i]] + ggplot2::theme(legend.text = ggplot2::element_text(size = 6)) +
                                   ggplot2::ylab("MSI")
    }

    if (show.on_device) {
        gridExtra::grid.arrange(grobs = plots)
    }

    invisible(plots)

}

#' @title Visualization of SPADE cluster or sample similarities using MDS
#'
#' @description 
#' Generate a Multidimensional Scaling (MDS) representation showing the similarities between SPADE results based on theirs abundances.
#' 
#' @details 
#' The 'space' parameter specifying if the cluster or sample similarities will be determined using MDS.
#' Available method for the 'dist.method' parameter are : "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#'
#' @param Results a SPADEResults or Results object
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the biologicial condition and secondly the individuals of the experiment
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param space a character specifying the space ("clusters" or "samples", "cluster" by default)
#' @param dist.method a character string containing the name of the distance measure to use
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a list of 'ggplot' objects
#' 
#' @import data.table MASS ggplot2 ggrepel grDevices
#' 
#' @export
MDSViewer <- function(Results,
                      use.percentages = TRUE,
                      assignments     = NULL,
                      clusters        = NULL,
                      space           = "clusters",
                      dist.method     = "euclidean",
                      show.on_device  = TRUE) {
    
    data        <- Results@cells.count
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

    }else if(all(clusters %in% rownames(data))){
        if (typeof(clusters) != "character"){
            stop("Error : the 'clusters' parameter must be a character vector")
        }
        clusters    <- unique(clusters)
        cells.count <- cells.count[clusters,]
        data        <- data[clusters,]
    }else{
        stop("Error : Unknown cluster ")
    }
    
    data <- cbind(cluster = rownames(data), data)               

    if(space == "samples"){
        if (!is.null(assignments)) {
            colnames(assignments) <- c("biological.conditions","individuals")
            data <- data[, rownames(assignments), drop = FALSE]
        }
        
        data <- t(data[, colnames(data) != "cluster"])
    }else if(space != "clusters"){
        stop("Error in \'space\' parameter")
    }
    
    dist   <- dist(data, method = dist.method)
    fit    <- MASS::isoMDS(dist, k = 2, trace = FALSE)
    stress <- fit$stress
    fit    <- fit$point

    x       <- fit[, 1]
    y       <- fit[, 2]
    
    data_i = data.frame(x = x,y = y)
    
    min.lim <- min(min(x), min(y)) * 1.1
    max.lim <- max(max(x), max(y)) * 1.1
    lim     <- max(abs(min.lim), abs(max.lim))
    min.lim <- -lim
    max.lim <- lim
    
    if (space == "samples") {

        title    <- "MDS at the sample level"
        subtitle <- paste0("Kruskal Stress : ", format(round(stress, 2), nsmall = 2))
        
        plot <- ggplot2::ggplot() +
                ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))

        if (!is.null(assignments)) {
            data_i                       <- cbind(data_i, assignments)
            data_i$individuals           <- as.factor(data_i$individuals)
            data_i$biological.conditions <- factor(data_i$biological.conditions, levels = gtools::mixedsort(unique(data_i$biological.conditions)))
            data.table_i                 <- data.table::data.table(data_i, key = "biological.conditions")
            hulls                        <- data.table_i[, .SD[grDevices::chull(x, y)], by = "biological.conditions"]

            plot <- plot + ggplot2::geom_polygon(data = hulls,
                                                 ggplot2::aes_string(x = "x", y = "y", group = "biological.conditions", fill = "biological.conditions"),
                                                 colour = "black",
                                                 alpha = 0.3) +
                           ggplot2::geom_point(data = data_i, ggplot2::aes_string(x = "x", y = "y", colour = "biological.conditions", shape = "individuals"), size = 4)
        } else {
            data_i$sample <- rownames(data_i)
            
            plot <- plot + ggrepel::geom_label_repel(data = data_i, ggplot2::aes_string(x = "x", y = "y", label = "sample", color = "sample"), size = 5) +
                           ggplot2::geom_point(data = data_i, ggplot2::aes_string(x = "x", y = "y", color = "sample"), size = 4)
                           
       }

        plot <- plot + ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                       ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                       ggplot2::xlim(min.lim, max.lim) +
                       ggplot2::ylim(min.lim, max.lim) +
                       ggplot2::coord_fixed() +
                       ggplot2::theme_bw() +
                       ggplot2::theme(panel.background = ggplot2::element_blank(),
                                      panel.border = ggplot2::element_rect(fill = NA),
                                      axis.text.x = ggplot2::element_blank(),
                                      axis.text.y = ggplot2::element_blank(),
                                      axis.title.x = ggplot2::element_blank(),
                                      axis.title.y = ggplot2::element_blank(),
                                      panel.grid.minor = ggplot2::element_blank(),
                                      panel.grid.major = ggplot2::element_blank(),
                                      axis.ticks = ggplot2::element_blank(),
                                      legend.position = "right")

    } else {
        data_i <- cbind(data_i, cluster = data[, "cluster"])
        data_i$cluster <- as.factor(data_i$cluster)

        title    <- "MDS at the cluster level"
        subtitle <- paste0("Kruskal Stress : ", signif(stress, 2))

        plot <- ggplot2::ggplot(data = data_i) +
                ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
                ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y"),
                                    size = 2) +
                ggrepel::geom_text_repel(ggplot2::aes_string(x = "x", y = "y", label = "cluster"),
                                         size = 5) +
                ggplot2::xlim(min.lim, max.lim) +
                ggplot2::ylim(min.lim, max.lim) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                               panel.border = ggplot2::element_rect(fill = NA),
                               axis.text.x = ggplot2::element_blank(),
                               axis.text.y = ggplot2::element_blank(),
                               axis.title.x = ggplot2::element_blank(),
                               axis.title.y = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               axis.ticks = ggplot2::element_blank(),
                               legend.position = "none")
    }

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}

#' @title biplotViewer
#'
#' @description 
#' Generates a biplot representation with two markers
#'
#' @details 
#' In such representation, each dot corresponds to a cell profile and dots are ploted in a 2-dimentional space corresponding to the selected markers. 
#' When too cells dots are displayed, it can require some seconds. In order to seep up the computation, it is possible to reduce the number of cells displayed (downsampling) using the `resample.ratio` parameter. 
#'
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param x.marker a character indicating the marker name of the first dimension
#' @param y.marker a character indicating the marker name of the second dimension
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be used)
#' @param sample.merge a logical specifying if the selected samples must be merged in a single biplot
#' @param resample.ratio a numeric ratio (between 0 and 1) specifying the downsample ratio to show less dots (or NULL)
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a 'ggplot' object
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
                         resample.ratio = NULL,
                         show.on_device = TRUE) {

    if (names(SPADEResults) == "Results"){
        stop("Error : biplotViewer required a SPADEResults object")
    }  
    if (is.null(SPADEResults@flowset)){
        stop("Error : The 'flowset' slot of the 'SPADEResults' object is not loaded, use the function load.flowSet before using the 'biplotViewer' function")
    }

    message("Biplot computation")

    flowset <- SPADEResults@flowset
    
    x.data <- c()
    y.data <- c()
    facet  <- c()
    flowset.samples <- flowCore::sampleNames(flowset)
    
    plots <- list()
    
    if(is.null(samples)){
        samples <- SPADEResults@sample.names
    }

    for (sample in flowset.samples){

        if (is.element(sample, samples)) {
            
            flowframe <- flowset[[which(sample == flowset.samples),]]
            exprs     <- flowframe@exprs
            if (!is.null(clusters)){
                exprs <- subset(exprs, exprs[, "cluster"] %in% clusters)
            }
            x.data <- c(x.data,exprs[, x.marker])
            y.data <- c(y.data,exprs[, y.marker])

            if (!sample.merge){
                facet  <- c(facet, rep(sample, length(exprs[, x.marker])))
            }
        }    
    }
    
    if(is.null(samples)){
        cells.count <- SPADEResults@cells.count
    }else{
        cells.count <- SPADEResults@cells.count[, samples, drop = FALSE]
    }
    cells.number.by.sample <- colSums(cells.count)

    if (sample.merge){
        data <- data.frame(x = x.data, y = y.data)
    } else {
        data <- data.frame(x = x.data, y = y.data, facet = paste0(facet, " (" , format(cells.number.by.sample[facet], big.mark = " "), " cells)"))
    }
    
    if (!is.null(resample.ratio)){
        if (resample.ratio > 0 && resample.ratio < 1){
            data <- data[sample(nrow(data), round((nrow(data) * resample.ratio))),]
        } else {
            stop("resample.ratio must be > 0 and < 1 or null")
        }
    }
    
    x.max              <- max(data["x"]) * 1.1
    y.max              <- max(data["y"]) * 1.1
    
    colramp          <- grDevices::colorRampPalette(c("yellow", "red"))
    data$cols        <- grDevices::densCols(data$x, data$y, colramp = colramp)
    
    cells.number <- sum(cells.number.by.sample)
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(paste0(" BiplotViewer (", format(cells.number, big.mark = " "), " cells)", sep = "")) +
            ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", colour = "cols"), size = 0.25) +
            ggplot2::stat_density2d(ggplot2::aes_string(x = "x", y = "y"), size = 0.2, colour = "blue", linetype = "dashed") +
            ggplot2::scale_color_identity() +
            ggplot2::xlab(x.marker) +
            ggplot2::ylab(y.marker) +
            ggplot2::coord_cartesian(xlim = c(-1, x.max), ylim = c(-1, y.max)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "black", linetype = "dotted"))
    if (!sample.merge){
        plot <- plot + ggplot2::facet_wrap(~ facet, scales = "free")
    }

    message("done")

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
    
}

#' @title Visualization of marker co-expressions
#'
#' @description 
#' Generate a distogram representation showing the marker co-expressions.
#' 
#' @details 
#' A Pearson correlation matrix is calculated between selected markers using selected clusters and samples.
#' High positive correlated markers are shown by a green tile at their perpendicular intersection. 
#' In the same way, absence of correlation are shown by black tiles and negative correlation by red tiles.
#'
#' @param Results a SPADEResults or Results object
#' @param clusters a character vector containing the clusters names to be use (by default all clusters will be used)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param markers a character vector specifying the markers to be displayed 
#' @param show.on_device a logical specifying if the respresentation will be displayed on device 
#'
#' @return a list of 'ggplot' objects
#' 
#' @import reshape2 ggplot2
#' 
#' @export
distogramViewer <- function(Results,
                            clusters       = NULL,
                            samples        = NULL,
                            markers        = NULL,
                            show.on_device = TRUE) {
    
    if(is.null(samples)){
        data        <- Results@marker.expressions
        cells.count <- Results@cells.count
    }else{
        data        <- subset(Results@marker.expressions, sample %in% samples, drop = FALSE)
        cells.count <- Results@cells.count[, samples, drop = FALSE]
    }
    
    if(!is.null(clusters)){
        if (typeof(clusters) != "character"){
            stop("Error : The clusters parameter must be a character vector")
        }
        clusters.select <- data[, "cluster"] %in% clusters
        data            <- data[clusters.select,]
        cells.count     <- cells.count[clusters,]
    }
    
    data <- data[ , -c(1, 2) ]
    if(!is.null(markers)){
        data <- data[, markers]
    }
    
    data <- na.omit(data)# NA values are removed, generate a warning ?
    
    cormat <- round(cor(data), 2)
    dist   <- as.dist(1 - cormat)
    hc     <- hclust(dist)
    cormat <-cormat[hc$order, hc$order]
    
    cormat[upper.tri(cormat, diag = TRUE)] <- NA
    
    markers          <- colnames(cormat)
    dimnames(cormat) <- NULL
    melted.cormat    <- reshape2::melt(cormat)#TODO maybe use a data.frame
    
    bold.markers <- "plain"
    if (names(Results) == "SPADEResults"){
        clustering.markers <- is.element(markers, Results@marker.names[Results@marker.clustering])
        bold.markers <- ifelse(clustering.markers, "bold", "plain")
    }
    
    plot <- ggplot2::ggplot(data = melted.cormat, ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")) + 
            ggplot2::ggtitle("Distogram of marker phenotypes correlations") +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black", 
                                          midpoint = 0, limit = c(-1, 1), na.value = 'white',
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
            ggplot2::theme(axis.line            = ggplot2::element_blank(),
                           axis.text.x          = ggplot2::element_blank(),
                           axis.text.y          = ggplot2::element_blank(),
                           axis.ticks           = ggplot2::element_blank(),
                           axis.title.x         = ggplot2::element_blank(),
                           axis.title.y         = ggplot2::element_blank(),
                           panel.background     = ggplot2::element_blank(),
                           panel.border         = ggplot2::element_blank(),
                           panel.grid.major     = ggplot2::element_blank(),
                           panel.grid.minor     = ggplot2::element_blank(),
                           plot.background      = ggplot2::element_blank(),
                           legend.justification = c(1, 0),
                           legend.position      = c(0.4, 0.7),
                           legend.direction     = "horizontal") +
            ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7,
                                                           barheight      = 1,
                                                           title.position = "top",
                                                           title.hjust    = 0.5))

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)

}