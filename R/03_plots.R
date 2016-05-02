#' @title abundantClustersViewer
#'
#' @description Generate a scatter plot representation displaying for each cluster its m
#' 
#' @details By default, only significant abundant clusters are labeled. Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#' 
#' @param AC an object of class AC (object returned by the 'computeAC()' function)
#' @param show.cluster.sizes a logicial specifing if points size are related to cell count or not (TRUE by default)
#' @param show.all.labels a logicial specifing if all cluster label must be show or just significant cluster
#' 
#' @return a ggplot object
#' 
#' @import ggplot2 ggrepel
#' 
#' @export
abundantClustersViewer <- function(AC,
                                   show.cluster.sizes = TRUE,
                                   show.all.labels    = FALSE){

    AC@result <- cbind (AC@result, cluster.size = AC@cluster.size)
    
    data.text <- AC@result
    if (!show.all.labels){
        data.text <- subset(AC@result, AC@result$significance)
    }

    plot <-  ggplot2::ggplot(data = AC@result) +
             ggplot2::ggtitle("Cells abundance of clusters")

    plot <- plot + ggplot2::geom_hline(yintercept = AC@th.mean,
                                       linetype   = "dashed",
                                       alpha      = 0.3,
                                       color      = "red",
                                       size       = 1) +  
                   ggplot2::geom_vline(xintercept = -log10(AC@th.pvalue),
                                       linetype   = "dashed",
                                       alpha      = 0.3,
                                       color      = "red",
                                       size       = 1)
    if (show.cluster.sizes > 0){
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


#' @title Volcano Plot Viewer
#'
#' @description Generate a Volcano plot representation based on Differentially Enriched Clusters (DEC) of The SPADE result object.
#' 
#' @details By default, only significant differentially enriched clusters are labeled. Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param SPADEResults the SPADEViewer result object
#' @param DEC an object of class DEC (object returned by the 'computeDEC()' function)
#' @param fc.log2 a logicial specifing if fold-change or log2(fold-change) is use 
#' @param show.cluster.sizes a logicial specifing if points size are related to cell count or not (TRUE by default)
#' @param show.all.labels a logicial specifing if all cluster labels must be show or just significant cluster
#'
#' @return a ggplot object
#'  
#' @import ggplot2
#' 
#' @export
volcanoViewer <- function(DEC                = NULL,
                          fc.log2            = TRUE,
                          show.cluster.sizes = TRUE,
                          show.all.labels    = FALSE){
    
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
    if (!show.all.labels){
        data.text <- subset(DEC@result, DEC@result$significance)
    }
    min           <- floor(min(DEC@result$fold.change))
    max           <- ceiling(max(DEC@result$fold.change))
    max           <- max(max,abs(min))
    x.breaks      <- c(round(c(-th.fc, th.fc),2), seq(-max, max, by = 1))
    title.details <- ifelse(DEC@use.percentages,"(with % of cells)","(with # of cells)")
    
    plot <- ggplot2::ggplot(data = DEC@result, ggplot2::aes_string(x = "fold.change", y = "-log10(pvalue)")) +
            ggplot2::ggtitle(paste0("Volcano plot showing differentially enriched clusters ",title.details)) +
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

#' @title correlatedClustersViewer
#'
#' @description 
#' 
#' @details By default, only significant correlated clusters are labeled. Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param CC an object of class CC (object returned by the 'computeCC()' function)
#' @param show.cluster.sizes a logicial specifing if points size are related to cell count or not (TRUE by default)
#' @param show.all.labels a logicial specifing if all cluster label must be show or just significant cluster
#' 
#' @return a ggplot object
#' 
#' @import ggplot2 ggrepel
#' 
#' @export
correlatedClustersViewer <- function(CC,
                                     show.cluster.sizes = TRUE,
                                     show.all.labels    = FALSE){
    
    CC@result <- cbind (CC@result, cluster.size = CC@cluster.size)
    
    data.text <- CC@result
    if (!show.all.labels){
        data.text <- subset(CC@result, CC@result$significance)
    }
    title.details <- ifelse(CC@use.percentages,"(in # of cells)","(in % of cells")
    plot <- ggplot2::ggplot(data = CC@result) +
            ggplot2::ggtitle(paste0("Correlation of clusters kinetics", title.details," with variables provided")) +
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
                    ggplot2::xlab(paste(AC@method,"coeficient of correlation")) +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::theme_bw()
    return(plot)
}

##' @title profilesViewer
##'
##' @description Generate a graph representation of PhenoProfiles classes
##' 
##' @details xxx
##'
##' @param profile.object a PhenoProfiles object or an EnrichmentProfiles object
##' @param show.unclassified a logical specifying if unclassified clusters must be shown or not (FALSE by default)
##' @param cluster.size a logicial specifing if points size are related to cell count or not (TRUE by default)
##' 
##' @return a ggplot object
##' 
##' @import ggplot2 ggnetwork
##' 
##' @export
##' 
#profilesViewer <- function (profile.object,
#                            show.unclassified = FALSE,
#                            cluster.size      = TRUE){
#    
#    classes           <- profile.object@classes
#    
#    if (nrow(classes) < 2){
#        stop("Error : Not enougth clusters classified to plot this profile (at least 2)")
#    }
#    nb.cluster        <- profile.object@cluster.number
#    
#    unclassified <- setdiff(1:nb.cluster,classes[,"cluster"])
#    rest <- data.frame(cluster = unclassified, classe = rep(NA,length(unclassified)))
#    classes <- rbind(classes,rest)
#    
#    classes <- classes[order(classes$cluster),]
#    classes.withoutNA <- na.omit(classes)
#    
##   print(classes.withoutNA)
#    edges <- data.frame()
#    for (i in 1:length(unique(classes.withoutNA$classe))){
#        
#        same.classes <- classes.withoutNA[classes.withoutNA$classe == i,]
#        same.classes <- same.classes[order(same.classes$cluster),]
#        
#        print(same.classes)
#        if (nrow(same.classes) > 1){
#            x      <- c()
#            y      <- c()
#            classe <- c()
#            
#            previous <- NA
#            for (j in 1:nrow(same.classes)){
#                
#                x        <- c(x,same.classes[j,"cluster"])
#                y        <- c(y,previous)
#                previous <- same.classes[j,"cluster"]
#                
#            }
#            y[1] <- previous 
#            print(paste0("classe :",i))
#            print(data.frame(x = x, y = y))
#            
#            if (nrow(edges) > 0){
#                edges <- rbind(edges,data.frame(x = x, y = y))
#            }else{
#                edges <- data.frame(x = x, y = y)
#            }        
#        }
#    }
#        
#    graph <- network::network.initialize(nb.cluster, directed = FALSE)
#    graph <- network::add.edges(graph,edges$x,edges$y)
#    
#    network::set.vertex.attribute(x = graph, attrname = "classe", value = classes$classe)
#    
#    print(graph)
#    
#    if (cluster.size){
#        network::set.vertex.attribute(x = graph, attrname = "size", value = profile.object@cluster.size)
#    }
#    
#    graph <- ggnetwork::ggnetwork(graph, layout = "kamadakawai")#kamadakawai  , kkconst = nb.cluster^3
#    
#    if (!show.unclassified){
#        graph <- graph[!is.na(graph$classe),]
#    }
#    print(graph)
#    plot <- ggplot2::ggplot(graph, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
#            ggnetwork::geom_edges(linetype = "twodash", color = "grey90", size = 1)
#    
#    if (cluster.size){
#        plot <- plot + ggnetwork::geom_nodes(ggplot2::aes_string(color = "as.factor(classe)", size = "size"))
#    }else{
#        plot <- plot + ggnetwork::geom_nodes(ggplot2::aes_string(color = "as.factor(classe)"), size = 8)
#    }
#    
#    plot <- plot + ggnetwork::geom_nodelabel_repel(ggplot2::aes_string(label = "vertex.names"),
#                                                   color       = "black",
#                                                   fontface    = "bold",
#                                                   box.padding = grid::unit(1, "lines")) +
#                   ggnetwork::theme_blank()
#        
#    return(plot)
#    
#}

#' @title profilesViewer
#'
#' @description Generate a graph representation of PhenoProfiles classes
#' 
#' @details xxx
#'
#' @param profile.object a PhenoProfiles object or an EnrichmentProfiles object
#' @param cluster.size a logicial specifing if points size are related to cell count or not (TRUE by default)
#' 
#' @return a ggplot object
#' 
#' @import ggplot2 ggnetwork
#' 
#' @export
#' 
profilesViewer <- function (profile.object,
                             cluster.size      = TRUE){

    classes <- profile.object@classes    
    if (nrow(classes) < 2){
        stop("Error : Not enougth clusters classified to plot this profile (at least 2)")
    }
    
    classes$ID.node   <- 1:nrow(profile.object@classes)
    
    edges <- data.frame()
    for (i in 1:length(unique(classes$classe))){
        
        same.classes <- classes[classes$classe == i,]
        
        #print(same.classes)
        if (nrow(same.classes) > 1){
            x      <- c()
            y      <- c()
            classe <- c()
            
            previous <- NA
            for (j in 1:nrow(same.classes)){
                
                x        <- c(x,same.classes[j,"ID.node"])
                y        <- c(y,previous)
                previous <- same.classes[j,"ID.node"]
                
            }
            y[1] <- previous 
            #print(paste0("classe :",i))
            #print(data.frame(x = x, y = y))
            
            if (nrow(edges) > 0){
                edges <- rbind(edges,data.frame(x = x, y = y))
            }else{
                edges <- data.frame(x = x, y = y)
            }
        }
    }
        
    graph <- network::network.initialize(nrow(profile.object@classes), directed = FALSE)
    graph <- network::add.edges(graph,edges$x,edges$y)
    
    network::set.vertex.attribute(x = graph, attrname = "classe",  value = classes$classe,  v = classes$ID.node)
    network::set.vertex.attribute(x = graph, attrname = "cluster", value = classes$cluster, v = classes$ID.node)
    
    if (cluster.size){
        network::set.vertex.attribute(x = graph, attrname = "size", value = profile.object@cluster.size, v = classes$ID.node)
    }
    #print(graph)
    graph <- ggnetwork::ggnetwork(graph, layout = "kamadakawai")
    
    #print(graph)
    plot <- ggplot2::ggplot(graph, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
            ggnetwork::geom_edges(linetype = "twodash", color = "grey90", size = 1)
    
    if (cluster.size){
        plot <- plot + ggnetwork::geom_nodes(ggplot2::aes_string(color = "as.factor(classe)", size = "size"))
    }else{
        plot <- plot + ggnetwork::geom_nodes(ggplot2::aes_string(color = "as.factor(classe)"), size = 8)
    }
    
    plot <- plot + ggnetwork::geom_nodelabel_repel(ggplot2::aes_string(label = "cluster"),
                                                   color       = "black",
                                                   fontface    = "bold",
                                                   box.padding = grid::unit(1, "lines")) +
                   ggnetwork::theme_blank()
    
    return(plot)
    
}

##' @title xxx
##'
##' @description xxx
##'
##' @details xxx
##'
##' @param xxx
##'
##' @return xxx
##' 
##' @name plot
##' @rdname plot-methods
##' @export 
#setGeneric("plot", function(x,...){ standardGeneric("plot") })

#' @rdname plot-methods
#' @export
setMethod("plot",c("DEC"),
        function(x,...){
            return(volcanoViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("AC"),
        function(x,...){
            return(abundantClustersViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("CC"),
        function(x,...){
            return(correlatedClustersViewer(x,...))
        }
)


#' @rdname plot-methods
#' @export
setMethod("plot",c("PhenoProfiles"),
        function(x,...){
            return(profilesViewer(x,...))
        }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("EnrichmentProfiles"),
        function(x,...){
            return(profilesViewer(x,...))
        }
)

#' @title ClusterViewer
#' 
#' @description Cluster viewer 
#' 
#' @details xxx
#' 
#' 
#' @param Results a SPADEResuts or Result object (without quantiles bounds if a Results object is provided)
#' @param samples a named vector providing the correspondence between samples name (in rowname) and the logical value TRUE to use these samples (all samples by default)
#' @param clusters a character vector containing the clusters to use for the representation
#' @param markers a pattern describing markers to observe
#' @param show.mean a character : "none" "both" "only", "both" by default
#' 
#' @return a list of ggplot objects 
#'
#' @import reshape2 ggplot2
#' 
#' @export
clusterViewer <- function(Results,
                          samples       = NULL,
                          clusters      = NULL,
                          markers       = NULL,
                          show.mean     = "both"){
    
    if(show.mean != "none" && show.mean != "both" && show.mean != "only"){
        stop("Error : show.mean must be one of those : 'none' 'both' 'only' ")
    }
    data <- c()
    if(is.null(samples)){ 
        data        <- Results@marker.expressions
        cells.count <- Results@cells.count
    }else{
        data        <- subset(Results@marker.expressions, sample %in% names(samples[ samples == TRUE]), drop = FALSE)
        cells.count <- Results@cells.count[,c("cluster",names(samples[ samples == TRUE]))]
    }

    if(!is.null(clusters)){
        clusters.select <- data[,"cluster"] %in% clusters
        data            <- data[clusters.select,]
    }else {
        clusters <- data[,"cluster"]
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
            #print(paste0(i,data[i,"marker"]))
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
        
        cells.number <- sum(cells.count[cells.count$cluster == current.cluster,colnames(cells.count) != "cluster"])#-which()
                
        plots[[i]] <- ggplot2::ggplot(data = data.temp) +
                      ggplot2::ggtitle(paste("cluster ",current.cluster," - Cluster Viewer (", cells.number," cells)" ,sep = ""))                      
        
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

#' @title Pheno Viewer
#' 
#' @description Heatmap Viewer aims to xxx
#' 
#' @details xxx
#' 

#' 
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param pheno.table a result of the function computePhenoTable
#' @param num a numeric indicating the precision of computed expression scores
#'
#' @return a list of ggplot objects
#'
#' @import reshape2
#' 
#' @export
phenoViewer <- function(SPADEResults,
                        pheno.table = NULL,
                        num         = 5){
    #print(SPADEResults)
    if (names(SPADEResults) == "Results"){
       stop("Error : phenoViewer required a SPADEResults object")
    }                
                    
    if(is.null(pheno.table)){
        pheno.table <- computePhenoTable(SPADEResults, num)
    }
    
    pheno.table <- reshape2::dcast(pheno.table, cluster ~ marker)
    pheno.table <- pheno.table[,2:ncol(pheno.table)]
    pheno.table <- t(pheno.table)
    pheno.table <- as.matrix(pheno.table)

    plot.elements <- ggheatmap(pheno.table, num = num, clustering.markers = SPADEResults@marker.names[SPADEResults@marker.clustering])

    heatmap <- ggheatmap.plot(plot.elements)
    
    return(heatmap)
    
}

#' @title SPADE Tree viewer
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param SPADEResults a SPADEResults object (Results object is not accepted)
#' @param samples a named vector providing the correspondence between samples name (in rowname) and the logical value TRUE to use these samples for cell counting (all samples by default)
#' @param stat.object an AC, DEC or CC object to highligth significant clusters in the SPADE tree
#' @param vertex_size a numeric vector of two values indicating the range of 
#'
#' @return a list of ggplot objects
#'
#' @import reshape2 ggplot2
#' 
#' @export
treeViewer <- function(SPADEResults,
                       samples       = NULL,
                       stat.object   = NULL,
                       vertex_size   = c(1,15)){
    
    if (names(SPADEResults) == "Results"){
        stop("Error : treeViewer required a SPADEResults object")
    }
                   
    data.filtered   <- SPADEResults@cells.count[, colnames(SPADEResults@cells.count) != "cluster"]
    
    if(!is.null(samples)){ 
        data.filtered   <- data.filtered[names(samples[ samples == TRUE ]), drop = FALSE]
    }else{
        data            <- data.filtered
    }
    
    vertex.size <- apply(data,1,sum)
    
    pos.vertex <- data.frame(id   = as.character(1:nrow(SPADEResults@graph.layout)),
                             x    = SPADEResults@graph.layout[,1],
                             y    = SPADEResults@graph.layout[,2],
                             size = vertex.size)
    edges    <- igraph::get.edgelist(SPADEResults@graph,names = FALSE)
    pos.edge <- data.frame(x    = SPADEResults@graph.layout[edges[,1],1],
                           xend = SPADEResults@graph.layout[edges[,2],1],
                           y    = SPADEResults@graph.layout[edges[,1],2],
                           yend = SPADEResults@graph.layout[edges[,2],2])

    plot <- ggplot2::ggplot() +
            ggplot2::ggtitle("Tree Viewer") +
            ggplot2::geom_segment(data = pos.edge, ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"))
    
    if (!is.null(stat.object)){
        stat.object.name <- names(stat.object)
        pos.vertex[,stat.object.name] <- stat.object@result$significance
        plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size", fill = stat.object.name), shape = 21)  
    }else{
        plot <- plot + ggplot2::geom_point(data = pos.vertex, ggplot2::aes_string(x = "x", y = "y", size = "size"), fill = "#009ACD", shape = 21)
    }
    plot <- plot + ggplot2::scale_fill_manual(values = c("deepskyblue3","firebrick3")) +
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

#' @title Kinetic Viewer
#' 
#' @description Kinetic viewer aim to represent the kinetics of xxx
#' 
#' @details xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param assignments a 2 column data.frame with the samples names in rownames providing firstly the timepoints (numeric) and secondly the individuals (caracter) of the experiment
#' @param clusters a numerical vector containing the clusters to use in the representation, by default all clusters will be use
#' @param use.percentages a logical specifying if the visualisation should be performed on percentage
#' 
#' @return a ggplot object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
kineticsViewer <- function(Results,
                           assignments,
                           clusters        = NULL,
                           use.percentages = TRUE){
    
    data <- Results@cells.count

    if(is.null(clusters)){
        clusters <- data[,"cluster"]
        message("All clusters will be cumpute")
    }
    else if(all(clusters %in% data[,"cluster"])){
        
        clusters <- unique(clusters)
        data     <- subset(data, cluster %in% clusters)
        
        message("These clusters will be cumpute :")
        message(paste0(clusters,collapse="\t"))
        
    }else{
        stop("Error : Unknown cluster ")
    }
    
    cells.count <- data[,c("cluster",rownames(assignments))]
    
    if(use.percentages){
        data.percent <- prop.table(as.matrix(data[,colnames(data) != "cluster"],2)) * 100
        data         <- data.frame(cluster = data[,"cluster"],data.percent)
        legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
            
    data.melted <- reshape2::melt(data, id = "cluster")
    
    colnames(data.melted) <- c ("cluster", "sample", "value")

    data.melted$individuals <- assignments[data.melted$sample,'individuals']
    data.melted$timepoints  <- assignments[data.melted$sample,'timepoints']

    #print(data.melted)
    

    
    plots <- list()
    for(current.cluster in clusters){
        #print(clusters)
        data.temp  <- data.melted[data.melted$cluster == current.cluster,]
        
        max.value  <- max(data.temp$value)
        max.value  <- max.value + max.value*0.1 + 1
        
        i <- length(plots) + 1
        
        cells.number <- sum(cells.count[cells.count$cluster == current.cluster,colnames(cells.count) != "cluster"])
        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "timepoints", y = "value", group = "individuals", color = "individuals")) +
                      ggplot2::ggtitle(paste("cluster ",current.cluster," - Kinetics Viewer (", cells.number," cells)",sep = "")) +
                      ggplot2::geom_line(size = 1, na.rm = TRUE) + #TODO add smooth curve with spline function
                      ggplot2::geom_point(ggplot2::aes_string(x = "timepoints", y = "value", color = "individuals"), na.rm = TRUE)
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

#' @title boxplotViewer
#'
#' @description Generate a boxplot representation to compare the cell enrichment of biological conditions for each cluster
#' 
#' @details xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param conditions conditions a named vector providing the correspondence between a sample name (in rownames) and the condition of this sample : NA to exclude a sample from tests
#' @param clusters a numerical vector containing the clusters to use in the representation
#' @param use.percentages a logical specifying if the visualisation should be performed on percentage
#' @param label a logical to show sample label or not (FALSE by default)
#' @param violin xxx
#'
#' @return a ggplot object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
boxplotViewer <- function(Results,
                          conditions,
                          clusters        = NULL,
                          use.percentages = TRUE,
                          label           = FALSE,
                          violin          = TRUE){
                      
      data <- Results@cells.count
                      
      if(is.null(clusters)){
         clusters <- data[,"cluster"]
         message("All clusters will be cumpute")
      }
      else if(all(clusters %in% data[,"cluster"])){
                          
         clusters <- unique(clusters)
         data     <- subset(data, cluster %in% clusters)
                          
         message("These clusters will be cumpute :")
         message(paste0(clusters,collapse="\t"))
                          
     }else{
         stop("Error : Unknown cluster ")
     }
                      
     if(use.percentages){
         data.percent <- prop.table(as.matrix(data[,colnames(data) != "cluster"],2)) * 100
         data         <- data.frame(cluster = data[,"cluster"],data.percent)
         legendy = "% of cells relative to parent"
    }else{
        legendy = "# of cells"
    } 
    
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
        
        cells.number <- sum(Results@cells.count[Results@cells.count$cluster == current.cluster,])
        plots[[i]] <- ggplot2::ggplot(data = data.temp, ggplot2::aes_string(x = "cond", y = "value")) +
                ggplot2::ggtitle(paste("cluster ",current.cluster," - Biological conditions viewer (", cells.number," cells)",sep = "")) +
                ggplot2::geom_boxplot() +
                ggplot2::geom_jitter(ggplot2::aes_string(color = "sample"), width = 0.2, show.legend = label)
        if(violin){
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

#' @title Distogram Viewer
#'
#' @description Generate a distogram representation showing the marker co-expression.
#' 
#' @details xxx
#'
#' @param Results a SPADEResults or Results object
#' 
#' @return a list of ggplot objects
#' 
#' @import reshape2 ggplot2
#' 
#' @export
distogramViewer <- function(Results){
    
    data <- Results@marker.expressions[,grep("sample|cluster",colnames(Results@marker.expressions), invert = TRUE)]
    
    data <- na.omit(data)# NA values are removed
    
    cormat <- round(cor(data),2)
    dist <- as.dist((1-cormat)/2)
    hc <- hclust(dist)
    cormat <-cormat[hc$order, hc$order]
    
    cormat[upper.tri(cormat,diag = TRUE)] <- NA
    
    markers <- colnames(cormat)
    
    dimnames(cormat) <- NULL
    
    melted.cormat <- reshape2::melt(cormat)#TODO maybe use a data.frame
    
    bold.markers <- "plain"
    if (names(Results) == "SPADEResults"){
        clustering.markers <- is.element(markers, Results@marker.names[Results@marker.clustering])
        bold.markers <- ifelse(clustering.markers,"bold","plain")
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

#' @title streamgraphViewer
#'
#' @description Generate a streamgraph representation showing the evolution of cell count in clusters across samples
#' 
#' @details xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param order a named vector a named vector providing the correspondence between a sample name (in rownames) and an integer ordering samples in numeric (NA to exclude this sample)
#' @param clusters a numerical vector containing the clusters to use in the representation, by default all clusters will be use
#'
#' @return a ggplot object
#' 
#' @import reshape2 ggplot2
#' 
#' @export
streamgraphViewer <- function(Results,
                              order        = NULL,
                              clusters     = NULL){
    
    cells.count    <- Results@cells.count
    
    print(clusters)
    
    if(!is.null(order)){
        cells.count    <- cells.count[,c("cluster",names(order[!is.na(order)]))]
    }
    
    if(is.null(clusters)){
        message("All clusters will be cumpute")
    }
    else if(all(clusters %in% cells.count[,"cluster"])){     

        cells.count  <- subset(cells.count, cluster %in% clusters)
        
        message("These clusters will be cumpute :")
        message(paste0(clusters,collapse="\t"))
        
    }else{
        stop("Error : Unknown cluster ")
    }

    melted.cells.count           <- reshape2::melt(cells.count, id = "cluster")
    colnames(melted.cells.count) <- c("cluster","samples","value")
    
    print(melted.cells.count)
    
    melted.cells.count$cluster <- as.factor(melted.cells.count$cluster)
    
    plot = ggplot2::ggplot(data = melted.cells.count, ggplot2::aes_string(x = "samples", y = "value", group = "cluster", fill = "cluster")) +
           ggplot2::ggtitle("Streamgraph showing the evolution of cell abondance in clusters across samples") +
           stat_steamgraph() +
           ggplot2::theme_bw() +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0),
                          legend.text = ggplot2::element_text(size = 6))
    
    return(plot)
          
}




#' @title MDS viewer.
#'
#' @description Generate a MDS representation based on the SPADE result object.
#' 
#' @details xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param use.percentages a logical specifying if the visualisation should be performed on percentage
#' @param assignments a 2 column data.frame with the samples names in rownames providing firstly the timepoints (numeric) and secondly the individuals (caracter) of the experiment
#' @param clusters specify the set of clusters to use
#' @param space a caracter specifying the space ("clusters" or "samples", cluster by default)
#' @param dist.method a character string containing the name of the distance measure to use
#' 
#' @return a list of ggplot objects
#' 
#' @import MASS ggplot2 ggrepel
#' 
#' @export
MDSViewer <- function(Results,
                      use.percentages = TRUE,
                      assignments,
                      clusters        = NULL,
                      space           = "clusters",
                      dist.method     = "euclidean"){
    
    data <- Results@cells.count
                  
       if(is.null(clusters)){
            clusters <- data[,"cluster"]
            message("All clusters will be cumpute")
        }else if(all(clusters %in% data[,"cluster"])){
             clusters <- unique(clusters)
             data     <- subset(data, cluster %in% clusters)
                      
             message("These clusters will be cumpute :")
             message(paste0(clusters,collapse="\t"))
                      
        }else{
             stop("Error : Unknown cluster ")
        }
                  
        if(use.percentages){
                data.percent <- prop.table(as.matrix(data[,colnames(data) != "cluster"],2)) * 100
                data         <- data.frame(cluster = data[,"cluster"],data.percent)
                legend.y = "% of cells relative to parent"
        }else{
                legend.y = "# of cells"
        }    
    
    if(space == "samples"){
        data       <- t(data[,colnames(data) != "cluster"])
    }else if(space != "clusters"){
        stop("Error in \"space\" parameter")
    }
    
    dist <- dist(data,method = dist.method)
    message("MDS computation")
    
    fit    <- MASS::isoMDS(dist, k = 2)
    stress <- fit$stress
    fit    <- fit$point
    message("done")
    
    x       <- fit[,1]
    y       <- fit[,2]
    
    datai = data.frame(x = x,y = y)
    
    min.lim <- min(min(x),min(y))*1.1
    max.lim <- max(max(x),max(y))*1.1
    lim     <- max(abs(min.lim),abs(max.lim))
    min.lim <- -lim
    max.lim <- lim
    
    
    if(space == "samples"){
        
        datai   <- cbind (datai, assignments)
        samples <- rownames(data)
        
        datai$individuals <- as.factor(datai$individuals)
        datai$timepoints  <- as.factor(datai$timepoints)
        
        plot <- ggplot2::ggplot(data = datai)  +
                ggplot2::ggtitle("MDS(samples)") +
                ggplot2::geom_hline(yintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_vline(xintercept = (min.lim+max.lim)/2,linetype = "dashed") +
                ggplot2::geom_polygon(ggplot2::aes_string(x = "x",y = "y", group = "individuals",fill = "individuals"),colour = "black", alpha = 0.3) +
                ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", colour = "individuals", shape = "timepoints"),size = 4) +
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
                ggplot2::annotate(geom  = "text", x = -Inf, y = -Inf, hjust = -1,vjust = -1, label = paste0("Kruskal Stress : ",round(stress,2)))
        
    }else{
        datai <- cbind (datai, cluster = data[,"cluster"])
        datai$cluster <- as.factor(datai$cluster)
        
        plot <- ggplot2::ggplot(data = datai)  +
                ggplot2::ggtitle("MDS(clusters)") +
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
                ggplot2::annotate(geom  = "text", x = -Inf, y = -Inf, hjust = -1,vjust = -1, label = paste0("Kruskal Stress : ",round(stress,2)))
    }
    
    return(plot)
}

#' @title CountViewer
#'
#' @description Generate a two dimensional vizualisation showing the number of cells (sum of selected samples) of each cluster.
#' 
#' @details xxx
#' 
#' @param Results a SPADEResults or Results object
#' @param samples a named vector providing the correspondence between samples name (in rowname) and the logical value TRUE to use these samples (all samples by default)
#' @param clusters a character vector containing the clusters to use for the representation
#' @param min.cells a numeric specifying the minimun number of cell (sum of all selected samples) to display a cluster
#' @param sort a logical specifying if clusters will be to be sorted (descending) based on the sum of all selected samples for each cluster.
#' @param show.samples a logical specifying if the number of cells for all selected samples will be displayed.
#' 
#' @return a ggplot object
#' 
#' @import ggplot2 reshape2
#' 
#' @export
CountViewer <- function(Results,
                        samples      = NULL,
                        clusters     = NULL,
                        min.cells    = 0,
                        sort         = TRUE,
                        show.samples = TRUE){

    if(is.null(samples)){ 
        data <- Results@cells.count
    }else{
        data  <- Results@cells.count[, c("cluster", names(samples[ samples == TRUE]))]
    }
   
    if(!is.null(clusters)){
        clusters.select <- data[,"cluster"] %in% clusters
        data            <- data[clusters.select,]
    }
    
    count <- data[,colnames(data) != "cluster", drop = FALSE]
    #print(dim(count))
    data  <- cbind(data, sum.of.samples = apply(count,1,sum))
    data  <- data[data$sum.of.samples > min.cells,]
    
    if (sort){
        data <- transform(data, cluster = reorder(cluster, -sum.of.samples))
    }else{
        data$cluster <- as.factor(data$cluster)
    }
    
    data.melted <- reshape2::melt(data, id = c("cluster"))  
    colnames(data.melted) <- c("cluster","sample","value")
    
    data.melted$total <- ifelse(data.melted[,"sample"] == "sum.of.samples","sum of selected samples","")
    
    plot <- ggplot2::ggplot(data = data.melted) +
            ggplot2::ggtitle("CountViewer showing the number of cells for each cluster") +
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

##' @title ClusterViewer
##' 
##' @description Cluster viewer 
##' 
##' @details xxx
##' 
##' @param SPADEResults a SPADEResuts object
##' @param samples a named vector providing the correspondence between samples name (in rowname) and the logical value TRUE to use these samples (all samples by default)
##' @param clusters a numerical vector containing the clusters to use in the representation
##' @param markers a pattern describing markers to observe
##' @param show.mean a character : "none" "both" "only", "both" by default
##' 
##' @return a list of ggplot objects 
##'
##' @import reshape2 ggplot2
##' 
##' @export
#clusterViewer2 <- function(SPADEResults,
#                           samples       = NULL,
#                           clusters      = NULL,
#                           markers       = NULL,
#                           show.mean     = "both"){
#    
#    if(show.mean != "none" && show.mean != "both" && show.mean != "only"){
#        stop("Error : show.mean must be one of those : 'none' 'both' 'only' ")
#    }
#    data <- c()
#    if(is.null(samples)){ 
#        data  <- SPADEResults@marker.expressions
#    }else{
#        data  <- subset(SPADEResults@marker.expressions, name %in% names(samples[ samples == TRUE]) )
#    }
#    clusters.select <- c()
#    if(!is.null(clusters)){
#        clusters.select <- data[,"cluster"] %in% clusters
#        data            <- data[clusters.select,]
#    }
#
#    if(!is.null(markers)){
#        data.marker <- data[,c("sample","cluster")]
#        data <- cbind(data.marker,data[,markers])
#    }
#    markers <- colnames(data[, grep ("cluster|sample",colnames(data), invert = TRUE)])
#    clustering.markers <- is.element(markers,SPADEResults@marker.names[SPADEResults@marker.clustering])
#    bold.markers <- ifelse(clustering.markers,"bold","plain")
#    
#    plots <- list()
#    
#    for(i in 1:length(clusters)){
#
#        data.temp  <- data[data["cluster"] == clusters[i], colnames(data) != "cluster"]
#
#        print(data.temp)
#        
#        plots[[i]] <- GGally::ggparcoord(data = data.temp,  showPoints = FALSE,scale = "globalminmax", scaleSummary = "mean", centerObsID = 1, #missing = "exclude",
#                                         splineFactor = FALSE, alphaLines = 1, boxplot = FALSE, groupColumn = "name", order = "allClass",
#                                         shadeBox = NULL, mapping = NULL, title = paste("cluster ",clusters[i]," - Cluster Viewer",sep = ""))
#                
#    }
#    
#    return(plots)    
#}
#
#
#
#
#
#
