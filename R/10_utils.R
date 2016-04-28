#' @title Internal - computePhenoTable
#' 
#' @description computePhenoTable
#' 
#' @details xxx
#' 
#' @param SPADEResults a SPADEResults object
#' @param num a numeric indicating the number of "phenotype"
#'  
#' @return xxx
#' 
#' @import gtools
#'
computePhenoTable <- function(SPADEResults, num = 5){
    
    message("[START] - computing PhenoTable")
    
    data.melted            <- data.frame(reshape2::melt(SPADEResults@marker.expressions[colnames(SPADEResults@marker.expressions)],id.vars = c("sample","cluster")))
    
    colnames(data.melted)  <- c("sample","cluster","marker","value")

    data.melted$marker     <- as.vector(data.melted$marker)

    means                  <- plyr::ddply(data.melted, c("cluster","marker"), function(df){mean(df$value, na.rm = TRUE)}) #NA values are removed
    colnames(means)        <- c("cluster","marker","value")

    for(i in 1:nrow(means)){
        cluster <- means[i,"cluster"]
        value   <- means[i,"value"]
        min     <- SPADEResults@quantiles[1,means[i,"marker"]]
        max     <- SPADEResults@quantiles[2,means[i,"marker"]]
        seq     <- seq(from = min,to = max,length.out = num)
        if(!is.na(value)){
            res              <- which.min(abs(value - seq))
            means[i,"value"] <- res
        }else{
            means[i,"value"] <- as.numeric(NA)
        }
    }
    
    means <- means[gtools::mixedsort(colnames(means))]
    message("[END] - computing PhenoTable")
    return(means)
    
}


#' title Internal - Create a list of elements allowing to build a heatmap   
#'
#' @description Use the function ggheatmap.plot to display this heatmap
#'
#' @details 
#'
#' @importFrom ggplot2 reshape2
#' 
#' @param mat a matrix
#' @param dendrogram.type a caracter spycifing the look of dendrograms ("rectangle" or "triangle", "rectangle" by default)
#' @param num xxx
#' @return a list of 3 plots (top dendrogram, right dendrogram, heatmap)
#'
#' @import ggplot2 reshape2 grDevices
ggheatmap <- function(mat, dendrogram.type = "rectangle", num = 5, dists = c("eucledian","eucledian")) {

    colnames(mat) <- 1:ncol(mat)
    
    row.hc <- hclust(dist(mat), "ward.D")# change to eucledian or correlation
    col.hc <- hclust(dist(t(mat)), "ward.D")# change to eucledian or correlation
    
    row.dendro <- ggdendro::dendro_data(as.dendrogram(row.hc),type = dendrogram.type)
    col.dendro <- ggdendro::dendro_data(as.dendrogram(col.hc),type = dendrogram.type)
    
    col.plot <- dendro(col.dendro, col=TRUE)
    row.plot <- dendro(row.dendro, row=TRUE)
    
    col.ord <- match(col.dendro$labels$label, colnames(mat))
    row.ord <- match(row.dendro$labels$label, rownames(mat))
    
    mat.ordered  <- mat[row.ord,col.ord]

    data.frame           <- as.data.frame(mat.ordered)
    data.frame$markers   <- rownames(mat.ordered)
    data.frame$markers   <- with(data.frame, factor(markers, levels=markers, ordered=TRUE))
    melted.data.frame    <- reshape2::melt(data.frame, id.vars="markers")
        
    colfunc <- grDevices::colorRampPalette(c("#FFFFFF","#ECE822","#F9A22B","#EE302D","#A32D33"))#white -> yellow -> orange -> red -> brown
    
    print(head(melted.data.frame))
    
    melted.data.frame$value <- as.factor(melted.data.frame$value)
    
    centre.plot <- ggplot2::ggplot(melted.data.frame, ggplot2::aes_string(x = "variable",y = "markers")) + 
                   ggplot2::geom_tile(ggplot2::aes_string(fill = "value"), colour = "black") +
                   ggplot2::scale_fill_manual(values = colfunc(num), guide = ggplot2::guide_legend(direction = "horizontal", ncol = 10, byrow = TRUE, label.theme = ggplot2::element_text(size = 10, angle = 0), label.position="bottom", label.hjust = 0.5, title.position = "top")) + #low = "white", mid = "yellow", high = "#4d2f00",# #4d2f00 -> dark brown    midpoint = round( (num + 1) / 2 ) , limit = c(1,num+1), na.value = "red"
                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 290, hjust = 0),                                  
                                  legend.text = ggplot2::element_text(size = 4),
                                  panel.background = ggplot2::element_rect("white"))
           
    ret <- list(col=col.plot,row=row.plot,centre=centre.plot)
    invisible(ret)
    
}

#' title Internal - Build a dendrograms plot
#'
#' @description xxx
#'
#' @details xxx
#'
#' @param ddata xxx
#' @param row xxx
#' @param col xxx
#' 
#' @return a dendroplot
#'
#' @import ggplot2 ggdendro
dendro <- function(ddata, row=!col, col=!row) {

    tangle <- if(row) { 0 } else { 90 }

    p <- ggplot2::ggplot() +
         ggplot2::geom_segment(data = ggdendro::segment(ddata),
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
        p <- p + ggplot2::scale_x_continuous(expand=c(0.005,0.005)) +
                 ggplot2::coord_flip()                
    } 
    else {
        p <- p +
                ggplot2::scale_x_continuous(expand=c(0.005,0.005))
    }
    return(p)
}

#' title Internal - g_legend 
#'
#' @description Extract legend from a plot
#'
#' @details xxx
#'
#' @param a.gplot a ggplot2 plot
#' 
#' @return a ggplot2 legend
#'
#' @import ggplot2 gtable
g_legend<-function(a.gplot){
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    tmp <- gtable::gtable_filter(tmp, "guide-box")
    
    return(tmp$grobs[[TRUE]])
}

#' title Internal - g_axis
#'
#' @description Extract axis from a plot
#'
#' @details x and y are mutuality excluded (both can't be TRUE) with priority to x
#'
#' @param a.gplot a ggplot2 plot
#' @param x a logical, TRUE to extract x axis
#' @param y a logical, TRUE to extract y axis
#' 
#' @return a ggplot2 axis
#'
#' @import ggplot2 gtable
g_axis<-function(a.gplot, x =! y, y =! x ){
    if (x){
        name <- "axis-l"
    }
    else {
        name <- "axis-b"
    }
    
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    tmp <- gtable::gtable_filter(tmp, name)
    return(tmp$grobs[[TRUE]])
}

#' title Internal - ggheatmap.plot
#'
#' @description ggheatmap.plot display the heatmap build by ggheatmap
#'
#' @details xxx
#'
#' @param L the list of ggplot object provided by ggheatmap
#' @param col.width size of horizontal dendrogram
#' @param row.width size of vertical dendrogram
#' 
#' @return a ggplot2 axis
#'
#' @import ggplot2 gridExtra grid
ggheatmap.plot <- function(L, col.width=0.15, row.width=0.15) {

    layout <- rbind(c(2,1,NA),
                    c(5,3,4),
                    c(NA,6,NA))
    
    center.withoutlegend = L$centre + ggplot2::theme(axis.line        = ggplot2::element_blank(),
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
    
    ret <- gridExtra::grid.arrange(L$col,g_legend(L$centre),center.withoutlegend,L$row,
                                   g_axis(L$centre, x = TRUE ),
                                   g_axis(L$centre, y = TRUE ),
                                   layout_matrix = layout,
                                   widths = grid::unit(c(col.width,1-(2*col.width),col.width), "null"),
                                   heights = grid::unit(c(row.width,1-(2*row.width),row.width), "null"))

    return(ret)
}

#' title Internal - Transforms data for a steam graph (from by ggTimeSeries : https://github.com/Ather-Energy/ggTimeSeries)
#'
#' @description xxx
#'
#' @details xxx
#' 
#' @param xxx
#' 
#' @return xxx
#'
#' @import data.table ggplot2
statSteamgraph <- ggplot2::ggproto("statSteamgraph",
                                   ggplot2::Stat,
                                   required_aes = c("x", "y", "group"),
                                   setup_params = function(data, params) {
            
            data.table::setDT(data)
            
            datasdorder = data[,
                               list(
                                    sdy = sd(y)
                                ),
                                    by = list(group)
                              ][,
                                 ranksdy := rank(sdy, ties.method = 'random')
                              ][
                                 ranksdy %%2 == 0,
                                 ranksdy := -ranksdy
                              ]
            
            data = merge(data,
                         datasdorder,
                         'group',
                         all = T)
            
            data.table::setkey(data, x, ranksdy)
            data[, ymax := cumsum(y) - (sum(y)/2), by = list(x, PANEL)]
            data[, ymin := ymax - y, by = list(x, PANEL)]
            
            list(
                    overalldata = data,
                    na.rm = T
            )
            
        },
        
        compute_group = function(data, scales, overalldata) {

            data = overalldata[group %in% data$group]
            data.frame(
                    x = data[, x],
                    ymin = data[, ymin],
                    ymax = data[, ymax]
            )
            
        }
)

#' title Internal - Plot a steamgraph (from by ggTimeSeries : https://github.com/Ather-Energy/ggTimeSeries)
#'
#' @description xxx
#'
#' @details xxx
#'
#' @param xxx
#' 
#' @return xxx
#'
#' @import ggplot2
stat_steamgraph <- function(mapping     = NULL,
                            data        = NULL,
                            show.legend = NA,
                            inherit.aes = TRUE,
                            na.rm       = TRUE,
                            ...) {  
    ggplot2::layer(stat = statSteamgraph, data = data, mapping = mapping, geom = 'ribbon',
                   position = 'identity', show.legend = show.legend, inherit.aes = inherit.aes,
                   params = list(na.rm = na.rm, ...)) 
}