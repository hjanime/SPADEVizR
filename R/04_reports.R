#' @title Generate a report including SPADEVizR plots.
#'
#' @description 
#' Generate a customizable PDF report based on SPADEVizR vizualisation features.
#' Available plots are :
#' \itemize {
#' \item "pheno" (included by default): Display an heatmap representation
#' \item "kinetic": Display a kinetic representation for each cluster. This plot required to provide the 'assignment' parameter.
#' \item "cluster" (included by default): Display a parallel coordinate representation showing for each cluster the marker median expression.
#' \item "kinetic_cluster": Display a kinetic representation and a parallel coordinate juxtaposed (are arranged one on the side of the other) for each cluster
#' \item "tree" (included by default): Display a tree representation showing combined SPADE trees. 
#' \item "disto" (included by default): Display a distogram representation showing the marker co-expressions.
#' \item "stream" (included by default): Display a 
#' \item "MDS_clusters" (included by default): Display a Multidimensional Scaling (MDS) representation showing the 
#' \item "MDS_samples": Display a Multidimensional Scaling (MDS) representation showing the. This plot required to provide the 'assignment' parameter.
#' }
#' 
#' @param Results a 'SPADEResults' or 'Result' object
#' @param PDFfile a character specifying the output path
#' @param plots.names a character vector specifying the names (see details) and the order of the desired plots
#' @param clusters a character vector of clusters to include in the report (all will be included by default)
#' @param markers a character vector of markers to include in the report (all will be included by default)
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param stat.objects a vector of stat.object to be displayed in the report (object of class 'DEC', 'AC' or 'CC')
#' @param profile.objects a vector of profile.objects to be displayed in the report (object of class 'PhenoProfiles' or 'EnrichmentProfiles')
#' @param width a numeric specifying the plot width
#' @param height a numeric specifying the plot height
#'
#' @import gridExtra
#' 
#' @export
generateReport <- function(results,
                           PDFfile,
                           plots.names     = c("pheno", "cluster", "tree", "disto", "stream", "MDS_clusters"),
                           clusters        = NULL,
                           markers         = NULL,
                           assignments     = NULL,
                           stat.objects    = list(),
                           profile.objects = list(),
                           width           = 29.7,
                           height          = 21){
    
    message("[BEGIN]-report")

    #print(results)
    
    plots  <- list()

    if (is.element(c("kinetic-cluster", "kinetic"), plots) && is.null(assignments)){
        stop("Error in generateReport : 'kinetic-cluster' and/or 'kinetic' report required assignments")
    }
    
    for(i in length(plots.names)){
        switch(plots.names[i],
               MDS_clusters    = {
                  plots <- c(plots, list(MDSViewer(results, space = "clusters", clusters = clusters)))
               },
               MDS_samples     = {
                  plots <- c(plots, list(MDSViewer(results, space = "samples", clusters = clusters, assignments = assignments)))
               },
               pheno           = {
                  plots <- c(plots, list(phenoViewer(results)))
               },
               tree            = {
                  plots <- c(plots, list(treeViewer(results)))
               },
               kinetic_cluster = {
                   kinetics.plots <- kineticsViewer(results, clusters = clusters, assignments = assignments)
                   cluster.plots  <- clusterViewer(results, clusters = clusters, markers = markers)
                   
                   for (i in 1:length(kinetics.plots)){
                       plots <- c(plots, list(gridExtra::arrangeGrob(kinetics.plots[[i]], cluster.plots[[i]], ncol = 2)))
                   }
               },
               kinetic         = {
                   plots <- c(plots, kineticsViewer(results, assignments = assignments, clusters = clusters))
               },
               cluster         = {
                   plots <- c(plots, clusterViewer(results,clusters = clusters, markers = markers))
               },
               disto           = {
                   plots <- c(plots, list(distogramViewer(results)))
               },
               stream          = {
                   plots <- c(plots, list(streamgraphViewer(results, clusters = clusters)))
               },
               count           = {
                   plots <- c(plots, list(countViewer(results, clusters = clusters)))
               })
       
    }

    for(stat.object in stat.objects){
        plots <- c(plots, list(plot(stat.object)))
    }
    for(profile.object in profile.objects){
        plots <- c(plots, list(plot(profile.object)))
    }
    
    print(length(plots))
    
    plot(plots)
    pages.plots <- gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol = 1)

    ggplot2::ggsave(PDFfile, pages.plots, width = width, height = height)
    message("[END]-report")

}

