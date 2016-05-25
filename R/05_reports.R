#' @title Generate a report including SPADEVizR plots.
#'
#' @description 
#' Generate a customizable PDF report based on SPADEVizR vizualisation features.
#' Available plots are :
#' \itemize{
#' \item {"count" (included by default):}{Display an representation showing the number of cells for each cluster}
#' \item {"tree" (included by default):}{Display a tree representation showing combined SPADE trees}
<<<<<<< HEAD
#' \item {"heatmap" (included by default):}{Display an heatmap representation}
#' \item {"boxplot":}{Display a boxplot representation. This plot required to provide the 'conditions' parameter}
#' \item {"kinetics":}{Display a kinetic representation for each cluster. This plot required to provide the 'assignments' parameter}
#' \item {"stream":}{Display a streamgraphViewer representation showing the evolution of cells abundance. The 'clusters' parameter is required}
#' \item {"pheno" (included by default):}{Display a parallel coordinate representation showing for each cluster the marker median expression}
#' \item {"MDSclusters" (included by default):}{Display the cluster similarities using MDS}
#' \item {"MDSsamples":}{Display the samples similarities using MDS}
=======
#' \item {"pheno" (included by default):}{Display an heatmap representation}
#' \item {"boxplot":}{Display a boxplot representation. This plot required to provide the 'conditions' parameter}
#' \item {"kinetics":}{Display a kinetic representation for each cluster. This plot required to provide the 'assignments' parameter}
#' \item {"stream":}{Display a streamgraphViewer representation showing the evolution of cells abundance}
#' \item {"cluster" (included by default):}{Display a parallel coordinate representation showing for each cluster the marker median expression}
#' \item {"MDSclusters" (included by default):}{Display the cluster similarities using MDS}
#' \item {"MDSsamples":}{Display the samples similarities using MDS. This plot required to provide the 'assignments' parameter}
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
#' \item {"disto" (included by default):}{Display a distogram representation showing the marker co-expressions}
#' \item {"kinetics_cluster":}{Display a kinetic representation and a parallel coordinate juxtaposed (are arranged one on the side of the other) for each cluster}
#' \item {"boxplot_cluster":}{Display a boxplot representation and a parallel coordinate juxtaposed (are arranged one on the side of the other) for each cluster}
#' }
#' 
#' @param Results a 'SPADEResults' or 'Result' object
#' @param PDFfile a character specifying the output path
#' @param plot.names a character vector specifying the names (see details) and the order of the desired plots
#' @param clusters a character vector of clusters to include in the report (all will be included by default)
#' @param markers a character vector of markers to include in the report (all will be included by default)
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param conditions conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample or NA to exclude
<<<<<<< HEAD
#' @param stat.objects a vector of plotable objects to be displayed in the report (object of class 'DEC', 'AC', 'CC' or 'CCR' accepted)
=======
#' @param stat.objects a vector of stat.object to be displayed in the report (object of class 'DEC', 'AC', 'CC' or 'CCR')
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
#' @param width a numeric specifying the plot width in centimeter
#' @param height a numeric specifying the plot height in centimeter
#'
#' @import gridExtra
#' 
#' @export
generateReport <- function(Results,
                           PDFfile         = "report.pdf",
<<<<<<< HEAD
                           plot.names      = c("count", "heatmap", "tree", "disto", "MDSclusters", "pheno"),
=======
                           plot.names      = c("count", "pheno", "tree", "disto", "MDSclusters", "cluster"),
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
                           clusters        = NULL,
                           markers         = NULL,
                           assignments     = NULL,
                           conditions      = NULL,
                           stat.objects    = list(),
                           width           = 29.7,
                           height          = 21){
    
    message("[BEGIN] - report")

    #print(results)
    
    plots  <- list()

<<<<<<< HEAD
    if (is.element(c("kinetics_pheno", "kinetics"), plot.names) && is.null(assignments)){
        stop("Error in generateReport : 'kinetics_pheno' and/or 'kinetics' report required assignments")
    }
    
    if (is.element(c("boxplot_pheno", "boxplot"), plot.names) && is.null(conditions)){
        stop("Error in generateReport : 'boxplot_pheno' and/or 'boxplot' report required conditions")
=======
    if (is.element(c("kinetics_cluster", "kinetics"), plot.names) && is.null(assignments)){
        stop("Error in generateReport : 'kinetics_cluster' and/or 'kinetics' report required assignments")
    }
    
    if (is.element(c("boxplot_cluster", "boxplot"), plot.names) && is.null(conditions)){
        stop("Error in generateReport : 'boxplot_cluster' and/or 'boxplot' report required conditions")
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
    }
    
    for(i in 1:length(plot.names)){
        switch(plot.names[i],
               "count"           = {
<<<<<<< HEAD
               plots <- c(plots, list(countViewer(Results, clusters = clusters, show.on_device = FALSE)))
               },
               "MDS_clusters"    = {
               plots <- c(plots, list(MDSViewer(Results, space = "clusters", clusters = clusters, show.on_device = FALSE)))
               },
               "MDS_samples"     = {
               plots <- c(plots, list(MDSViewer(Results, space = "samples", clusters = clusters, assignments = assignments, show.on_device = FALSE)))
               },
               "heatmap"           = {
               plots <- c(plots, list(heatmapViewer(Results, show.on_device = FALSE)))
               },
               "tree"            = {
               plots <- c(plots, list(treeViewer(Results, show.on_device = FALSE)))
               },
               "kinetics_pheno" = {
                   kinetics.plots <- kineticsViewer(Results, clusters = clusters, assignments = assignments, show.on_device = FALSE)
                   cluster.plots  <- phenoViewer(Results, clusters = clusters, markers = markers, show.on_device = FALSE)
=======
                   plots <- c(plots, list(countViewer(Results, clusters = clusters)))
               },
               "MDS_clusters"    = {
                  plots <- c(plots, list(MDSViewer(Results, space = "clusters", clusters = clusters)))
               },
               "MDS_samples"     = {
                  plots <- c(plots, list(MDSViewer(Results, space = "samples", clusters = clusters, assignments = assignments)))
               },
               "pheno"           = {
                  plots <- c(plots, list(phenoViewer(Results)))
               },
               "tree"            = {
                  plots <- c(plots, list(treeViewer(Results)))
               },
               "kinetics_cluster" = {
                   kinetics.plots <- kineticsViewer(Results, clusters = clusters, assignments = assignments)
                   cluster.plots  <- clusterViewer(Results, clusters = clusters, markers = markers)
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
                   
                   for (i in 1:length(kinetics.plots)){
                       plots <- c(plots, list(gridExtra::arrangeGrob(kinetics.plots[[i]], cluster.plots[[i]], ncol = 2)))
                   }
               },
<<<<<<< HEAD
               "boxplot_pheno" = {
                   boxplot.plots <- boxplotViewer(Results, clusters = clusters, conditions = conditions, show.on_device = FALSE)
                   cluster.plots <- phenoViewer(Results, clusters = clusters, markers = markers, show.on_device = FALSE)
=======
               "boxplot_cluster" = {
                   boxplot.plots <- boxplotViewer(Results, clusters = clusters, conditions = conditions)
                   cluster.plots <- clusterViewer(Results, clusters = clusters, markers = markers)
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
                   
                   for (i in 1:length(boxplot.plots)){
                       plots <- c(plots, list(gridExtra::arrangeGrob(boxplot.plots[[i]], cluster.plots[[i]], ncol = 2)))
                   }
               },
               "boxplot"         = {
<<<<<<< HEAD
                   plots <- c(plots, boxplotViewer(Results, clusters = clusters, conditions = conditions, show.on_device = FALSE))
               },
               "kinetics"         = {
                   plots <- c(plots, kineticsViewer(Results, assignments = assignments, clusters = clusters, show.on_device = FALSE))
               },
               "cluster"         = {
                   plots <- c(plots, phenoViewer(Results, clusters = clusters, markers = markers, show.on_device = FALSE))
               },
               "disto"           = {
                   plots <- c(plots, list(distogramViewer(Results, show.on_device = FALSE)))
               },
               "stream"          = {
                   plots <- c(plots, list(streamgraphViewer(Results, clusters = clusters, show.on_device = FALSE)))
               },
               "count"           = {
                   plots <- c(plots, list(countViewer(Results, clusters = clusters, show.on_device = FALSE)))
=======
                   plots <- c(plots, boxplotViewer(Results, clusters = clusters, conditions = conditions))
               },
               "kinetics"         = {
                   plots <- c(plots, kineticsViewer(Results, assignments = assignments, clusters = clusters))
               },
               "cluster"         = {
                   plots <- c(plots, clusterViewer(Results,clusters = clusters, markers = markers))
               },
               "disto"           = {
                   plots <- c(plots, list(distogramViewer(Results)))
               },
               "stream"          = {
                   plots <- c(plots, list(streamgraphViewer(Results, clusters = clusters)))
               },
               "count"           = {
                   plots <- c(plots, list(countViewer(Results, clusters = clusters)))
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
               })
       
    }

    for(stat.object in c(stat.objects)){
<<<<<<< HEAD
        plots <- c(plots, list(plot(stat.object, show.on_device = FALSE)))
=======
        plots <- c(plots, list(plot(stat.object)))
>>>>>>> 43ace80f59a8e1b7a65c6f88e2a134a1e91a2f23
    }

    pages.plots <- gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol = 1)

    ggplot2::ggsave(PDFfile, pages.plots, width = width, height = height, unit = "cm")
    message("[END] - report")

}

