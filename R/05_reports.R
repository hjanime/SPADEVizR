#' @title Generate a report including SPADEVizR plots.
#'
#' @description 
#' Generate a customizable PDF report based on SPADEVizR vizualisation features.
#' Available plots are :
#' \itemize{
#' \item {"count" (included by default):}{Display an representation showing the number of cells for each cluster}
#' \item {"tree" (included by default):}{Display a tree representation showing combined SPADE trees}
#' \item {"heatmap" (included by default):}{Display an heatmap representation}
#' \item {"boxplot":}{Display a boxplot representation. This plot required to provide the 'conditions' parameter}
#' \item {"kinetics":}{Display a kinetic representation for each cluster. This plot required to provide the 'assignments' parameter}
#' \item {"stream":}{Display a streamgraphViewer representation showing the evolution of cells abundance. The 'clusters' parameter is required}
#' \item {"pheno" (included by default):}{Display a parallel coordinate representation showing for each cluster the marker median expression}
#' \item {"MDSclusters" (included by default):}{Display the cluster similarities using MDS}
#' \item {"MDSsamples":}{Display the samples similarities using MDS}
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
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param conditions conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample or NA to exclude
#' @param stat.objects a vector of plotable objects to be displayed in the report (object of class 'DAC', 'AC', 'CC' or 'CCR' accepted)
#' @param width a numeric specifying the plot width in centimeter
#' @param height a numeric specifying the plot height in centimeter
#'
#' @import gridExtra
#' 
#' @export
generateReport <- function(Results,
                           PDFfile         = "report.pdf",
                           plot.names      = c("count", "heatmap", "tree", "disto", "MDSclusters", "pheno"),
                           clusters        = NULL,
                           markers         = NULL,
                           samples         = NULL,
                           assignments     = NULL,
                           conditions      = NULL,
                           stat.objects    = list(),
                           width           = 30,
                           height          = 15,
                           verbose         = TRUE) {
    
    message("[BEGIN] - report")

    #print(results)
    
    plots  <- list()

    if (is.element(c("kinetics_pheno", "kinetics"), plot.names) && is.null(assignments)){
        stop("Error in generateReport : 'kinetics_pheno' and/or 'kinetics' report required assignments")
    }
    
    if (is.element(c("boxplot_pheno", "boxplot"), plot.names) && is.null(conditions)){
        stop("Error in generateReport : 'boxplot_pheno' and/or 'boxplot' report required conditions")
    }

    nb.plot <- length(plot.names) + length(stat.objects)
    i       <- 0   

    if (length(plot.names)) {
        for (i in 1:length(plot.names)) {
            if (verbose) {
                message(paste0("\t Generate: ", plot.names[i],", ", i, " on ", nb.plot))
            }
            switch(plot.names[i],
               "count" = {
                   plots <- c(plots, list(countViewer(Results, samples = samples, clusters = clusters, show.on_device = FALSE)))
               },
               "MDS_clusters" = {
                   plots <- c(plots, list(MDSViewer(Results, space = "clusters", clusters = clusters, show.on_device = FALSE)))
               },
               "MDS_samples" = {
                   plots <- c(plots, list(MDSViewer(Results, space = "samples", clusters = clusters, assignments = assignments, show.on_device = FALSE)))
               },
               "heatmap" = {
                   plots <- c(plots, list(heatmapViewer(Results, show.on_device = FALSE)))
               },
               "tree" = {
                   plots <- c(plots, list(treeViewer(Results, samples = samples, show.on_device = FALSE)))
               },
               "kinetics_pheno" = {
                   kinetics.plots <- kineticsViewer(Results, clusters = clusters, assignments = assignments, show.on_device = FALSE)
                   cluster.plots  <- phenoViewer(Results, samples = samples, clusters = clusters, markers = markers, show.on_device = FALSE)

                   for (j in 1:length(kinetics.plots)) {
                       if (verbose) {
                           message(paste0("\tCluster ", j, " on ", length(kinetics.plots)))
                       }
                       plots <- c(plots, list(gridExtra::arrangeGrob(kinetics.plots[[j]], cluster.plots[[j]], ncol = 2)))
                   }
               },
               "boxplot_pheno" = {
                   boxplot.plots <- boxplotViewer(Results, clusters = clusters, conditions = conditions, show.on_device = FALSE)
                   cluster.plots <- phenoViewer(Results, samples = samples, clusters = clusters, markers = markers, show.on_device = FALSE)

                   for (j in 1:length(boxplot.plots)) {
                       if (verbose) {
                            message(paste0("\tCluster ", j, " on ", length(boxplot.plots)))
                       }
                       plots <- c(plots, list(gridExtra::arrangeGrob(boxplot.plots[[j]], cluster.plots[[j]], ncol = 2)))
                   }
               },
               "boxplot" = {
                       plots <- c(plots, boxplotViewer(Results, clusters = clusters, conditions = conditions, verbose = verbose, show.on_device = FALSE))
               },
               "kinetics" = {
                       plots <- c(plots, kineticsViewer(Results, assignments = assignments, clusters = clusters, verbose = verbose, show.on_device = FALSE))
               },
               "pheno" = {
                       plots <- c(plots, phenoViewer(Results, samples = samples, clusters = clusters, markers = markers, verbose = verbose, show.on_device = FALSE))
               },
               "disto" = {
                       plots <- c(plots, list(distogramViewer(Results, samples = samples, show.on_device = FALSE)))
               },
               "stream" = {
                       plots <- c(plots, list(streamgraphViewer(Results, samples = samples, clusters = clusters, show.on_device = FALSE)))
               })

        }
    }
    for (stat.object in stat.objects) {
        if (verbose) {
            i <- i + 1
            message(paste0("\t Generate: ", names(stat.object), ", ", i, " on ", nb.plot))
        }
        plots <- c(plots, list(plot(stat.object, show.on_device = FALSE)))
    }

    pages.plots <- gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol = 1)

    ggplot2::ggsave(PDFfile, pages.plots, width = width, height = height, unit = "cm")
    message("[END] - report")

}

