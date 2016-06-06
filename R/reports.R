#' @title Generate a report including SPADEVizR plots
#'
#' @description 
#' Generate a customizable PDF report based on SPADEVizR visualization features.
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
#' \item {"kinetics_pheno":}{Display a kinetic representation and a parallel coordinate juxtaposed (are arranged one on the side of the other) for each cluster}
#' \item {"boxplot_pheno":}{Display a boxplot representation and a parallel coordinate juxtaposed (are arranged one on the side of the other) for each cluster}
#' \item {AC, DAC, CC and CCR objects}
#' }
#' 
#' @param Results a 'SPADEResults' or 'Result' object
#' @param PDFfile a character specifying the output path
#' @param select.plots a vector combining character and stat objects ('AC', 'DAC', 'CC' and 'CCR) specifying the order of the desired plots (see details) 
#' @param clusters a character vector of clusters to include in the report (all will be included by default)
#' @param markers a character vector of markers to include in the report (all will be included by default)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param assignments a 2 column data.frame with the samples names in row names providing firstly the time-points (numeric) and secondly the individuals (character) of the experiment
#' @param conditions conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample or NA to exclude
#' @param stat.objects a list containing one or several AC, DEC, CC or CCR objects to plot in the report
#' @param width a numeric specifying the plot width in centimeter
#' @param height a numeric specifying the plot height in centimeter
#' @param verbose a boolean specifying if some verbose messages must be displayed during the generation of the report
#'
#' @return none
#'
#' @export
#' 
#' @import gridExtra
generateReport <- function(Results,
                           PDFfile         = "report.pdf",
                           select.plots    = c("count", "heatmap", "tree", "disto", "MDSclusters", "pheno"),
                           clusters        = NULL,
                           markers         = NULL,
                           samples         = NULL,
                           assignments     = NULL,
                           conditions      = NULL,
                           stat.objects    = list(),
                           width           = 50,
                           height          = 30,
                           verbose         = TRUE) {
    
    message("[BEGIN] - report")

    #print(results)
    
    plots  <- list()

    if (is.element(c("kinetics_pheno", "kinetics"), select.plots) && is.null(assignments)){
        stop("Error in generateReport : 'kinetics_pheno' and/or 'kinetics' report required assignments")
    }
    
    if (is.element(c("boxplot_pheno", "boxplot"), select.plots) && is.null(conditions)){
        stop("Error in generateReport : 'boxplot_pheno' and/or 'boxplot' report required conditions")
    }

    nb.plot <- length(select.plots)

    if (length(select.plots)) {
        for (i in 1:length(select.plots)) {

            current.name <- ifelse(typeof(select.plots[[i]]) == "character", select.plots[[i]], "object")

            if (verbose) {
                message(paste0("     Generate: ",
                               ifelse(current.name == "object", names(select.plots[[i]]), current.name),
                               ", ", i, " on ", nb.plot))
            }
            switch(current.name,
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
                           message(paste0("    Cluster ", j, " on ", length(kinetics.plots)))
                       }
                       plots <- c(plots, list(gridExtra::arrangeGrob(kinetics.plots[[j]], cluster.plots[[j]], ncol = 2)))
                   }
               },
               "boxplot_pheno" = {
                   boxplot.plots <- boxplotViewer(Results, clusters = clusters, conditions = conditions, show.on_device = FALSE)
                   cluster.plots <- phenoViewer(Results, samples = samples, clusters = clusters, markers = markers, show.on_device = FALSE)

                   for (j in 1:length(boxplot.plots)) {
                       if (verbose) {
                            message(paste0("    Cluster ", j, " on ", length(boxplot.plots)))
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
               },
               "object" = {
                       plots <- c(plots, list(plot(select.plots[[i]], show.on_device = FALSE)))
               })

        }
    }

    pages.plots <- gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol = 1)

    ggplot2::ggsave(PDFfile, pages.plots, width = width, height = height, unit = "cm")
    message("[END] - report")

}

