#' @title Generate report and plot.
#'
#' @description Generate a report based on the SPADE result object.
#'
#' @param Results a SPADEResuts or Result object
#' @param reports a vector of the plot names to add in the report (following the vector order). Plot names are : "pheno", "MDS", "cluster", "tree", "disto","stream". By default all plots will be add
#' @param clusters a character vector of ID cluster to be reported
#' @param markers a character vector of markers to be reported
#' @param stat.objects a list of stat.object (object of type DEC, AC or CC)
#' @param profile.objects a list of profile.objects (object of type PhenoProfiles or EnrichmentProfiles)
#' @param PDFfile a character specifying the output path
#' @param width a numeric specifying the plot width
#' @param height a numeric specifying the plot height
#'
#' @import gridExtra
#' 
#' @export
generateReport <- function(Results,
                           PDFfile,
                           reports         = c("pheno", "kinetic", "cluster", "kinetic_cluster", "tree", "disto", "stream", "MDS"),
                           clusters        = NULL,
                           markers         = NULL,
                           assignments     = NULL,
                           stat.objects    = list(),
                           profile.objects = list(),
                           width           = 29.7,
                           height          = 21){
    
    message("[BEGIN]-report")

    print(Results)
    
    plots  <- list()

    if (is.element(c("kinetic-cluster", "kinetic"), reports) && is.null(assignments)){
        stop("Error in generateReport : 'kinetic-cluster' and/or 'kinetic' report required assignments")
    }
    
    for(i in 1:length(reports)){

        print(reports[i])
        switch(reports[i],
               MDS             = {
                  plots <- c(plots, list(MDSViewer(Results, space = "clusters", clusters = clusters)))
                  if (!is.null(assignments)){
                      plots <- c(plots, list(MDSViewer(Results, space = "samples", clusters = clusters, assignments = assignments)))
                  }
               },
               pheno           = {
                  plots <- c(plots, list(ggheatmap.plot(phenoViewer(Results))))
               },
               tree            = {
                  plots <- c(plots, list(treeViewer(Results)))
               },
               kinetic_cluster = {
                   kinetics.plots <- kineticsViewer(Results, clusters = clusters, assignments = assignments)
                   cluster.plots  <- clusterViewer(Results, clusters = clusters, markers = markers)
                   
                   for (i in 1:length(kinetics.plots)){
                       plots <- c(plots, list(gridExtra::grid.arrange(kinetics.plots[[i]], cluster.plots[[i]], ncol = 2)))
                   }
               },
               kinetic         = {
                  plots <- c(plots,kineticsViewer(results, assignments = assignments, clusters = clusters))
               },
               cluster         = {
                  plots <- c(plots,clusterViewer(Results,clusters = clusters, markers = markers))
               },
               disto           = {
                   plots <- c(plots, list(distogramViewer(Results)))
               },
               stream          = {
                   plots <- c(plots, list(streamgraphViewer(Results, clusters = clusters)))
               })
       
    }
    print(length(stat.objects))
    for(stat.object in stat.objects){
        plots <- c(plots, list(plot(stat.object)))
    }
    for(profile.object in profile.objects){
        plots <- c(plots, list(plot(profile.object)))
    }
    
    print(length(plots))
    pages.plots <- gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol=1)

    ggplot2::ggsave(PDFfile, pages.plots, width = width, height = height)
    message("[END]-report")
    
}

