#' @title Textual previews for all SPADEVizR objects
#'
#' @description Prints a previews for a SPADEVizR object.
#'
#' @param x a SPADEVizR object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","Results",
        function(x){
            cat("Object class: Results\n")
            cat(paste0("Markers : "))
            cat(paste0(x@marker.names, collapse="; "))
            cat("\n")
            cat(paste0("Samples : "))
            cat("\n")
            cat(paste0(" ", x@sample.names, collapse="\n"))
            cat("\n")
        }
)


#' @rdname print-methods
#' @export
setMethod("print","SPADEResults",
	function(x){
	cat("Object class: SPADEResults\n")
	cat(paste0("Markers : "))
	cat(paste0(x@marker.names, collapse="; "))
	cat("\n")
	cat(paste0("Clustering Markers : "))
	cat(paste0(x@marker.names[x@marker.clustering], collapse="; "))
	cat("\n")
	cat(paste0("Samples : "))
    cat("\n")
	cat(paste0(" ", x@sample.names, collapse="\n"))
	cat("\n")
	}
)

#' @rdname print-methods
#' @export
setMethod("print","AC",
        function(x){
            cat("Object class: Abundant Clusters (AC)\n")
            cat(paste0("Samples: "))
            cat(paste0(x@sample.names, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent: ", x@use.percentages))
            cat("\n")
            cat(paste0("Statistical test used is: ", x@method))
            cat("\n")
            cat(paste0("Adjusted: ", x@method.adjust))
            cat("\n")
            cat(paste0("P-value threshold: "))
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("Mean threshold: "))
            cat(paste0(" ", x@th.mean, collapse="\n"))
            cat("\n")    
        }
)

#' @rdname print-methods
#' @export
setMethod("print","DEC",
        function(x){
            cat("Object class: Differentially Enriched Clusters (DEC)\n")
            cat(paste0("Sample of Condition 1: "))
            cat(paste0(x@sample.cond1, collapse="; "))
            cat("\n")
            cat(paste0("Sample of Condition 2: "))
            cat(paste0(x@sample.cond2, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent: ", x@use.percentages))
            cat("\n")
            cat(paste0("Statistical test used is: ", x@method))
            cat("\n")
            cat(paste0("Adjusted: ", x@method.adjust))
            cat("\n")
            cat(paste0("Paired: ", x@method.paired))
            cat("\n")
            cat(paste0("P-value threshold: "))
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("Fold-change threshold: "))
            cat(paste0(" ", x@th.fc, collapse="\n"))
            cat("\n") 
        }
)

#' @rdname print-methods
#' @export
setMethod("print","CC",
        function(x){
            cat("Object class: Correlated Clusters (CC)\n")
            cat(paste0("Samples: "))
            cat(paste0(x@sample.names, collapse="; "))
            cat("\n")
            cat(paste0("Phenotypic variables: "))
            cat(paste0(x@variable, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent: ", x@use.percentages))
            cat("\n")
            cat(paste0("Statistical test used is: ", x@method))
            cat("\n")
            cat(paste0("Adjusted : ", x@method.adjust))
            cat("\n")
            cat(paste0("P-value threshold: "))
            cat("\n")
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("Correlation threshold: "))
            cat("\n")
            cat(paste0(" ", x@th.correlation, collapse="\n"))
            cat("\n")
        }
)

#' @rdname print-methods
#' @export
setMethod("print","PhenoProfiles",
        function(x){
            cat("Object class: PhenoProfiles\n")
            cat(paste0("Number of class: "))
            cat(paste0(x@class.number, collapse = "; "))
            cat("\n")
            cat(paste0("Classification method used: "))
            cat(paste0(x@method, collapse="; "))
            cat("\n")
            cat(paste0("Parameter used"))
            cat(paste0(names(x@method.parameter), " = ", x@method.parameter, collapse="; "))
            cat("\n")
        }
)

#' @rdname print-methods
#' @export
setMethod("print","EnrichmentProfiles",
        function(x){
            cat("Object class: EnrichmentProfiles\n")
            cat(paste0("Number of class: "))
            cat(paste0(x@class.number, collapse = "; "))
            cat("\n")
            cat(paste0("Classification method used: "))
            cat(paste0(x@method, collapse="; "))
            cat("\n")
            cat(paste0("Parameter used"))
            cat(paste0(names(x@method.parameter), " = ", x@method.parameter, collapse="; "))
            cat("\n")
        }
)



#' @title Textual previews for SPADEVizR objects
#'
#' @description Show a previews for a SPADEVizR object.
#'
#' @param object a SPADEVizR object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","Results",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","SPADEResults",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","AC",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","DEC",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","CC",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","PhenoProfiles",
        definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","EnrichmentProfiles",
        definition = function(object){print(object)}
)