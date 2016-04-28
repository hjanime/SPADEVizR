#' @title Textual preview for Results objects
#'
#' @description Prints a preview for a Results object.
#'
#' @param x a Results object
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


#' @title Textual preview for SPADEResults objects
#'
#' @description Prints a preview for a SPADEResults object.
#'
#' @param x a SPADEResults object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

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

#' @title Textual preview for Abundant Clusters (AC) object
#'
#' @description Prints a preview for a Abundant Clusters (AC) object.
#'
#' @param x a Abundant Clusters (AC) object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","AC",
        function(x){
            cat("Object class: Abundant Clusters (AC)\n")
            cat(paste0("Samples : "))
            cat(paste0(x@sample.names, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent : ", x@use.percentages))
            cat("\n")
            cat(paste0("Statiscal test used is : ", x@method))
            cat("\n")
            cat(paste0("Adjusted : ", x@method.adjust))
            cat("\n")
            cat(paste0("pvalue threshold : "))
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("mean threshold : "))
            cat(paste0(" ", x@th.mean, collapse="\n"))
            cat("\n")
            cat(paste0("see @result slot to get further informations...")) 
            cat("\n")     
        }
)


#' @title Textual preview for Differentially Enriched Clusters (DEC) object
#'
#' @description Prints a preview for a Differentially Enriched Clusters (DEC) object.
#'
#' @param x a Differentially Enriched Clusters (DEC) object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","DEC",
        function(x){
            cat("Object class: Differentially Enriched Clusters (DEC)\n")
            cat(paste0("Sample of Condition 1 : "))
            cat(paste0(x@sample.cond1, collapse="; "))
            cat("\n")
            cat(paste0("Sample of Condition 2 : "))
            cat(paste0(x@sample.cond2, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent : ", x@use.percentages))
            cat("\n")
            cat(paste0("Statiscal test used is : ", x@method))
            cat("\n")
            cat(paste0("Adjusted : ", x@method.adjust))
            cat("\n")
            cat(paste0("Paired : ", x@method.paired))
            cat("\n")
            cat(paste0("pvalue threshold : "))
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("foldchange threshold : "))
            cat(paste0(" ", x@th.fc, collapse="\n"))
            cat("\n")
            cat(paste0("see @result slot to get further informations...")) 
            cat("\n")   
        }
)

#' @title Textual preview for Correlated Clusters (CC) objects
#'
#' @description Prints a preview for a Correlated Clusters (CC) object.
#'
#' @param x a Correlated Clusters (CC) object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","CC",
        function(x){
            cat("Object class: Correlated Clusters (CC)\n")
            cat(paste0("Samples : "))
            cat(paste0(x@sample.names, collapse="; "))
            cat("\n")
            cat(paste0("Phenotypic variables : "))
            cat(paste0(x@variable, collapse="; "))
            cat("\n")
            cat(paste0("Use matrix of percent : ", x@use.percentages))
            cat("\n")
            cat(paste0("Statiscal test used is : ", x@method))
            cat("\n")
            cat(paste0("Adjusted : ", x@method.adjust))
            cat("\n")
            cat(paste0("pvalue threshold : "))
            cat("\n")
            cat(paste0(" ", x@th.pvalue, collapse="\n"))
            cat("\n")
            cat(paste0("correlation threshold : "))
            cat("\n")
            cat(paste0(" ", x@th.correlation, collapse="\n"))
            cat("\n")
            cat(paste0("see @result slot to get further informations...")) 
            cat("\n")
        }
)

#' @title Textual preview for PhenoProfiles objects
#'
#' @description Prints a preview for a PhenoProfiles object.
#'
#' @param x a PhenoProfiles object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","PhenoProfiles",
        function(x){
            cat("Object class: PhenoProfiles\n")
            cat(paste0("Clustering method used : "))
            cat(paste0(x@method, collapse="; "))
            cat("\n")
            cat(paste0("Parameter used : "))
            cat(paste0(names(x@method.parameter), " = ", x@method.parameter, collapse="; "))
            cat("\n")
            cat(paste0("Number of cluster : "))
            cat(paste0(x@cluster.number, collapse="; "))
            cat("\n")
            cat(paste0("Number of class : "))
            cat(paste0(x@class.number, collapse="; "))
            cat("\n")
            cat(paste0("see @classes slot to get further informations...")) 
            cat("\n")
        }
)


#' @title Textual preview for EnrichmentProfiles objects
#'
#' @description Prints a preview for a EnrichmentProfiles object.
#'
#' @param x a EnrichmentProfiles object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","EnrichmentProfiles",
        function(x){
            cat("Object class: EnrichmentProfiles\n")
            cat(paste0("Clustering method used : "))
            cat(paste0(x@method, collapse="; "))
            cat("\n")
            cat(paste0("Parameter used : "))
            cat(paste0(names(x@method.parameter), " = ", x@method.parameter, collapse="; "))
            cat("\n")
            cat(paste0("Number of cluster : "))
            cat(paste0(x@cluster.number, collapse="; "))
            cat("\n")
            cat(paste0("Number of class : "))
            cat(paste0(x@class.number, collapse="; "))
            cat("\n")
            cat(paste0("see @classes slot to get further informations...")) 
            cat("\n")
        }
)



#' @title Textual preview for Results objects
#'
#' @description Show a preview a Results object.
#'
#' @param x a Results object
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

#' @title Textual preview for SPADEResults objects
#'
#' @description Show a preview a SPADEResults object.
#'
#' @param x a SPADEResults object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","SPADEResults",
        definition = function(object){print(object)}
)

#' @title Textual preview for Abundant Clusters (AC) object
#'
#' @description Show a preview a for Abundant Clusters (AC) object
#'
#' @param x a Abundant Clusters (AC) object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","AC",
        definition = function(object){print(object)}
)

#' @title Textual preview for Differentially Enriched Clusters (DEC) object
#'
#' @description Show a preview a for Differentially Enriched Clusters (DEC) object
#'
#' @param x a Differentially Enriched Clusters (DEC) object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","DEC",
        definition = function(object){print(object)}
)


#' @title Textual preview for Correlated Clusters (CC) objects
#'
#' @description Show a preview a Correlated Clusters (CC) object.
#'
#' @param x a Correlated Clusters (CC) object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","CC",
        definition = function(object){print(object)}
)

#' @title Textual preview for PhenoProfiles objects
#'
#' @description Show a preview a PhenoProfiles object.
#'
#' @param x a PhenoProfiles object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","PhenoProfiles",
        definition = function(object){print(object)}
)


#' @title Textual preview for EnrichmentProfiles objects
#'
#' @description Show a preview a EnrichmentProfiles object.
#'
#' @param x a EnrichmentProfiles object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","EnrichmentProfiles",
        definition = function(object){print(object)}
)