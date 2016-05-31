#' @title Exportation of SPADEVizR objects
#'
#' @description Exports a SPADEVizR object into a tab separated file.
#'
#' @param object a SPADEVizR object
#' @param filename a character indicating the location of output file
#'
#' @return none
#'
#' @name export
#' @rdname export-methods
setGeneric("export", function(object, filename = "export.txt") { standardGeneric("export") })

#' @rdname export-methods
#' @export
setMethod("export",c("Results"),
        function(object,filename){
            
            print(filename)
            cat(file = filename, paste0("#Object Results", names(object),"\n"))
            cat(file = filename, "#sample.names\t", paste0("\"", object@sample.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#marker.names\t", paste0("\"", object@marker.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#cluster.number\t", paste0(object@cluster.number, collapse = "\t"), sep = "", "\n", append = TRUE)
            
            if (names(object) == "SPADEResults"){
                cat(file = filename, "#use.raw.medians\t", paste0("\"", object@use.raw.medians, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
                cat(file = filename, "#marker.clustering\t", paste0(object@marker.clustering, collapse = "\t"), sep = "", "\n", append = TRUE)
                cat(file = filename, "#fcs.files\t", paste0("\"",object@fcs.files, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
            }
            cat(file = filename, "#cells.count below :\n", append = TRUE)
            write.table(object@cells.count, file = filename, append = TRUE, sep = "\t", col.names = NA)
            cat(file = filename, "#marker.expressions below :\n", append = TRUE)
            write.table(object@marker.expressions, file = filename, append = TRUE, sep = "\t", col.names = NA)
            if (names(object) == "SPADEResults"){
                cat(file = filename, "#quantiles below :\n", append = TRUE)
                write.table(object@quantiles, file = filename, append = TRUE, sep = "\t", col.names = NA)
            }

        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("AC"),
        function(object,filename){
            cat(file = filename, "#Object AC\n")
            cat(file = filename, "#sample.names\t", paste0("\"", object@sample.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#cluster.size\t", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages\t", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method\t", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust\t", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue\t", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.mean\t", paste0(object@th.mean, collapse = "\t"), "\n", sep = "", append = TRUE) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "\t",  :  appending column names to file
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "\t", col.names = NA)
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("DAC"),
        function(object,filename){
            cat(file = filename, "#Object DAC\n")
            cat(file = filename, "#sample.cond1\t", paste0("\"", object@sample.cond1, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#sample.cond2\t", paste0("\"", object@sample.cond2, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#cluster.size\t", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages\t", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method\t", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust\t", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.paired\t", paste0(object@method.paired, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue\t", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.fc\t", paste0(object@th.fc, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "\t", col.names = NA) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "\t",  :  appending column names to file
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("CC"),
        function(object,filename){
            cat(file = filename, "#Object CC\n")
            cat(file = filename, "#sample.names\t", paste0("\"", object@sample.names,"\"", collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#variable\t", paste0(object@variable, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#cluster.size\t", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages\t", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method\t", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust\t", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue\t", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.correlation\t", paste0(object@th.correlation, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "\t", col.names = NA) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "\t",  :  appending column names to file
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("CCR"),
        function(object,filename){
            cat(file = filename, "#Object CCR\n")
            cat(file = filename, "#type\t", paste0(object@type, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#class.number\t", paste0(object@class.number, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method\t", paste0("\"", object@method, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.parameter\t", paste0(object@method.parameter, collapse = "\t"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#classes below :\n", append = TRUE)
            write.table(object@classes, file = filename, append = TRUE, sep = "\t", col.names = NA)
        }
)
