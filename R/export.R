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
            cat(file = filename, "#sample.names    ", paste0("\"", object@sample.names, "\"", collapse = "	"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#marker.names    ", paste0("\"", object@marker.names, "\"", collapse = "	"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#cluster.number    ", paste0(object@cluster.number, collapse = "	"), sep = "", "\n", append = TRUE)
            
            if (names(object) == "SPADEResults"){
                cat(file = filename, "#use.raw.medians    ", paste0("\"", object@use.raw.medians, "\"", collapse = "	"), sep = "", "\n", append = TRUE)
                cat(file = filename, "#marker.clustering    ", paste0(object@marker.clustering, collapse = "	"), sep = "", "\n", append = TRUE)
                cat(file = filename, "#fcs.files    ", paste0("\"",object@fcs.files, "\"", collapse = "	"), sep = "", "\n", append = TRUE)
            }
            cat(file = filename, "#cells.count below :\n", append = TRUE)
            write.table(object@cells.count, file = filename, append = TRUE, sep = "\t", col.names = NA)
            cat(file = filename, "#marker.expressions below :\n", append = TRUE)
            write.table(object@marker.expressions, file = filename, append = TRUE, sep = "	", col.names = NA)
            if (names(object) == "SPADEResults"){
                cat(file = filename, "#quantiles below :\n", append = TRUE)
                write.table(object@quantiles, file = filename, append = TRUE, sep = "	", col.names = NA)
            }

        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("AC"),
        function(object,filename){
            cat(file = filename, "#Object AC\n")
            cat(file = filename, "#sample.names    ", paste0("\"", object@sample.names, "\"", collapse = "	"), sep = "", "\n", append = TRUE)
            cat(file = filename, "#cluster.size    ", paste0(object@cluster.size, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages    ", paste0(object@use.percentages, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method    ", paste0(object@method, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust    ", paste0(object@method.adjust, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue    ", paste0(object@th.pvalue, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.mean    ", paste0(object@th.mean, collapse = "	"), "\n", sep = "", append = TRUE) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "	",  :  appending column names to file
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "	", col.names = NA)
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("DAC"),
        function(object,filename){
            cat(file = filename, "#Object DAC\n")
            cat(file = filename, "#sample.cond1    ", paste0("\"", object@sample.cond1, "\"", collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#sample.cond2    ", paste0("\"", object@sample.cond2, "\"", collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#cluster.size    ", paste0(object@cluster.size, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages    ", paste0(object@use.percentages, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method    ", paste0(object@method, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust    ", paste0(object@method.adjust, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.paired    ", paste0(object@method.paired, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue    ", paste0(object@th.pvalue, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.fc    ", paste0(object@th.fc, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "	", col.names = NA) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "	",  :  appending column names to file
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("CC"),
        function(object,filename){
            cat(file = filename, "#Object CC\n")
            cat(file = filename, "#sample.names    ", paste0("\"", object@sample.names,"\"", collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#variable    ", paste0(object@variable, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#cluster.size    ", paste0(object@cluster.size, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#use.percentages    ", paste0(object@use.percentages, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method    ", paste0(object@method, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.adjust    ", paste0(object@method.adjust, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.pvalue    ", paste0(object@th.pvalue, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#th.correlation    ", paste0(object@th.correlation, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#result below :\n", append = TRUE)
            write.table(object@result, file = filename, append = TRUE, sep = "	", col.names = NA) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "	",  :  appending column names to file
        }
)

#' @rdname export-methods
#' @export
setMethod("export",c("CCR"),
        function(object,filename){
            cat(file = filename, "#Object CCR\n")
            cat(file = filename, "#type    ", paste0(object@type, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#class.number    ", paste0(object@class.number, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method    ", paste0("\"", object@method, "\"", collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#method.parameter    ", paste0(object@method.parameter, collapse = "	"), "\n", sep = "", append = TRUE)
            cat(file = filename, "#classes below :\n", append = TRUE)
            write.table(object@classes, file = filename, append = TRUE, sep = "	", col.names = NA) # Warning message: In write.table(object@result, file = filename, append = TRUE, sep = "	",  :  appending column names to file
        }
)
