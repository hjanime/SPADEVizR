setOldClass("igraph") # give access to igraph class

#' @title Results class definition
#' 
#' @description Results is a S4 object containing cytometric results.
#' 
#' @details Results is the core (super classe) of SPADEResults object but allow importation of data comming from other software than SPADE. It aims to store the data comming from the importResults() function. 
#' 
#' @slot marker.expressions a dataframe containing the number of cells for each cluster of each sample
#' @slot cells.count a numerical dataframe containing marker median expressions for each cluster of each sample
#' @slot sample.names a character vector containing the sample names
#' @slot marker.names a character vector containing the markers names
#' @slot cluster.number a numeric specyfing the number of clusters
#' 
#' @name Results-class
#' @rdname Results-class
#' @exportClass Results
Results <- setClass("Results",
                    slots=c(cells.count        = "data.frame",
                            marker.expressions = "data.frame",
                            sample.names       = "character",
                            marker.names       = "character",
                            cluster.number     = "numeric"),
                    validity = function(object){
                        if((length(object@marker.names) + 2) != ncol(object@marker.expressions)){
                            message(paste0("Object Results, Error : marker.names length (",
                                            length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistents with marker.expressions size (number of colunms : ",
                                            ncol(object@marker.expressions),")"))
                            return (FALSE)
                        }
                        if(nrow(object@cells.count) != object@cluster.number){
                            message(paste0("Object Results, Error : cluster.number (",
                                            object@cluster.number,
                                            ") is inconsistent with cells.count matrix size (",
                                            nrow(object@cells.count),")"))
                            return (FALSE)
                        }
                    })
            
#' @title SPADEResults class definition
#' 
#' @description SPADEResults is a S4 object containing SPADE results.
#' 
#' @details SPADEResults extend Results object adding specific SPADE data. SPADEResults object is the central element used by the most of SPADEViR functions, it aims to store the SPADE data comming from the importSPADEResults() function. 
#' 
#' @slot use.raw.medians a logical specifying if the marker expressions are obtained from "raw_median" or "transformed" columns
#' @slot dictionary a two column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)
#' @slot marker.clustering a logical vector specifying marker that have been used during the clustering precedure
#' @slot fcs.files a character vector containing the absolute path of the original FCS files
#' @slot quantiles a numerical matrix containing the quantiles for each each markers of cluster
#' @slot graph a igraph object containing the SPADE tree
#' @slot graph.layout a numerical matrix containing the layout of the SPADE tree
#' 
#' @import igraph
#' 
#' @name SPADEResults-class
#' @rdname SPADEResults-class
#' @exportClass SPADEResults
SPADEResults <- setClass("SPADEResults",
                         contains = "Results",
	                     slots=c(use.raw.medians    = "logical",#check the best name
                                 dictionary         = "data.frame",
                                 marker.clustering  = "logical",
                                 flowset            = "flowSet",
                		         fcs.files          = "character",#TODO think about storing flowset rather than fcs.files
                                 quantiles          = "data.frame",
                		         graph              = "igraph",
                		         graph.layout       = "matrix"),
                         validity = function(object){
                                                          
                             if((length(object@sample.names) * object@cluster.number) != nrow(object@marker.expressions)){
                                 message(paste0("Object SPADEResults, Error : sample.names length (",
                                                length(object@sample.names),") and cluster.number ("
                                                ,object@cluster.number,
                                                ") are inconsistents with marker.expressions size (number of row : ",
                                                nrow(object@marker.expressions),")"))
                                 return (FALSE)
                             }
                             if((length(object@marker.names) + 2) != ncol(object@marker.expressions)){
                                 message(paste0("Object SPADEResults, Error : marker.names length (",
                                                length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistents with marker.expressions size (number of colunms : ",
                                                ncol(object@marker.expressions),")"))
                                 return (FALSE)
                             }
                             if(nrow(object@cells.count) != object@cluster.number){
                                 message(paste0("Object SPADEResults, Error : cluster.number (",
                                                 object@cluster.number,
                                                 ") is inconsistent with cells.count matrix size (",
                                                 nrow(object@cells.count),")"))
                                 return (FALSE)
                             }
                             if(ncol(object@cells.count) != (length(object@sample.names)+1)){
                                 message(paste0("Object SPADEResults, Error : number of samples (",
                                                 length(object@sample.names),
                                                 " + 1 for cluster) is inconsistent with cells.count matrix size (",
                                                 ncol(object@cells.count),")"))
                                 return (FALSE)
                             }
                             if (length(object@marker.clustering) > length(object@marker.names)){
                                 message(paste0("Object SPADEResults, Error : marker.clustering length (",
                                                 length(object@marker.clustering),
                                                 ") can't be higher than marker.names length (",
                                                 length(object@marker.names),")"))
                                 return (FALSE)
                             }
                             if (length(object@marker.clustering) == 0){
                                 warning("Object SPADEResults, Warning : marker.clustering length is 0")
                             }
                             for(fcs.file in object@fcs.files){
                                 if(!file.exists(fcs.file)){
                                     message(paste0("Object SPADEResults, Error : file not exist :",fcs.file))
                                     return (FALSE)
                                 }
                             }
                             if (nrow(object@quantiles) != 2){
                                 message(paste0("Object SPADEResults, Error : quantiles number of rows (",
                                                nrow(object@quantiles),
                                                ") is incorrect (only 2 rows accepted)"))
                                 return (FALSE)
                             }
                             if (ncol(object@quantiles) != length(object@marker.names)){
                                 message(paste0("Object SPADEResults, Error : quantiles number of columns (",
                                                 ncol(object@quantiles),
                                                 ") is inconsistent with marker.names length (",
                                                 length(object@marker.names),")"))
                                 print(colnames(object@quantiles))
                                 print(object@marker.names)
                                 return (FALSE)
                             }                             
                            
                             return(TRUE)
                         }
)

#' @title Abundant Clusters (AC) class definition
#' 
#' @description AC is a S4 object containing the result of computeAC function. 
#' 
#' @details AC is a printable and a plotable object calling the abundantClustersViewer() fonction.
#' 
#' @slot sample.names a character vector containing the samples names used
#' @slot cluster.size a numeric vector containing number of cells ( sum of all samples ) for each cluster
#' @slot use.percentages a logical specifying if computations was performed on percentage of cell abondance
#' @slot method a character containing the name of the statistical test used to identify the AC
#' @slot method.adjust a character containing the name of the multiple correction method used
#' @slot th.pvalue a numeric vector with pvalue threshold 
#' @slot th.mean a numeric vector with mean threshold
#' @slot result a data.frame with clusters in row. The first and second colunm contain the mean and the standard deviation and the third contain the pvalue. Finnaly the last columns is a logicial ("significance") describing if the 2 thresholds was reached or not  
#' 
#' @name AC-class
#' @rdname AC-class
#' @exportClass AC
AC <- setClass("AC",
               slots=c(sample.names    = "character",
                       cluster.size    = "numeric",
                       use.percentages = "logical",
                       method          = "character",
                       method.adjust   = "character",
                       th.pvalue       = "numeric",   
                       th.mean         = "numeric",
                       result          = "data.frame"),
               validity = function(object){
            
                   if(length(object@sample.names) == 0){
                       message("Object AC, Error : sample.names length can't be equal to 0")
                       return(FALSE)
                   }
                   if(!any(object@method %in% c("t.test","wilcox.test"))){
                       message("Object AC, Error : method must match 't.test' or 'wilcox.test'")
                       message(paste0("method founded is : ", object@method))
                       return(FALSE)
                   }
                   if(!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))){
                       message("Object AC, Error : method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
                       message(paste0("method.adjust founded is : ", object@method.adjust))
                       return(FALSE)
                   }
                   if(object@th.pvalue < 0 || object@th.pvalue > 1){
                       message("Object AC, Error : th.pvalue must be include into [0,1] interval")
                       message(paste0("th.pvalue founded is : ", object@th.pvalue))
                       return(FALSE)
                   }                
                   if(!identical(colnames(object@result),c("cluster","mean","sd","pvalue","significance"))){
                      print(colnames(object@result))
                      message("Object AC, Error in result slot, result must have this colmuns : 'cluster','mean','sd','pvalue','significance'")
                      message("Colmuns founded are : ")
                      message(paste(colnames(object@result)))
                      return(FALSE)
                   }
                   
                   
            
                   return(TRUE)
               }
)

#' @title Differentially Enriched Clusters (DEC) class definition
#' 
#' @description DEC is a S4 object containing the result of computeDEC function.
#' 
#' @details DEC is a printable and a plotable object calling the volcanoViewer() fonction.  
#' 
#' @slot sample.cond1 a character specifying the names of the samples corresponding the first condition
#' @slot sample.cond2 a character specifying the names of the samples corresponding the second condition
#' @slot cluster.size a numeric vector containing number of cells ( sum of all samples ) for each cluster
#' @slot use.percentages a logical specifying if computations was performed on percentage of cell abondance
#' @slot method a character containing the name of the statistical test used to identify the DEC
#' @slot method.adjust a character containing the name of the multiple correction method used 
#' @slot method.paired a logical indicating if the test has beeen performed in a paired manner
#' @slot th.pvalue a numeric vector with pvalue threshold
#' @slot th.fc a numeric vector with foldchange threshold
#' @slot result a dataframe with clusters in row. The first and second colunm contain the mean and the standard deviation for the first condition, the third and fourth the mean and standard deviation for the second condition, the sixth the fold-change and the seventh the pvalue and finnaly the significance according to pvalue threshold and fold change threshold  
#'
#' @name DEC-class
#' @rdname DEC-class
#' @exportClass DEC
DEC <- setClass("DEC",
        slots=c(sample.cond1    = "character",
                sample.cond2    = "character",
                cluster.size    = "numeric",
                use.percentages = "logical",
                method          = "character",
                method.adjust   = "character",
                method.paired   = "logical",
                th.pvalue       = "numeric",
                th.fc           = "numeric",
                result          = "data.frame"),
        validity = function(object){

            if(length(object@sample.cond1) == 0){
                message("Object DEC, Error : sample.cond1 length can't be equal to 0")
                return(FALSE)
            }
            if(length(object@sample.cond2) == 0){
                message("Object DEC, Error : sample.cond2 length can't be equal to 0")
                return(FALSE)
            }
            if(!any(object@method %in% c("t.test","wilcox.test"))){
                message("Object DEC, Error : method must match 't.test' or 'wilcox.test'")
                message(paste0("method founded is : ", object@method))
                return(FALSE)
            }
            if(!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))){
                message("Object DEC, Error : method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
                message(paste0("method.adjust founded is : ", object@method.adjust))
                return(FALSE)
            }
            if(object@th.pvalue < 0 || object@th.pvalue > 1){
                message("Object DEC, Error : th.pvalue must be include into [0,1] interval")
                message(paste0("th.pvalue founded is : ", object@th.pvalue))
                return(FALSE)
            }
            if(object@th.fc == 0){
                message("Object DEC, Error : th.fc can't be equal to 0")
                message(paste0("th.fc founded is : ", object@th.fc))
                return(FALSE)
            }
            if(!identical(colnames(object@result),c("cluster","mean.cond1","sd.cond1","mean.cond2","sd.cond2","fold.change","pvalue","significance"))){
                print(colnames(object@result))
                message("Object DEC, Error in result slot, result must have this colmuns : 'cluster','mean.cond1','sd.cond1','mean.cond2','sd.cond2','fold.change','pvalue','significance'")
                message("Colmuns found are : ")
                message(paste(colnames(object@result)))
                return(FALSE)
            }
                       
            return(TRUE)
        }
)


#' @title Correlated Clusters (CC) class definition
#' 
#' @description CC is a S4 object containing the result of computeCC function.
#' 
#' @details CC is a printable and a plotable object calling the correlatedClustersViewer() fonction 
#' 
#' @slot cluster.size a numeric vector containing number of cells ( sum of all samples ) for each cluster
#' @slot variable a numeric vector containing the expression values of the variable
#' @slot use.percentages a logical specifying if computations was performed on percentage of cell abondance
#' @slot method a character containing the name of the statistical test used to identify the CC
#' @slot method.adjust a character containing the name of the multiple correction method used 
#' @slot th.pvalue a numeric vector with pvalue threshold
#' @slot th.correlation a numeric vector with correlation threshold (r)
#' @slot result a three colmuns dataframe with clusters in row. The first colunm contains the coefficiant of correlation (r), the second contains the associated pvalue and the third a logical (significance) specifying if the two thresholds was reached or not.
#' 
#' @name CC-class
#' @rdname CC-class
#' @exportClass CC
CC <- setClass("CC",
               slots=c(sample.names    = "character",
                       variable        = "numeric",
                       cluster.size    = "numeric",
                       use.percentages = "logical",
                       method          = "character",
                       method.adjust   = "character",
                       th.pvalue       = "numeric",
                       th.correlation  = "numeric",
                       result          = "data.frame"),
               validity = function(object){#TODO complete this

                   if(length(object@variable) != length(object@sample.names)){
                       message(paste0("Object CC, Error : variable length (",
                                       ,length(object@variable),
                                       ") is inconsistents with sample names length (",
                                       length(object@sample.names),")"))
                       return(FALSE)
                   }
                   if(length(object@sample.names) == 0){
                       message("Object CC, Error : sample.names length can't be equal to 0")
                       return(FALSE)
                   }
                   if(!any(object@method %in% c("pearson", "kendall", "spearman"))){
                       message("Object CC, Error : method must match 'pearson', 'kendall', 'spearman'")
                       message(paste0("method founded is : ", object@method))
                       return(FALSE)
                   }
                   if(!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))){
                       message("Object CC, Error : method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
                       message(paste0("method.adjust founded is : ", object@method))
                       return(FALSE)
                   }
                   if(object@th.pvalue < 0 || object@th.pvalue > 1){
                       message("Object CC, Error : th.pvalue must be include into [0,1] interval")
                       message(paste0("th.pvalue founded is : ", object@th.pvalue))
                       return(FALSE)
                   }
                   if(object@th.correlation < 0 || object@th.correlation > 1){
                       message("Object CC, Error : th.correlation must be include into [0,1] interval")
                       message(paste0("th.correlation founded is : ", object@th.pvalue))
                       return(FALSE)
                   }
                   if(!identical(colnames(object@result),c("cluster","correlation","pvalue","significance"))){
                       print(colnames(object@result))
                       message("Object CC, Error in result slot, result must have this colmuns : 'cluster','correlation','pvalue','significance'")
                       message("Colmuns found are : ")
                       message(paste(colnames(object@result)))
                       return(FALSE)
                   }
                    
                   return(TRUE)
               }
)


#' @title PhenoProfiles class definition
#' 
#' @description PhenoProfiles is a S4 object containing the result of classifyPhenoProfiles() function.
#' It classify clusters by their markers properties provided by the fonction computePhenoTable()
#' 	
#' @details PhenoProfiles is a printable and a plotable object calling the xxx() fonction 
#' 
#' @slot method a character specifying the method used to classify cluster
#' @slot method.parameter a list of parameters used by the selected method
#' @slot cluster.size a numeric vector with the number of cell in each cluster (sum of all samples)
#' @slot cluster.number a numeric providing the number of cluster
#' @slot class.number a numeric providing the number of classes
#' @slot classes a two column dataframe with the cluster in first colunm and corresponding classe in the second colunm 
#' 
#' @name PhenoProfiles-class
#' @rdname PhenoProfiles-class
#' @exportClass PhenoProfiles
PhenoProfiles <- setClass("PhenoProfiles",
                          slots=c(method              = "character",
                                  method.parameter         = "list",
                                  cluster.size             = "numeric",
                                  cluster.number           = "numeric",
                                  class.number             = "numeric",
                                  classes                  = "data.frame"),
                          validity = function(object){#TODO complete this
                                
                             return(TRUE)
                         }
)

#' @title EnrichmentProfiles class definition
#' 
#' @description EnrichmentProfiles is a S4 object containing the result of classifyEnrichmentProfiles() function.
#' It classify clusters by their cells enrichment properties across the samples.
#'  
#' @details EnrichmentProfiles is a printable and a plotable object calling the profileViewer() fonction 
#' 
#' @slot method a character specifying the method used to classify cluster
#' @slot method.parameter a list of parameters used by the selected method
#' @slot cluster.size a numeric vector with the number of cell in each cluster (sum of all samples)
#' @slot cluster.number a numeric providing the number of cluster
#' @slot class.number a numeric providing the number of classes
#' @slot classes a two column dataframe with the cluster in first colunm and corresponding classe in the second colunm 
#' 
#' @name EnrichmentProfiles-class
#' @rdname EnrichmentProfiles-class
#' @exportClass EnrichmentProfiles
EnrichmentProfiles <- setClass("EnrichmentProfiles",
                               slots=c(method                   = "character",
                                       method.parameter         = "list",
                                       cluster.size             = "numeric",
                                       cluster.number           = "numeric",
                                       class.number             = "numeric",
                                       classes                  = "data.frame"),
                               validity = function(object){
                                    
                                    #TODO complete this
                                    return(TRUE)
                               }
)

#' @title Results object's name
#'
#' @description Named a Results object.
#'
#' @param x a Results object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","Results",
        definition = function(x){return("Results")}
)


#' @title SPADEResults object name
#'
#' @description Named a SPADEResults object.
#'
#' @param x a SPADEResults object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","SPADEResults",
        definition = function(x){return("SPADEResults")}
)


#' @title AC object name
#'
#' @description Named a AC object.
#'
#' @param x a AC object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","AC",
        definition = function(x){return("AC")}
)

#' @title DEC object name
#'
#' @description Named a DEC object.
#'
#' @param x a DEC object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","DEC",
        definition = function(x){return("DEC")}
)

#' @title CC object name
#'
#' @description Named a CC object.
#'
#' @param x a CC object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","CC",
        definition = function(x){return("CC")}
)


#' @title PhenoProfiles object name
#'
#' @description Named a PhenoProfiles object.
#'
#' @param x a PhenoProfiles object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","PhenoProfiles",
        definition = function(x){return("PhenoProfiles")}
)

#' @title EnrichmentProfiles object name
#'
#' @description Named a EnrichmentProfiles object.
#'
#' @param x a EnrichmentProfiles object
#' 
#' @return none
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","EnrichmentProfiles",
        definition = function(x){return("EnrichmentProfiles")}
)
