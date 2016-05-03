setOldClass("igraph") # give access to igraph class

#' @title Results class definition
#' 
#' @description 
#' The Results object is a S4 object containing cell clustering results obtained from various automatic gating algorithms. 
#' 
#' 
#' This object mainly stores the count matrix (i.e. the number of cells associated with each cluster of each sample) and the cell cluster phenotypes (i.e. the marker median expressions for each cluster). 
#' It is to note that the Results object is a super class of the SPADEResult object.
#' 
#' @details 
#' The 'cells.count' dataframe stores the number of cells associated with each cluster of each sample. This dataframe has in row the clusters and in column the samples.
#' 
#' The 'marker.expressions' dataframe stores the marker median expressions for each cluster. This dataframe has in the first the sample names, in the second column the cluster names, and the maker median expressions in the others columns.
#' 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'importX()' function.
#' 
#' @slot cells.count a dataframe containing the number of cells for each cluster of each sample
#' @slot marker.expressions a numerical dataframe containing marker median expressions for each cluster of each sample
#' @slot sample.names a character vector containing the sample names
#' @slot marker.names a character vector containing the markers names
#' @slot cluster.number a numeric specifying the number of cell clusters
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
#' @description 
#' The 'SPADEResults' object is a S4 object containing cell clustering results obtained from SPADE.
#' 
#' 
#' This object inherits from the 'Result' object and stores the count matrix (i.e. the number of cells associated with each cluster of each sample) and the cell cluster phenotypes (i.e. the marker median expressions for each cluster). 
#' In addition to the 'Result' object, the 'SPADEResults' object contains information about SPADE clustering results, such as the SPADE tree, the clustering makers and the FCS files.
#' 
#' @details 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'importSPADEResult()' function. 
#' 
#' @slot use.raw.medians a logical specifying if the marker expressions correspond to the raw or transformed data
#' @slot dictionary a two column data.frame providing the correspondence between the original marker names (first column) and the real marker names (second column)
#' @slot marker.clustering a logical vector specifying marker that have been used during the clustering precedure
#' @slot fcs.files a character vector containing the absolute path of the original FCS files
#' @slot quantiles a numeric data.frame containing the quantiles for each each markers of cluster
#' @slot graph a igraph object containing the SPADE tree
#' @slot graph.layout a numeric matrix containing the layout of the SPADE tree
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
                             if(ncol(object@cells.count) != length(object@sample.names)){
                                 message(paste0("Object SPADEResults, Error : number of samples (",
                                                 length(object@sample.names),
                                                 ") is inconsistent with cells.count matrix size (",
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

#' @title Abundant Clusters (AC class) definition
#' 
#' @description 
#' The 'AC' object is a S4 object containing the information related to the abundant clusters in a given biological condition. 
#' Moreover this object contains all parameters used in the statistical analysis.  
#' 
#' @details 
#' A cluster is considered as a significant abundant cluster if its associated p-value and mean are below the specific thresholds 'th.pvalue' and 'th.mean'.  
#' 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyAC()' function. 
#' 
#' @slot sample.names a character vector containing the samples used to compute the abundant clusters
#' @slot cluster.size a numeric vector containing the number of cells ( -- sum of all samples -- ) for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the abundant clusters
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot th.mean a numeric value specifying the mean threshold
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot result a data.frame containing for each cluster (first column): the mean (second column) and the standard deviation (third column) of the biological condition, the associated p-value (fourth column) and a logical (fifth column) specifying if the cluster is significantly abundant.
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
                       th.mean         = "numeric",
                       th.pvalue       = "numeric",   
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

#' @title Differentially Enriched Clusters (DEC class) definition
#' 
#' @description 
#' The 'DEC' object is a S4 object containing the information related to the differentially enriched clusters between two given biological conditions. 
#' Moreover this object contains all parameters used in the statistical analysis.  
#' 
#' @details 
#' A cluster is considered as a differentially enriched cluster if its associated p-value and fold-change are below the specific thresholds 'th.pvalue' and 'th.fc'.  
#' 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyDEC()' function. 
#' 
#' @slot sample.cond1 a character specifying the names of the samples of the first biological condition
#' @slot sample.cond2 a character specifying the names of the samples of the second biological condition
#' @slot cluster.size a numeric vector containing number of cells ( -- sum of all samples -- ) for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the DEC
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot method.paired a logical indicating if the statistical test have been performed in a paired manner
#' @slot th.fc a numeric value specifying the fold-change threshold
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot result a data.frame containing for each cluster (first column): the fold-change (second column) and the standard deviation (third column) for the first biological condition, the fold-change (fourth column) and the standard deviation (fifth column) for the second biological condition, the associated p-value (sixth column) and a logical (seventh column) specifying if the cluster is significantly differentially enriched.
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
                th.fc           = "numeric",
                th.pvalue       = "numeric",
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


#' @title Correlated Clusters (CC class)  definition
#' 
#' @description 
#' The 'CC' object is a S4 object containing coefficient of correlation associated between each cluster and a phenotypic variable.
#' Moreover this object contains all parameters used in the statistical analysis.
#' 
#' @details
#' A cluster is considered as a significant correlated cluster if its associated p-value and correlation threshold are below the specific thresholds 'th.pvalue' and 'th.correlation'.  
#' 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyCC()' function. 
#' 
#' @slot sample.names a character vector containing the samples used to compute correlated clusters
#' @slot variable a numeric vector containing the expression values of the associated variable
#' @slot cluster.size a numeric vector containing number of cells ( -- sum of all samples -- ) for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the CC
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot th.correlation a numeric value specifying the correlation threshold (R)
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot result a data.frame containing for each cluster (first column): the coefficiant of correlation R (second column) , the associated p-value (third column) and a logical (fourth column) specifying if the cluster is significantly correlated.
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
                       th.correlation  = "numeric",
                       th.pvalue       = "numeric",
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


#' @title Classification of clusters based on their phenotype profiles (PhenoProfiles class) definition
#' 
#' @description 
#' The 'PhenoProfiles' is a S4 object containing the information related to the cluster classification based on theirs marker expressions.
#' 
#' 
#' This object contains all information about the classification method and parameters used.  
#' 	
#' @details 
#' Five methods are available to classify cellular clusters: 'hierarchical_k', 'hierarchical_h', 'kmeans', 'eigencell' and 'clique'. Each method can paramaterized using the 'method.parameter' parameter.
#'  
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'classifyPhenoProfiles()' function. 
#' 
#' @slot cluster.number a numeric providing the number of cluster
#' @slot method a character specifying the method used to classify cluster
#' @slot method.parameter a named list of parameters used by the classification method
#' @slot cluster.size a numeric vector containing the number of cells associated with each cluster (-- sum of all samples --)

#' @slot class.number a numeric value specifying the number of clusters
#' @slot classes a two column dataframe with the cluster in first colunm and corresponding classe in the second colunm
#' 
#' @name PhenoProfiles-class
#' @rdname PhenoProfiles-class
#' @exportClass PhenoProfiles
PhenoProfiles <- setClass("PhenoProfiles",
                          slots=c(method           = "character",
                                  method.parameter = "numeric",
                                  cluster.size     = "numeric",
                                  cluster.number   = "numeric",
                                  class.number     = "numeric",
                                  classes          = "data.frame"),
                          validity = function(object){#TODO complete this
                                
                             return(TRUE)
                         }
)

#' @title Classification of clusters based on their enrichment profiles (EnrichmentProfiles class) definition
#' 
#' @description 
#' The 'EnrichmentProfiles' is a S4 object containing the information related to the cluster classification based on theirs enrichment profiles.
#' 
#' 
#' This object contains all information about the classification method and parameters used.   
#' 
#' @details
#' Five methods are available to classify cellular clusters: 'hierarchical_k', 'hierarchical_h', 'kmeans', 'eigencell' and 'clique'. Each method can paramaterized using the 'method.parameter' parameter.
#' 
#' The 'print()' and 'show()' can be used to display a summury of this object. Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'classifyEnrichmentProfiles()' function. 
#' 
#' @slot method a character specifying the method used to classify cluster
#' @slot method.parameter a numeric parameters associated with the chosen method
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
                                       method.parameter         = "numeric",
                                       cluster.size             = "numeric",
                                       cluster.number           = "numeric",
                                       class.number             = "numeric",
                                       classes                  = "data.frame"),
                               validity = function(object){
                                    
                                    #TODO complete this
                                    return(TRUE)
                               }
)

#' @title Definition of class names 
#'
#' @description 
#' Provides the name of each SPADEVizR object
#'
#' @param x a SPADEVizR object
#' 
#' @return a character providing the name of the object
#'  
#' @name names
#' @rdname names-methods
NULL

#' @rdname names-methods
#' @export
setMethod("names","Results",
        definition = function(x){return("Results")}
)

#' @rdname names-methods
#' @export
setMethod("names","SPADEResults",
        definition = function(x){return("SPADEResults")}
)

#' @rdname names-methods
#' @export
setMethod("names","AC",
        definition = function(x){return("AC")}
)

#' @rdname names-methods
#' @export
setMethod("names","DEC",
        definition = function(x){return("DEC")}
)

#' @rdname names-methods
#' @export
setMethod("names","CC",
        definition = function(x){return("CC")}
)

#' @rdname names-methods
#' @export
setMethod("names","PhenoProfiles",
        definition = function(x){return("PhenoProfiles")}
)

#' @rdname names-methods
#' @export
setMethod("names","EnrichmentProfiles",
        definition = function(x){return("EnrichmentProfiles")}
)
