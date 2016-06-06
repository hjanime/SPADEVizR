#' @title Computation of the fraction of cluster with low number of cells
#' 
#' @description 
#' Computes the fraction of clusters having a number of cells less than a specific threshold
#' 
#' @details 
#' xxx
#' 
#' @param SPADEResults a SPADEResults object
#' @param num a numeric value specifying the cell threshold
#'  
#' @return a numeric value of the fraction of clusters having less cells that the specific threshold
#' 
computeNumberofClustersHavingLessThatThreshold <- function(SPADEResults, th = 5){

    data        <- SPADEResults@marker.expressions

    return(means)
    
}

#' @title Computation of the fraction of cluster with multiple population
#' 
#' @description 
#' Computes the fraction of clusters having a number of cells less than a specific threshold
#' 
#' @details 
#' xxx
#' 
#' @param SPADEResults a SPADEResults object
#' @param class a numeric value specifying the number of classes to found
#'  
#' @return a numeric value of the fraction of clusters having less cells that the specific threshold
computeNumberofClustersHavingMultiPop <- function(SPADEResults, class = 5){

    data        <- SPADEResults@marker.expressions

    return(means)
    
}