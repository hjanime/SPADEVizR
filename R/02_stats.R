#' @title Identification of the Abundant Clusters
#' 
#' @description 
#' This function is used to identify the abundant clusters. That is to say clusters that have cell abundance statistically greater than a specific threshold.
#' 
#' @param Results a 'Results' or 'SPADEResults' object
#' @param condition a named vector providing the correspondence between a sample name (in row names) and the logical value TRUE to test abundance for this sample or FALSE otherwise
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character specifying the statistical method used to identify the abundant clusters. The parameter can take the values "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method) 
#' @param th.pvalue a numeric specifying the p-value threshold
#' @param th.mean a numeric specifying the abundance mean threshold
#' 
#' @return a S4 object of class 'AC'
#' 
#' @export
identifyAC <- function(Results,
                      condition,
                      use.percentages = TRUE,
                      method          = "t.test",
                      method.adjust   = NULL,
                      th.pvalue       = 0.05,
                      th.mean         = 0){
    
    message("[START] - computing ACs")
    
    cells.count   <- Results@cells.count[,names(condition[ condition == TRUE ]), drop = FALSE]
    data          <- cells.count
    if(use.percentages){
        data   <- prop.table(as.matrix(data),2)
        data   <- data * 100
    }else{
        data   <- data
    }
            
    message("Sampled used :")
    message(paste0(colnames(data),"\n"))
    
    pv <- apply(data, 1, function(x){
                                          return(do.call(method,
                                                         args = list(x           = x,
                                                                     alternative = "greater",
                                                                     mu          = th.mean
                                                         ))$p.value
                                                 )
                                    })
    
    if(!is.null(method.adjust)){
        pv <- p.adjust(pv, method = method.adjust)
    }
    
    result <- data.frame(cluster = rownames(Results@cells.count),
                         mean    = apply(data,1,mean),
                         sd      = apply(data,1,sd),
                         pvalue  = pv)
   
   result$significance <- ifelse(result$pvalue < th.pvalue, TRUE , FALSE)
   result$significance <- ifelse(abs(result$mean) > th.mean, result$significance , FALSE) 

    AC <- new ("AC",
               sample.names    = colnames(data),
               cluster.size    = apply(cells.count,1,sum),
               use.percentages = use.percentages,
               method          = method,
               method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust),#TODO think about another way
               th.mean         = th.mean,
               th.pvalue       = th.pvalue,
               result          = result)
    
    message("[END] - computing ACs")
	
    return(AC)
}

#' @title Identification of the Differentially Enriched Clusters
#' 
#' @description
#' This function is used to identify differentially enriched clusters. 
#' That is to say clusters that are differentially abundant between two biologicals conditions.
#' 
#' @param Results a 'Results' or 'SPADEResults' object
#' @param conditions a named vector providing the correspondence between a sample name (in row names) and the condition of this sample : NA to exclude a sample from tests, 1 or 2 to attribute this sample, respectively to the first or second condition
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character specifying the name of the statistical test to use "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method) 
#' @param method.paired a logical indicating if the statistical test must be performed in a paired manner
#' @param th.pvalue a numeric specifying the p-value threshold
#' @param th.fc a numeric specifying the fold-change threshold
#'
#' @return a S4 object of class 'DEC'
#' 
#' @export
identifyDEC <- function(Results,
                       conditions,
                       use.percentages = TRUE,
                       method          = "t.test",
                       method.adjust   = NULL,
                       method.paired   = FALSE,
                       th.pvalue       = 0.05,
                       th.fc           = 1){
    
    message("[START] - computing DECs\n")
    
    cells.count   <- Results@cells.count[, colnames(Results@cells.count) != "cluster"]
    data          <- cells.count
    
    data.cond1   <- data[,names(conditions[(!is.na(conditions) & conditions == 1)]), drop = FALSE]
    data.cond2   <- data[,names(conditions[(!is.na(conditions) & conditions == 2)]), drop = FALSE]
    
    data         <- cbind(data.cond1,data.cond2) 
    
    if(use.percentages){
        data   <- prop.table(as.matrix(data),2)
        data   <- data * 100
    }else{
        data   <- data
    }
        
    message("cond1:")
    message(paste0(colnames(data.cond1),"\n"))
    message("cond2:")
    message(paste0(colnames(data.cond2),"\n"))
    
    s1 <- ncol(data.cond1)
    
    pv <- apply(data, 1, function(x){
                return(do.call(method, args   = list(x = x[1:s1],
                                       y      = x[-(1:s1)],
                                       paired = method.paired)
                        )$p.value
                )
            }
    )
    
    if(!is.null(method.adjust)){
        pv <- p.adjust(pv, method = method.adjust)
    }
    
    fc <- apply(data, 1, function(x){
                fc <- mean(x[1:s1])/mean(x[-(1:s1)])
                if(!is.na(fc) && fc < 1){
                    fc <- (-1/fc)
                }
                return(fc)
            }
    )

    result <- data.frame(cluster     = rownames(Results@cells.count),
                         mean.cond1  = apply(data.cond1,1,mean),
                         sd.cond1    = apply(data.cond1,1,sd),
                         mean.cond2  = apply(data.cond2,1,mean),
                         sd.cond2    = apply(data.cond2,1,sd),
                         fold.change = fc,
                         pvalue      = pv)
             
    result$significance <- ifelse(result$pvalue < th.pvalue, TRUE , FALSE)
    result$significance <- ifelse(abs(result$fold.change) > th.fc, result$significance , FALSE)       
        
    thresholds <- c(pvalue = th.pvalue,fc = th.fc)

    DEC <- new ("DEC",
                sample.cond1    = colnames(data.cond1),
                sample.cond2    = colnames(data.cond2),
                cluster.size    = apply(cells.count,1,sum),
                use.percentages = use.percentages,
                method          = method,
                method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust), #TODO think about another way
                method.paired   = method.paired,
                th.fc           = th.fc,
                th.pvalue       = th.pvalue,
                result          = result)
    
    message("[END] - computing DECs")
    
    return(DEC)
}

#' @title Identification of the correlation of SPADE cluster with a phenotype
#' 
#' @description 
#' This function is used to identify correlated clusters. 
#' That is to say clusters that correlate with a phenotypic variable.
#' 
#' @param Results a 'Results' or 'SPADEResults' object
#' @param variable a numerical named vector providing the correspondence between a sample name (in rownames) and the specific phenotype or NA to ignore a sample
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character indicating the correlation method to use : "pearson", "spearman"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method) 
#' @param th.pvalue a numeric specifying the p-value threshold
#' @param th.correlation a numeric specifying the absolute value of the correlation coefficient threshold
#' 
#' @return a S4 object of class 'CC'
#'
#' @export
identifyCC <- function(Results,
                      variable,
                      use.percentages = TRUE,
                      method          = "pearson",
                      method.adjust   = NULL,
                      th.pvalue       = 0.05,
                      th.correlation  = 0.75){
                  
    message("[START] - computing CCs")

    cells.count <- Results@cells.count
    variable    <- na.omit(variable) 
    cells.count <- cells.count[, names(variable), drop = FALSE]
    data        <- cells.count
    
    if(use.percentages){
        data   <- prop.table(as.matrix(data), 2)
        data   <- data * 100
    }else{
        data   <- data
    }
    
    n            <- nrow(data)
    cor.estimate <- vector(mode = "numeric", length = n)
    cor.pvalue   <- vector(mode = "numeric", length = n)
    
    for(i in 1:n){
        cor             <- cor.test(variable,data[i,1:ncol(data)], method = method)
        cor.estimate[i] <- cor$estimate
        cor.pvalue[i]   <- cor$p.value
    }
    
    if(!is.null(method.adjust)){
        cor.pvalue <- p.adjust(cor.pvalue, method = method.adjust)
    }
    
    result <- data.frame(cluster     = rownames(Results@cells.count),
                         correlation = cor.estimate,
                         pvalue      = cor.pvalue)
        
    result$significance <- ifelse(result$pvalue < th.pvalue, TRUE, FALSE)
    result$significance <- ifelse(abs(result$correlation) > th.correlation, result$significance , FALSE)       

    CC <- new("CC",
              sample.names    = colnames(data),
              variable        = variable,
              cluster.size    = apply(cells.count,1,sum),
              use.percentages = use.percentages,
              method          = method,
              method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust),#TODO think about another way
              th.correlation  = th.correlation,
              th.pvalue       = th.pvalue,
              result          = result)
              
    message("[END] - computing CCs")
	
    return(CC)
}

#' @title Classification of clusters based on their phenotype profiles
#' 
#' @description 
#' Classifies clusters based on their phenotype profiles.
#' 
#' @details 
#' The classification is done on cell abundances of each clusters and could be performed using 5 methods:
#' \itemize{
#' \item "hierarchical_k" 
#' This method first compute the pearson correlation matrix and then use this matrix to performs a hierarchical classification. 
#' The hierarchical classification is cutted in order to return the desired number of classes. 
#' This number of classes must be provided as a numeric integer using the 'method.parameter' parameter.
#' It is to note that negative correlations are considered as uncorrelated
#' \item "hierarchical_h" (default method)
#' This method works in the same way than 'hierarchical_k' but the height where the hierarchical tree is specified. 
#' This heigth is a correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' \item "eigencell" 
#' This method compute the performs a eigen vector decomposition and then calculate the correlations between cluster expressions and these vectors.
#' Clusters which correlate above a specific threshold with the same eigen vector are classified together.
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.8) provided using the 'method.parameter' parameter.
#' \item "clique" 
#' This method first compute the Pearson correlation matrix and then use this matrix to generate an undirected graph.
#' In this graph, an edge is drawn between two nodes if the correlation coefficient in the adjacency matrix is above a specific threshold. 
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' After building the graph, the method looking for the largest cliques wich are considered as classes of nodes. Cliques correspond to subgraph in which every two distinct vertices are adjacent.
#' }
#' 
#' @param Results a Results or SPADEResults object
#' @param method a character specifying the clustering method among one of those : "hierarchical_h", "hierarchical_k","k-means","eigencell","clique"
#' @param method.parameter a numeric specifying the numeric value required by the selected method 
#' 
#' @return a S4 object of class 'PhenoProfiles'
#'
#' @export
classifyPhenoProfiles <- function (Results,
                                   method           = "hierarchical_h",
                                   method.parameter = NULL){
    default.eigencell.correlation.th    <- 0.8
    default.clique.correlation.th       <- 0.7
    default.hierarchical.correlation.th <- 0.7                           
    
    message("[START] - computing classifyPhenoProfiles")
    
    if (!is.element(method,c("hierarchical_h", "hierarchical_k", "k-means", "eigencell", "clique"))){
        stop("Error : In classifyPhenoProfiles, method must be one of those : 'hierarchical_h','hierarchical_k','k-means','eigencell','clique'")
    }
    
    table <- computePhenoTable(Results)
    
    #print(table)
    
    table.wide <- reshape2::dcast(table, cluster~marker)
    table.wide <- table.wide[, colnames(table.wide) != "cluster"]
    rownames(table.wide) <- rownames(Results@cells.count)   

    table.wide <- na.omit(table.wide) 

    #print(table.wide)
    
    switch(method,
            "hierarchical_h" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.hierarchical.correlation.th
                }
                classes <- computeHierarchicalClustering(table.wide, class.number = NULL, hierarchical.correlation.th = method.parameter)
            },
            "hierarchical_k" = {
                if (is.null(method.parameter)){
                    stop("Error, class.number can't be null with 'hierarchical_k' method")
                }
                classes <- computeHierarchicalClustering(table.wide, class.number = method.parameter, hierarchical.correlation.th = NULL)
            },
            "k-means" = {
                if (is.null(method.parameter)){
                    stop("Error, class.number can't be null with 'k-means' method")
                }
                classes <- computeKmeans(table.wide, method.parameter)
            },
            "eigencell" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.eigencell.correlation.th
                }
                classes <- computeEigenCellClusters(table.wide, method.parameter)
            },
            "clique" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.clique.correlation.th
                }
                classes <- computeClique(table.wide, method.parameter)
            })
    
    classes$class <- as.numeric(classes$class)
    classes       <- classes[ order(classes$class), ]
    
    cluster.size           <- apply(Results@cells.count, 1, sum)
    names(cluster.size) <- rownames(Results@cells.count)
    
    pheno <- new ("PhenoProfiles",
                  method                   = method,
                  method.parameter         = method.parameter,
                  cluster.size             = cluster.size,
                  cluster.number           = Results@cluster.number,
                  class.number             = length(unique(classes[,2])),
                  classes                  = classes)
    
    message("[END] - computing classifyPhenoProfiles")
    return(pheno)    
}

#' @title Classification of clusters based on their enrichment profiles
#' 
#' @description 
#' Classifies clusters based on their enrichment profiles.
#' 
#' @details 
#' The classification is done on cell abundances of each clusters and could be performed using 5 methods:
#' \itemize{
#' \item "hierarchical_k" 
#' This method first compute the pearson correlation matrix and then use this matrix to performs a hierarchical classification. 
#' The hierarchical classification is cutted in order to return the desired number of classes. 
#' This number of classes must be provided as a numeric integer using the 'method.parameter' parameter.
#' It is to note that negative correlations are considered as uncorrelated
#' \item "hierarchical_h" (default method)
#' This method works in the same way than 'hierarchical_k' but the height where the hierarchical tree is specified. 
#' This heigth is a correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' \item "eigencell" 
#' This method compute the performs a eigen vector decomposition and then calculate the correlations between cluster expressions and these vectors.
#' Clusters which correlate above a specific threshold with the same eigen vector are classified together.
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.8) provided using the 'method.parameter' parameter.
#' \item "clique" 
#' This method first compute the Pearson correlation matrix and then use this matrix to generate an undirected graph.
#' In this graph, an edge is drawn between two nodes if the correlation coefficient in the adjacency matrix is above a specific threshold. 
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' After building the graph, the method looking for the largest cliques wich are considered as classes of nodes. Cliques correspond to subgraph in which every two distinct vertices are adjacent.
#' }
#' 
#' @param Results a Results or SPADEResults object
#' @param method a character specifying the clustering method among one of those : 'hierarchical_h','hierarchical_k','k-means','eigencell','clique' (hierarchical_h by default)
#' @param method.parameter a numeric specifying the numeric value required by the selected method
#' 
#' @return a S4 object of class 'EnrichmentProfiles'
#'
#' @export
classifyEnrichmentProfiles <- function(Results,
                                       method           = "hierarchical_h",
                                       method.parameter = NULL){ # think about select sample ?
 
    default.eigencell.correlation.th    <- 0.8
    default.clique.correlation.th       <- 0.7
    default.hierarchical.correlation.th <- 0.7                                 
     
    message("[START] - computing classifyEnrichmentProfiles")
    
    if (!is.element(method,c("hierarchical_h", "hierarchical_k", "k-means", "eigencell", "clique"))){
        stop("Error : In classifyEnrichmentProfiles, method must be one of those : 'hierarchical_h','hierarchical_k','k-means','eigencell','clique'")
    }
    
    data     <- Results@cells.count
    
    switch(method,
            "hierarchical_h" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.hierarchical.correlation.th
                }
                classes <- computeHierarchicalClustering(data, class.number = NULL, hierarchical.correlation.th = method.parameter)
            },
            "hierarchical_k" = {
                if (is.null(method.parameter)){
                    stop("Error, class.number can't be null with 'hierarchical_k' method")
                }
                classes <- computeHierarchicalClustering(data, class.number = method.parameter)
            },
            "k-means" = {
                if (is.null(method.parameter)){
                    stop("Error, class.number can't be null with 'k-means' method")
                }
                classes <- computeKmeans(data, method.parameter)
            },
            "eigencell" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.eigencell.correlation.th
                }
                classes <- computeEigenCellClusters(data, method.parameter)
            },
            "clique" = {
                if (is.null(method.parameter)){
                    method.parameter <- default.clique.correlation.th
                }
                classes <- computeClique(data, method.parameter)
            })
   classes$class <- as.numeric(classes$class)
   classes       <- classes[ order(classes$class), ]
   
   cluster.size        <- apply(Results@cells.count,1,sum)
   names(cluster.size) <- rownames(Results@cells.count)
   
   enrich <- new ("EnrichmentProfiles",
                  method                   = method,
                  method.parameter         = method.parameter,
                  cluster.size             = cluster.size,
                  cluster.number           = Results@cluster.number,
                  class.number             = length(unique(classes[,2])),
                  classes                  = classes)
                 
    message("[END] - computing classifyEnrichmentProfiles")
    return(enrich)  
}

#' @title Internal - Hierarchical classification
#' 
#' @description 
#' This function is used internally to classify clusters enrichment profiles or phenotype profiles using a hierarchical algorithm. 
#' 
#' @details 
#' This function compute the Pearson correlation matrix associated to the provided matrix. 
#' It is to note that negative correlations are considered as uncorrelated.
#' This correlation matrix is used to performs a hierarchical classification.
#' If 'class.number' is NULL, classification will be determined base on the cut height correlation threshold.
#' 
#' @param data a matrix with all clusters in rownames
#' @param class.number a numeric specifying the number of classes
#' @param hierarchical.correlation.th a numeric value specifying the cut height
#' 
#' @return a dataframe containing for each cluster, its name and class
#' 
computeHierarchicalClustering <- function (data,
                                           class.number                = NULL,
                                           hierarchical.correlation.th = 0.8){

    cor.data <- cor(t(data))
    #print(head(cor.data))
    
    cor.data[cor.data < 0] <- 0
    cor.data               <- 1 - cor.data
    
    dist.cor.data <- as.dist(cor.data)
    tree          <- hclust(dist.cor.data)
    
    if (!is.null(class.number)){
        res <- cutree(tree, k = class.number)
    }else{
        res <- cutree(tree, h = hierarchical.correlation.th)
    }
    
    res <- data.frame (cluster = as.character(names(res)), class = res)
    
    return(res)
}


#' @title Internal - Kmeans classification
#' 
#' @description 
#' This function is used internally to classify clusters enrichment profilies or phenotype profiles using a k-means algorithm. 
#' 
#' @details 
#' xxx
#' 
#' @param data a numeric matrix with cluster names in rownames
#' @param k a numeric specifying the desired number of classes
#' 
#' @return a dataframe containing for each cluster, its name and class
#' 
computeKmeans <- function(data,
                          k = NULL){
    kmeans <- kmeans(data, centers = k)    
    result <- data.frame(cluster = as.character(rownames(data)), class = as.numeric(kmeans$cluster))
        
    return(result)
}


#' @title Internal - Eigen vector classification
#' 
#' @description 
#' This function is used internally to classify clusters enrichment profiles or phenotype profiles using eigen vector decomposition. 
#' 
#' @details 
#' xxx
#' 
#' @param data a numeric matrix with all clusters in rownames
#' @param eigencell.correlation.th a numeric value indicating the correlation coefficient threshold
#' 
#' @return a dataframe containing for each cluster, its name and class
#' 
computeEigenCellClusters <- function(data, 
                                     eigencell.correlation.th = 0.80){
    
    svd     <- svd(data)
    eigenCC <- t(svd$v)
    res     <- c()
    for(i in c(1:nrow(eigenCC))){
        for(j in c(1:nrow(data))){
            cor <- cor(as.numeric(eigenCC[i,]),as.numeric(data[j,]),method = "spearman")
            if(cor > eigencell.correlation.th){
                res <- rbind(res,cbind(as.character(rownames(data[j,])),i))
            }
        }
    }

    res <- as.data.frame(res)
    colnames(res) <- c("cluster","class")
    
    classes.uniq <- unique(res[,"class"])
    classes <- data.frame(class = classes.uniq, class.ID =  1:length(classes.uniq))
    
    for (i in 1:nrow(res)){
        res[i, "class"] <- classes[classes$class == res[i, "class"], "class.ID"]
    }
    
    unclassified <- setdiff(rownames(data), res$cluster)
    print(unclassified)
    res <- rbind(res, data.frame(cluster = unclassified, class = rep(NA, length(unclassified))))

    return(res)
}

#' @title Internal - Clique percolation classification
#' 
#' @description 
#' This function is used internally to classify clusters enrichment profiles or phenotype profiles using a clique percolation algorithm. 
#' 
#' @details 
#' xxx
#' 
#' @param data a numeric matrix with all clusters in rownames
#' @param clique.correlation.th a numeric value indicating the correlation coefficient threshold
#' 
#' @return a dataframe containing for each cluster, its name and class
#' 
computeClique <- function(data,
                          clique.correlation.th = 0.7){
    
    res     <- c()
    for(i in 1:(nrow(data)-1)){
        for(j in (i+1):nrow(data)){
            cor <- cor(as.numeric(data[i,]), as.numeric(data[j,]), method = "pearson")
            if(cor>clique.correlation.th){
                res <- rbind(res, cbind(j,i,cor))
            }
        }
    }
    colnames(res) <- c("cluster", "cluster", "cor")
    res <- data.frame(res)
    
    graph <- igraph::graph.data.frame(res, directed=FALSE)

    lists <- igraph::largest.cliques(graph)
    
    res <- data.frame()
    
    for(i in 1:length(lists)){
        cluster <- as.character(rownames(data[names(lists[[i]]),]))
        res <- rbind(res, cbind(cluster, i))
    }

    colnames(res) <- c("cluster", "class")
    
    unclassified <- setdiff(rownames(data), res$cluster)
    res <- rbind(res, data.frame(cluster = unclassified, class = rep(NA, length(unclassified))))

    return(res)
    
}