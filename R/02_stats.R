#' @title Compute the Abundant Clusters
#' 
#' @description Abundant Clusters are clusters wich have an enrichement significatyvely different of zeor
#'  
#' @details xxx
#' 
#' @param Results a Results or SPADEResults object
#' @param condition a named vector providing the correspondence between a sample name (in rowname) and the logical value TRUE to test abondance for this sample or FALSE otherwise
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character containing the statistical method to use for the ACs detection. The parameter can take the values "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the pvalues should be corrected using the argument "method" available for function p.adjust, among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr"
#' @param th.pvalue a numeric specifying the pvalue threshold (0.05 by default)
#' @param th.mean a numeric specifying the mean threshold (0 by default)
#' 
#' @return a AC object
#' 
#' @export
computeAC <- function(Results,
                      condition,
                      use.percentages = TRUE,
                      method          = "t.test",
                      method.adjust   = NULL,
                      th.pvalue       = 0.05,
                      th.mean         = 0){
    
    message("[START] - computing ACs")
    
    data.filtered   <- Results@cells.count[, colnames(Results@cells.count) != "cluster"]

    data.filtered   <- data.filtered[,names(condition[ condition == TRUE ]), drop = FALSE]
    
    if(use.percentages){
        data   <- prop.table(as.matrix(data.filtered),2)
        data   <- data * 100
    }else{
        data   <- data.filtered
    }
            
    message("Sampled used :")
    message(paste0(colnames(data),"\n"))
    
    pv <- apply(data, 1, function(x){
                                          return(do.call(method,
                                                         args = list(x           = x,
                                                                     alternative = "greater",
                                                                     mu          = 0
                                                         ))$p.value
                                                 )
                                    })
    
    if(!is.null(method.adjust)){
        pv <- p.adjust(pv, method = method.adjust)
    }
    
    result <- data.frame(cluster = Results@cells.count[,"cluster"],
                         mean    = apply(data,1,mean),
                         sd      = apply(data,1,sd),
                         pvalue  = pv)
    
   result$significance <- ifelse(result$pvalue < th.pvalue, TRUE , FALSE)
   result$significance <- ifelse(abs(result$mean) > th.mean, result$significance , FALSE) 

    AC <- new ("AC",
               sample.names    = colnames(data),
               cluster.size    = apply(data.filtered,1,sum),
               use.percentages = use.percentages,
               method          = method,
               method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust),#TODO think about another way
               th.mean         = th.mean,
               th.pvalue       = th.pvalue,
               result          = result)
       
    message("[END] - computing ACs")
	
    return(AC)
}

#' @title Compute the Differentially Enriched Clusters
#' 
#' @description Differentially Enriched Clusters are clusters for which the means of cell number are significantly different between two conditions.
#'
#' @details xxx
#' 
#' @param result a Results or SPADEResults object
#' @param conditions a named vector providing the correspondence between a sample name (in rownames) and the condition of this sample : NA to exclude a sample from tests, 1 or 2 to attribute this sample, respectively to the first or second condition
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character specifying the name of the statistical test to use "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the pvalues should be corrected using the argument "method" available for function p.adjust, among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr"
#' @param method.paired a logical indicating whether data measurement are paired (by default FALSE)
#' @param th.pvalue a numeric specifying the pvalue threshold (0.05 by default)
#' @param th.fc a numeric specifying the foldchange threshold (1 by default)
#'
#' @return a DEC object
#' 
#' @export
computeDEC <- function(Results,
                       conditions,
                       use.percentages = TRUE,
                       method          = "t.test",
                       method.adjust   = NULL,
                       method.paired   = FALSE,
                       th.pvalue       = 0.05,
                       th.fc           = 1){
    
    message("[START] - computing DECs\n")
    
    data.filtered   <- Results@cells.count[, colnames(Results@cells.count) != "cluster"]
    
    data.cond1   <- data.filtered[,names(conditions[(!is.na(conditions) & conditions == 1)]), drop = FALSE]
    data.cond2   <- data.filtered[,names(conditions[(!is.na(conditions) & conditions == 2)]), drop = FALSE]
    
    data.filtered         <- cbind(data.cond1,data.cond2) 
    
    if(use.percentages){
        data   <- prop.table(as.matrix(data.filtered),2)
        data   <- data * 100
    }else{
        data   <- data.filtered
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

    result <- data.frame(cluster     = Results@cells.count[,"cluster"],
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
                cluster.size    = apply(data.filtered,1,sum),
                use.percentages = use.percentages,
                method          = method,
                method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust), #TODO think about another way
                method.paired   = method.paired,
                th.pvalue       = th.pvalue,
                th.fc           = th.fc,
                result          = result)
    
    message("[END] - computing DECs")
    
    return(DEC)
}

#' @title Compute the correlation of SPADE cluster with a cynetics phenotype
#' 
#' @description Correlated Clusters are clusters with count or % correlated with a specific phenotipical variable (such as titer)
#' 
#' @details xxx
#' 
#' @param Results a Results or SPADEResults object
#' @param variable a numerical named vector providing the correspondence between a sample name (in rowname) and the specific phenotype or NA to ignore a sample
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character indicating the correlation method to use : "pearson", "spearman"
#' @param method.adjust a character specifying if the pvalues should be corrected using the argument "method" available for function p.adjust, among : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr"
#' @param th.pvalue a numeric specifying the pvalue threshold (0.05 by default)
#' @param th.correlation a numeric specifying the r two sided threshold (0.5 by default) in [0,1]
#' 
#' @return a CC object
#'
#' @export
computeCC <- function(Results,
                      variable,
                      use.percentages = TRUE,
                      method     = "pearson",
                      method.adjust   = NULL,
                      th.pvalue       = 0.05,
                      th.correlation  = 0.5){
                  
    message("[START] - computing CCs")

    data.filtered <- Results@cells.count[, colnames(Results@cells.count) != "cluster"]
    variable      <- na.omit(variable) 
    data.filtered <- data.filtered[,names(variable), drop = FALSE]
    
    if(use.percentages){
        data   <- prop.table(as.matrix(data.filtered),2)
        data   <- data * 100
    }else{
        data   <- data.filtered
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
    
    result <- data.frame(cluster         = Results@cells.count[,"cluster"],
                         correlation     = cor.estimate,
                         pvalue          = cor.pvalue)
        
    result$significance <- ifelse(result$pvalue < th.pvalue, TRUE, FALSE)
    result$significance <- ifelse(abs(result$correlation) > th.correlation, result$significance , FALSE)       

    CC <- new("CC",
              sample.names    = colnames(data),
              variable        = variable,
              cluster.size    = apply(data.filtered,1,sum),
              use.percentages = use.percentages,
              method          = method,
              method.adjust   = ifelse(is.null(method.adjust),"none",method.adjust),#TODO think about another way
              th.pvalue       = th.pvalue,
              th.correlation  = th.correlation,
              result          = result)
              
    message("[END] - computing CCs")
	
    return(CC)
}

#' @title classifyPhenoProfiles xxx
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param Results a Results or SPADEResults object
#' @param method a character specifying the clustering method among one of those : 'hierarchical','kmeans','eigencell','clique'
#' @param class.number a numeric specifying the number of classes needed when the method parameter choosen is either 'hierarchical_k' or 'kmeans'
#' @param eigencell.correlation.th a numeric (ignored if method is not 'eigencell') specifying the correlation threshold (in [0,1], 0.8 by default)  in case of eigencell clustering
#' @param clique.correlation.th a numeric (ignored if method is not 'clique') specifying the correlation threshold (in [0,1], 0.7 by default) in case of clique clustering
#' @param hierarchical.correlation.th a numeric (ignored if method is not 'hierarchical_h') in [0,1]) specifying the threshold of correlation (in [0,1], 0.8 by default) use to cut the hirerchical tree
#' 
#' @return a PhenoProfiles object
#'
#' @export
classifyPhenoProfiles <- function (Results,
                                   method                      = "hierarchical_h",
                                   class.number                = NULL,
                                   eigencell.correlation.th    = 0.8,
                                   clique.correlation.th       = 0.7,
                                   hierarchical.correlation.th = 0.8){
    message("[START] - computing classifyPhenoProfiles")
    
    if (!is.element(method,c("hierarchical_h", "hierarchical_k", "kmeans", "eigencell", "clique"))){
        stop("Error : In classifyPhenoProfiles, method must be one of those : 'hierarchical_h','hierarchical_k','kmeans','eigencell','clique'")
    }
    
    table <- computePhenoTable(Results)
    
    print(table)
    
    table.wide <- reshape2::dcast(table, cluster~marker)
    table.wide <- table.wide[, colnames(table.wide) != "cluster"]
    rownames(table.wide) <- Results@cells.count[,"cluster"]   

    table.wide <- na.omit(table.wide) 

    method.parameter <- list()
    
    print(table.wide)
    
    switch(method,
            hierarchical_h = {
                classes <- computeHierarchicalClustering(data, class.number = NULL, hierarchical.correlation.th)
                method.parameter <- list(hierarchical.correlation.th = hierarchical.correlation.th)
            },
            hierarchical_k = {
                if (!is.null(class.number)){
                classes <- computeHierarchicalClustering(data, class.number = class.number)
                method.parameter <- list(class.number = class.number)
                }else{
                    Stop("Error, class.number can't be null with hierarchical_k")
                }
            },
            kmeans = {
                if (!is.null(class.number)){
                    classes <- computeKMeans(table.wide, k = class.number)
                    print(classes)
                    method.parameter <- list(class.number = class.number)
                }else{
                    Stop("Error, class.number can't be null with KMeans")
                }
            },
            eigencell = {
                classes <- computeEigenCellClusters(table.wide, eigencell.correlation.th)
                method.parameter <- list(eigencell.correlation.th = eigencell.correlation.th)
            },
            clique = {
                classes <- computeClique(table.wide, clique.correlation.th)
                method.parameter <- list(clique.correlation.th = clique.correlation.th)
            })
    
    cluster.count <- subset(Results@cells.count, cluster %in% classes$cluster)
    cluster.count <- cluster.count[, colnames(Results@cells.count) != "cluster"]
    cluster.size <- apply(cluster.count, 1, sum)
    
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

#' @title classifyEnrichmentProfiles xxx
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param Results a Results or SPADEResults object
#' @param method a character specifying the clustering method among one of those : 'hierarchical','kmeans','eigencell','clique'
#' @param class.number a numeric specifying the number of classes needed when the method parameter choosen is either 'hierarchical_k' or 'kmeans'
#' @param eigencell.correlation.th a numeric (ignored if method is not 'eigencell') specifying the correlation threshold (in [0,1], 0.8 by default)  in case of eigencell clustering
#' @param clique.correlation.th a numeric (ignored if method is not 'clique') specifying the correlation threshold (in [0,1], 0.7 by default) in case of clique clustering
#' @param hierarchical.correlation.th a numeric (ignored if method is not 'hierarchical_h') in [0,1]) specifying the threshold of correlation (in [0,1], 0.8 by default) use to cut the hirerchical tree
#' 
#' @return a EnrichmentProfiles object 
#'
#' @export
classifyEnrichmentProfiles <- function(Results,
                                       method                      = "hierarchical_h",
                                       class.number                = NULL,
                                       eigencell.correlation.th    = 0.8,
                                       clique.correlation.th       = 0.7,
                                       hierarchical.correlation.th = 0.8){ # think about select sample ?
    message("[START] - computing classifyEnrichmentProfiles")
    
    if (!is.element(method,c("hierarchical_h", "hierarchical_k", "kmeans", "eigencell", "clique"))){
        stop("Error : In classifyEnrichmentProfiles, method must be one of those : 'hierarchical_h','hierarchical_k','kmeans','eigencell','clique'")
    }
    
    data     <- Results@cells.count[, colnames(Results@cells.count) != "cluster"]
    rownames(data) <- Results@cells.count[,"cluster"] 
    
    method.parameter <- list()
    
    switch(method,
           hierarchical_h = {
               classes <- computeHierarchicalClustering(data, class.number = NULL, hierarchical.correlation.th)
               method.parameter <- list(hierarchical.correlation.th = hierarchical.correlation.th)
           },
           hierarchical_k = {
               if (!is.null(class.number)){
                   classes <- computeHierarchicalClustering(data, class.number = class.number)
                   method.parameter <- list(class.number = class.number)
               }else{
                   Stop("Error, class.number can't be null with hierarchical_k")
               }
           },
           kmeans = {
               classes <- computeKMeans(data, k = class.number)
               method.parameter <- list(class.number = class.number)
           },
           eigencell = {
               classes <- computeEigenCellClusters(data, eigencell.correlation.th)
               method.parameter <- list(eigencell.correlation.th = eigencell.correlation.th)
           },
           clique = {
               classes <- computeClique(data, clique.correlation.th)
               method.parameter <- list(clique.correlation.th = clique.correlation.th)
           })

   cluster.count <- subset(Results@cells.count, cluster %in% classes$cluster)
   cluster.count <- cluster.count[, colnames(Results@cells.count) != "cluster"]
   cluster.size <- apply(cluster.count,1,sum)
   
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

#' @title Internal - computeHierarchicalClustering
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param data a matrix with all clusters in rownames
#' @param class.number the number of classe to cluster, if class.number is NULL numer of classe will be determined base on hierarchical correlation threshold
#' @param hierarchical.correlation.th
#' 
#' @return a dataframe containing the cluster ID, the classe of this cluster and a clustering score.
#' 
computeHierarchicalClustering <- function (data,
                                           class.number                = NULL,
                                           hierarchical.correlation.th = 0.8){

    cor.data <- cor(t(data))
    print(head(cor.data))
    
    cor.data[cor.data < 0] <- 0
    cor.data               <- 1 - cor.data
    
    dist.cor.data <- as.dist(cor.data)
    tree          <- hclust(dist.cor.data)
    
    if (!is.null(class.number)){
        res <- cutree(tree, k = class.number)
    }else{
        res <- cutree(tree, h = hierarchical.correlation.th)
    }
    
    res <- data.frame (cluster = as.character(names(res)), classe = res)
    res <- res[order(res$cluster),]
    
    return(res)
}


#' @title Internal - computeKMeans
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param data a matrix with all clusters in rownames
#' @param k number of classes
#' 
#' @return a dataframe containing the cluster ID, the classe of this cluster and a clustering score.
#' 
computeKMeans <- function(data,
                          k){
    kmeans <- kmeans(data, centers = k)    
    result <- data.frame(cluster = as.character(rownames(data)), classe = as.numeric(kmeans$cluster), stringsAsFactors = FALSE )
        
    return(result)
}


#' @title Internal - computeEigenCellClusters
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param data a matrix with all clusters in rownames
#' @param eigencell.correlation.th
#' 
#' @return a dataframe containing the cluster ID, the classe of this cluster and a clustering score.
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
    colnames(res) <- c("cluster","classe")
    
    classes.uniq <- unique(res[,"classe"])
    classes <- data.frame(classe = classes.uniq, classe.ID =  1:length(classes.uniq))
    
    for (i in 1:nrow(res)){
        res[i,"classe"] <- classes[classes$classe == res[i,"classe"],"classe.ID"]
    }
    
    res <- res[order(res$cluster),]
    return(res)
}

#' @title Internal - computeClique
#' 
#' @description xxx
#' 
#' @details xxx
#' 
#' @param data a matrix with all clusters in rownames
#' @param eigencell.correlation.th
#' 
#' @return a dataframe containing the cluster ID, the classe of this cluster and a clustering score.
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

    colnames(res) <- c("cluster", "classe")
    res <- res[order(res$cluster),]

    return(res)
    
}