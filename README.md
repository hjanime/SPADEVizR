# SPADEVizR: an R package for the visualization and analysis of SPADE clustering results
Guillaume Gautreau, David Pejoski, Ludovic Platon, Brice Targat, Roger Le Grand, Anne-Sophie Beignon and Nicolas Tchitchek  



![](logoSPADEVizR.png)

# Table of Contents
1. [Package overview](#package_overview)
2. [Package installation](#package_installation)
3. [Importing automatic gating results](#loading_data)
	1. [Importing automatic gating results from SPADE](#loading_SPADE_data)
	2. [Importing automatic gating results from other algorithms](#loading_other_data)
4. [Object structures](#object_structures)
	1. [Results object](#object_structure_results)
	2. [SPADEResults object](#object_structure_SPADE_results)
	3. [Abundant clusters (AC)](#object_structure_AC)
	4. [Differentially enriched clusters (CC)](#object_structure_DEC)
	5. [Computation of correlated clusters (computeCC function)](#object_structure_CC)
	6. [Classification of clusters based on their phenotype profiles (PhenoProfiles object)](#object_structure_PhenoProfiles)
	7. [Classification of clusters based on their enrichment profiles (EnrichmentProfiles object)](#object_structure_EnrichmentProfiles)
5. [Statistical analyses](#stat_functions)
	1. [Computation of abundant clusters (computeAC function)](#stat_function_computeAC)
	2. [Computation of differentially enriched clusters (computeDEC function)](#stat_function_computeDEC)
	3. [computeCC (Correlated Clusters)](#stat_function_computecc)
	4. [Classification of phenotype profiles (classifyPhenoProfiles function)](#stat_function_classifyPhenoProfiles)
	5. [Classification of enrichment profiles (classifyEnrichmentProfiles function)](#stat_function_classifyEnrichmentProfiles)
6. [Vizualisation of statistical results](#viewer_functions)
	1. [Abundant Clusters Viewer (AC plot)](#abundant_clusters_viewer_function)
	2. [Volcano Viewer (DEC plot)](#volcano_viewer_function)
	3. [Correlated Clusters Viewer (CC plot)](#correlated_clusters_viewer_function)
	4. [Profile Viewer](#profile_viewer_function)
7. [Miscellaneous vizualisations](#viewer_functions)
	1. [Cluster Viewer](#cluster_viewer_function)
	2. [Pheno Viewer](#pheno_viewer_function)
	3. [Tree Viewer](#tree_viewer_function)
	3. [MDS Viewer (Multidimensional Scaling)](#MDS_viewer_function)
	4. [Distogram Viewer](#distogram_viewer_function)
	5. [Streamgraph Viewer](#streamgraph_viewer_function)
	6. [Boxplot Viewer](#boxplot_viewer_function)
8. [Generate report](#generate_report_function)
9. [Export](#generate_report_function)
10. [License](#license)
11. [References](#references)

# <a name="package_overview"/> 1. Package overview

--Flow and mass cytometry are experimental techniques used to characterize cells at a single-cell level. Both techniques use labelled antibodies to measure cell markers. Flow cytometry uses antibodies conjugated with fluorochromes to quantify stained cells with a laser detection system, and recently introduced mass cytometry ([CyTOF](https://www.fluidigm.com/products/cytof) [[1](http://www.ncbi.nlm.nih.gov/pubmed/21551058)]) uses antibodies conjugated with metals to quantify stained cells with a mass spectrometer. While flow cytometry can currently handle up to 18 fluorochromes, mass cytometry can currently handle up to 40 metals. Mass cytometry offers important perspectives as it can potentially evaluate more than 100 metals. Additionally, another recently introduced cytometry technique, called hyperspectral cytometry technique [[2](http://www.ncbi.nlm.nih.gov/pubmed/24271566)], combines ultrafast optical spectroscopy with flow cytometry and can also increase up to 40 the number of simultaneously usable fluorochromes by avoiding spectral overlaps. Such significant technological improvements make flow and mass cytometry suitable for massive datamining and bioinformatics developments, in order to better explore and analyze complex cell interaction systems.

Cytometry data can be manually analyzed using a gating strategy (hierarchical gating usually followed by Boolean gating) or using automatic gating/clustering algorithms. In both approaches, the aim is to identify cell populations having similar profiles, based on the expressions of selected markers. Hierarchical gating is done using biplot representations, where each axis represents the intensity of a marker of interest and where dots represent the cells. Through iterative steps, operators select cell subsets and explore their marker expression profiles. Boolean gating is done by quantifying cells present within all possible combinations of gates for a population of interest. Even if Boolean gating can be done in an automatic way, gates still have to be delineated for each cell marker. FlowJo FlowJo [[3](http://www.flowjo.com)] and CytoBank [[4](https://www.cytobank.org)] are among the most common software for manual gating. Both hierarchical and Boolean gating can be very fastidious and biased. The difficulty mainly comes from the high dimensionality of data that is to say from the high number of markers and cell populations. This task becomes even more complicated as the number of markers increases or as the studies become more complex. On the other hand, [SPADE](http://cytospade.org/) [[5](http://www.ncbi.nlm.nih.gov/pubmed/21964415)] is an automatic gating methods using an efficient algorithm to computationally identify cell populations having similar profiles and provides then less biased results.--

SPADE offers new opportunities to deal with high dimensional cytometric data in order to identify cell populations but does not provide strong visualization features and statistical approachs. SPADEVizR aims to improve the biological analysis of SPADE results with new visualization approaches. We extended the original SPADE visualization outputs with visualization techniques such as parallel coordinates, multidimensional scaling, volcano plots or streamgraph representations as ggplot2[6]. Moreover SPADEVizR allows to identify automatically abundont cluster that is to say for a given condition, if some clusters are statistically abundant or not. Differentially Enriched Clusters, that is to say  Identify differentially enriched cell populations (clusters) between 2 biological conditions and cluster with phenot... Futhermore to automatically classify cell populations that have a similar phenotypes or a similar enrichments profiles.

Beyond representations, SPADEVizR allows to automatically: 

* Identify differentially enriched cell populations (clusters) between 2 biological conditions. 
* Determine for a given condition, if some clusters are statistically abundant or not. 
* Find correlations between cluster kinetics and external phenotypical data.
* Classify clusters by their similar phenotypes or their similar enrichments profiles.

# <a name="package_installation"/> 2. Package installation
The `ggplot2`, `ggrepel`, `grid`, `gridExtra`, `igraph`, `MASS`, `reshape2`, `gtools`, `ggdendro` R packages as well as the `flowCore` [7] Bioconductor packages are required for running SPADEVizR. These packages can be installed using the following commands:

```r
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('ggtools')
install.packages('ggdendro')
install.packages('grid')
install.packages('gridExtra')
install.packages('igraph')
install.packages('MASS')
install.packages('scales')
install.packages('reshape2')
#install.packages('ggtimeseries')

source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates=TRUE)
biocLite("flowCore",suppressUpdates=TRUE)
```

SPADEVizR is available on [GitHub](https://github.com), at https://github.com/tchitchek-lab/SPADEVizR. Its installation can be done via the `devtools` package using the following commands:

```r
install.packages('devtools')
library("devtools")
install_github('tchitchek-lab/SPADEVizR')
```

Once installed, SPADEVizR can be loaded using the following command:

```r
library("SPADEVizR")
```



# <a name="loading_data"/> 3. Importing automatic gating results

## <a name="loading_SPADE_data"/> 3.1 Importing automatic gating results from SPADE

The `importSPADEResults()` function allows you to import cell clustering results generated by the SPADE algorythm. The function returns a `Results` object. Such import can be done using the following command:


```
## [START] - extracting SPADE results
## ImMemoryB-#00008_[MARKERSET10]_K070_P025
## FCS files import:
## 	archsin transform...
## 	compute quantiles...
## 	reading SPADE results...
## [END] - extracting SPADE results
```

```
## [START] - extracting SPADE results
## ImMemoryB-#00008_[MARKERSET10]_K070_P025
## FCS files import:
## 	archsin transform...
## 	compute quantiles...
## 	reading SPADE results...
## [END] - extracting SPADE results
```

*Optional:* The markers of SPADE results can be renamed using a dataframe (called dictionary). Such kind of dataframe must have in the first column the original marker names and have in the second column the new marker names. For instance, a dictionary can be loaded using the following command: 


```r
dict <- read.table("C:/Users/gg248485/Downloads/headerMAC.txt",sep="\t",header = TRUE)
head(dict, n = 4)
```

```
##         metal      marker
## 1        Time        time
## 2 Cell_length cell_length
## 3   (Rh103)Di empty-Rh103
## 4   (Xe131)Di empty-Xe131
```

Once a dictionary has been defined, SPADE cell clustering results can be loaded in the same way as above and by using the `dictionary` parameter. Specific markers can be excluded from the import procedure by providing their names to the `exclude.markers` parameter. the `quantile.heuristic` parameter (set by default to `TRUE`) can be used to approximate the computation of maker range quantiles. For instance, an import of a SPADE using a dictionary and by excluding the "empty-Rh103", "empty-Rh103", "empty-Rh103" and "empty-Rh103" markers can be done using the following command:

```r
results  <- importSPADEResults("C:/Users/gg248485/Desktop/SPADEVizR.data/HuMa-MAC",
							   dict               = dict,
							   quantile.heuristic = TRUE,
							   exclude.markers    = c("empty-Rh103"))
## [START] - extracting SPADE results
## HuMa-MAC
## FCS files import:
## 	archsin transform...
## 	compute quantiles...
## 	reading SPADE results...
## [END] - extracting SPADE results
```

## <a name="loading_other_data"/> 3.2 Importing automatic gating results from other algorithms
Moreover, SPADEVizR functionnalities can be used from results obtained from any cell clustering algorithms, using the `importX()` function. This function return a `Results` object. The function takes 2 dataframes in parameters: `cells.count` and `marker.expressions`.

 * `cells.count` is a dataframe containing the number of cells for each cluster of each sample. This dataframe must be formated as bellow: 

cluster |sample1 | sample2| ...
--------|--------|--------|----
c1      |749     | 5421   | ...
c2      |450     | 412    | ...
...     |...     | ...    | ...

 * `marker.expressions` is a dataframe containing the marker median expressions for each cluster of each sample. This dataframe must be formated as bellow:

sample   | cluster  | marker1  | marker2| ...
---------|----------|----------|--------|-----
sample1  |    c1    | 0.2      | 0.3    | ...
sample1  |    c2    | 0.1      | 0.3    | ...
sample2  |    c1    | 0.5      | 2.3    | ...
sample2  |    c2    | 1        | 1.3    | ...
...      |   ...    | ...      | ...    | ...


For instance, an import of cell clustering results obtained from as previously described can be done using the follwing command:

```r
#marker.expressions <- read.delim("PATH")
#cells.count        <- read.delim("PATH")
#results_other      <- importX(cells.count =  cells.count, marker.expressions = marker.expressions)
```

`Results` objects can be used by all functions excepting the `TreeViewer`, `PhenoViewer` (Heatmap viewer)and `classifyPhenoProfiles` which only accept a `SPADEResult` object.

# <a name="object_structures"/> 4. Object structures
The `print` and `show` functions are available for all objects of this package.

## <a name="object_structure_results"/> 4.1 Results object 
The `Results` object is a S4 object containing the count matrix and the cluster phenotypes. It is to note that `Results` is a super classe of the `SPADEResult` (defined in next subsection).

Different slots are available for a given `Results` object:

Slot               | Description
-------------------|----------------------------------------------------------------------------------------
cells.count        | a dataframe containing the number of cells for each cluster of each sample
marker.expressions | a numerical dataframe containing marker median expressions for each cluster of each sample
sample.names       | a character vector containing the sample names
marker.names       | a character vector containing the markers names
cluster.number     | a numeric specifying the number of cell clusters

## <a name="object_structure_SPADE_results"/> 4.2 SPADEResults object 
The `SPADEResults` object is a S4 object containing the clustering results from SPADE. It is to note that this object extend the `Results` object and contains additional slots related to SPADE data.

Different slots are available for a given `SPADEResults` object:

Slot          | Description            | Inherited
<<<<<<< HEAD
--------------|-------------------------------------------------------|--------
=======
--------------|------------------------|----------
>>>>>>> 7772a92c5202c701a0487ecf7d3d1ed73cc3e048
cells.count        | a dataframe containing the number of cells for each cluster of each sample                 | &#9745;
marker.expressions | a numerical dataframe containing marker median expressions for each cluster of each sample | &#9745;
sample.names       | a character vector containing the sample names                                             | &#9745;
marker.names       | a character vector containing the markers names                                            | &#9745;
cluster.number     | a numeric specyfing the number of clusters                                                 | &#9745;
use.raw.medians    | a logical specifying if the marker expressions correspond to the raw or transformed data   | &#9746;
marker.clustering  | a logical vector specifying markers that have been used during the clustering procedure    | &#9746;
fcs.files          | a character vector containing the absolute path of the original FCS files                  | &#9746;
quantiles          | a numeric matrix containing the quantiles for each marker of each cluster                  | &#9746;
graph              | an igraph object containing the SPADE tree                                                 | &#9746;
graph.layout       | a numeric matrix containing the layout of the SPADE tree                                   | &#9746;

## <a name="object_structure_AC"/> 4.3 Abundant Clusters (AC object)
The `AC` object is a S4 object containing the main information related to the abundant clusters, that is to say xxx, identify by the [`computeAC`](#stat_function_computeAC) function.  

Different slots are availables for a given `AC` object:

Slot       | Description
-----------|----------------------------------------------------------------------------------------
sample.names       | a character vector containing the samples used to compute the abundant clusters
cluster.size       | a numeric vector containing the number of cells ( -- sum of all samples -- ) for each cluster
use.percentages    | a logical specifying if computation was performed on percentage of cell abundance
method.name        | a character containing the name of the statistical test used to identify the abundant clusters
method.adjust      | a character containing the name of the multiple correction method used (if any)
th.mean            | a numeric value specifying the mean threshold
th.pvalue          | a numeric value specifying the p-value threshold
result             | a data.frame containing for each cluster (first column): the mean (second column) and the standard deviation (third column) of the biological condition, the associated p-value (fourth column) and a logical (fifth column) specifying if the cluster is significantly abundant.

## <a name="object_structure_DEC"/> 4.4 Differentially Enriched Clusters (DEC object)
The `DEC` object is a S4 object containing the main information related to the differentially enriched clusters, that is to say xxx, identify by the [`computeDEC`](#stat_function_computeDEC)  

Different slots are available for a given `DEC` object:

Slot       | Description
-----------|----------------------------------------------------------------------------------------
sample.cond1       | a character specifying the names of the samples of the first biological condition
sample.cond2       | a character specifying the names of the samples of the second biological condition
cluster.size       | a numeric vector containing number of cells ( -- sum of all samples -- ) for each cluster
use.percentages    | a logical specifying if computation was performed on percentage of cell abundance
method.name        | a character containing the name of the statistical test used to identify the DEC
method.adjust      | a character containing the name of the multiple correction method used (if any)
method.paired      | a logical indicating if the statistical test have been performed in a paired manner
th.fc              | a numeric value specifying the foldchange threshold
th.pvalue          | a numeric value specifying the p-value threshold
result             | a data.frame containing for each cluster (first column): the fold change (second column) and the standard deviation (third column) for the first biological condition, the fold change (fourth column) and the standard deviation (fifth column) for the second biological condition, the associated p-value (sixth column) and a logical (seventh column) specifying if the cluster is significantly differentially enriched.

## <a name="object_structure_CC"/> 4.4 Correlated Clusters (CC object)
The `CC` object is a S4 object containing object containing the main information related to the correlated clusters, that is to say xxx, identify by the [`computeCC`](#stat_function_computeCC)  

Different slots are available for a given `CC` object:

Slot       | Description
-----------|----------------------------------------------------------------------------------------
sample.names       | a character vector containing the samples used to compute correlated clusters
variables          | a numeric vector containing the expression values of the associated variable
cluster.size       | a numeric vector containing number of cells ( -- sum of all samples -- ) for each cluster
use.percentages    | a logical specifying if computation was performed on percentage of cell abundance
method.name        | a character containing the name of the statistical test used to identify the CC
method.adjust      | a character containing the name of the multiple correction method used (if any)
th.correlation     | a numeric value specifying the correlation threshold (R)
th.pvalue          | a numeric value specifying the p-value threshold
result             | a data.frame containing for each cluster (first column): the coefficiant of correlation R (second column) , the associated p-value (third column) and a logical (fourth column) specifying if the cluster is significantly correlated.

## <a name="object_structure_PhenoProfiles"/> 4.6 Classification of clusters based on their phenotype profiles (PhenoProfiles object)
The `PhenoProfiles` object is a S4 object containing the main information related of the cluster classification obtained by the [`classifyPhenoProfiles`](#stat_function_classifyPhenoProfiles)

Different slots are available for a given `PhenoProfiles` object:

Slot               | Description
-------------------|----------------------------------------------------------------------------------------
class.number       | a numeric value specifying the number of classes
method.name        | a character specifying the method used to classify cluster
method.parameter   | a named list of parameters used by the classification method
cluster.size       | a numeric vector containing the number of cells associated with each cluster (-- sum of all samples --)
( XXX) cluster.number| a numeric value specifying the number of clusters
classes            | a two column dataframe with the cluster in first colunm and corresponding classe in the second colunm

## <a name="object_structure_EnrichmentProfiles"/> 4.7 Classification of clusters based on their enrichment profiles (EnrichmentProfiles object)
The `EnrichmentProfiles` object is a S4 object containing the main information related of the cluster classification obtained by the [`classifyEnrichmentProfiles`](#stat_function_classifyEnrichmentProfiles)  

Different slots are available for a given `PhenoProfiles` object:

Slot               | Description
-------------------|----------------------------------------------------------------------------------------
class.number       | a numeric value specifying the number of classes
method.name        | a character specifying the method used to classify cluster
method.parameter   | a list of parameters used by the selected method
cluster.size       | a numeric vector with the number of cell in each cluster (-- sum of all samples --)
( XXX) cluster.number | a numeric value specifying the number of clusters
classes            | a two column dataframe with the cluster in first colunm and corresponding classe in the second colunm 


# <a name="stat_functions"/> 5. Statistical analysis

## <a name="stat_function_computeAC"/> 5.1 Computation of abundant clusters (computeAC function)
The `computeAC` function identify clusters with a number of associated cells statistically different from zero in a biological condition. 

The `computeAC` function takes as parameter a `SPADEResult` object and a named logical vector `condition` specifying the samples to use in the statistical computation. This named vector must provide the correspondence between samples (in names) and logical values (`TRUE` or `FALSE`). Significant abundant clusters are caracterized by two thresholds: the p-value threshold (`th.pvalue` parameter, set by default to `0.05`) and the mean abundance threshold (`th.mean` parameter, set by default to `0`).

The`computeAC` function returns a plotable `AC` object (see [abundantClustersViewer](#abundant_clusters_viewer_function)).

```r
condition <- c(TRUE,TRUE,FALSE,FALSE,TRUE)
names(condition) <- results@sample.names
resultsAC <- computeAC(results,condition = condition, th.pvalue = 0.01, th.mean = 2)
## [START] - computing ACs
## Sampled used :
## BB078 BPD19 concat_RightCells
## BB231 BPD19 concat_RightCells
## BD620 BPD19 concat_RightCells
## [END] - computing ACs
print(resultsAC)
## Object class: Abundant Clusters (AC)
## Samples : BB078 BPD19 concat_RightCells; BB231 BPD19 concat_RightCells; BD620 BPD19 concat_RightCells
## Use matrix of percent : TRUE
## Statiscal test used is : t.test
## Adjusted : none
## pvalue threshold :  0.01
## mean threshold :  2
## see @result slot to get further informations...
```

## <a name="stat_function_computeDEC"/> 5.2 Computation of differentially enriched clusters (computeDEC function)
The `computeDEC` function identify clusters with a number of associated cells statistically different between two biological conditions. 

The `computeDEC` function takes as parameter, a `SPADEResult` object and a named numeric vector `conditions` specifying the samples to consider in the two conditions. This named vector must provide the correspondence between samples (in names) and conditions (`1` to specify the first biological condition, `2` to indicate the second biological condition and `NA` otherwise). Differentially enriched clusters are caracterized by two thresholds: the p-value threshold (`th.pvalue` parameter, set by default to `0.05`) and the fold change threshold (`th.fc` parameter, set by default to `1`).

The `computeDEC` function returns a plotable `DEC` object (see [volcanoViewer](#volcano_viewer_function)).

```r
conditions <- c(1,1,2,2,1)
names(conditions) <- results@sample.names
resultsDEC <- computeDEC(results,conditions = conditions, th.pvalue = 0.01, th.fc = 3)
## [START] - computing DECs
## cond1:
## BB078 BPD19 concat_RightCells
## BB231 BPD19 concat_RightCells
## BD620 BPD19 concat_RightCells
## cond2:
## BC641 BPD19 concat_RightCells
## BD619 BPD19 concat_RightCells
## [END] - computing DECs
print(resultsDEC)
## Object class: Differentially Enriched Clusters (DEC)
## Sample of Condition 1 : BB078 BPD19 concat_RightCells; BB231 BPD19 concat_RightCells; BD620 BPD19 concat_RightCells
## Sample of Condition 2 : BC641 BPD19 concat_RightCells; BD619 BPD19 concat_RightCells
## Use matrix of percent : TRUE
## Statiscal test used is : t.test
## Adjusted : none
## Paired : FALSE
## pvalue threshold :  0.01
## foldchange threshold :  3
## see @result slot to get further informations...
```

## <a name="stat_function_computecc"/> 5.3 Computation of correlated clusters (computeCC function)
The `computeCC` function identify clusters correlated with an external biological variable. 

The `computeCC` function takes as parameter, a `SPADEResult` object and a named numeric vector `variable` specifying the expression values of the external biological variable. This named vector must provide the correspondence between samples (in names) and the expression values (`NA` to exclude this sample from analysis). Significant correlated clusters are caracterized by two thresholds: the p-value threshold (`th.pvalue` parameter, set by default to `0.05`) and the correlation (R) threshold (`th.correlation` parameter, set by default to `0.7`).

The `computeCC` function returns a plotable `CC` object (see [Correlated Clusters Viewer](#correlated_clusters_viewer_function))

```r
variable <- c(8,1.7,4,23,10)
names(variable) <- results@sample.names
resultsCC <- computeCC(results, variable = variable, th.pvalue = 0.01, th.correlation = 0.7)
## [START] - computing CCs
## [END] - computing CCs
print(resultsCC)
## Object class: Correlated Clusters (CC)
## Samples : BB078 BPD19 concat_RightCells; BB231 BPD19 concat_RightCells; BC641 BPD19 concat_RightCells; BD619 BPD19 concat_RightCells; BD620 BPD19 concat_RightCells
## Phenotypic variables : 8; 1.7; 4; 23; 10
## Use matrix of percent : TRUE
## Statiscal test used is : pearson
## Adjusted : none
## pvalue threshold : 
##  0.01
## correlation threshold : 
##  0.7
## see @result slot to get further informations...
```

## <a name="stat_function_classifyPhenoProfiles"/> 5.4 Classification of clusters based on their phenotype profiles (classifyPhenoProfiles function)
The `classifyPhenoProfiles` function classifies each cluster in different classes based on their marker expressions. Differents clustering methods are available among hierarchical clustering, kmeans, eigen vector decomposition, clique clustering. This function can only handle `SPADEResults` objects (but not a `Results` object). 

The `classifyPhenoProfiles` function takes a `SPADEResults` object and a `method` parameter among:

 * "hierarchical_k" (by default): performs a hierarchical clustering by specifying the desired number of classes. The number of desired classes needs to be specify using the `class.number` parameter.
 * "hierarchical_h" performs a hierarchical clustering by specifying the cut height. The height where the hierarchical tree should be cut needs to be specify using the `hierarchical.correlation.th` parameter.
 * "kmeans": performs a kmeans clustering. The number of desired classes needs to specify using the `class.number` parameter.
 * "eigencell": performs a eigen vector decomposition. The correlation coefficient cutoff needs to be specify using the `eigencell.correlation.th` parameter.
 * "clique": performs a clique clustering. The correlation coefficient cutoff needs to be specify using the `clique.correlation.th` parameter.

For instance, a classification of cell clusters based on their phenotypes can be performed using the following command: 

```r
#PhenoProfiles <- classifyPhenoProfiles(results, method.name = "kmeans", class.number = 10)
#print(PhenoProfiles)
```
The `classifyPhenoProfiles` function returns a `PhenoProfiles` object containing mainly the cluster associations.

## <a name="stat_function_classifyEnrichmentProfiles"/> 5.5 Classification of clusters based on enrichment profiles (classifyEnrichmentProfiles function)
The `classifyEnrichmentProfiles` function classifies each cluster in different classes based on their the number of cells associated with each sample (enrichment profiles). Differents clustering methods are available among hierarchical clustering, kmeans, eigen vector decomposition, clique clustering.

The `classifyEnrichmentProfiles` function takes a `Results` or `SPADEResult` object and a `method` parameter among:

 * "hierarchical_k" (by default): performs a hierarchical clustering by specifying the desired number of classes. The number of desired classes needs to be specify using the `class.number` parameter.
 * "hierarchical_h" performs a hierarchical clustering by specifying the cut height. The height where the hierarchical tree should be cut needs to be specify using the `hierarchical.correlation.th` parameter.
 * "kmeans": performs a kmeans clustering. The number of desired classes needs to be specify using the `class.number` parameter.
 * "eigencell": performs a eigen vector decomposition. The correlation coefficient cutoff needs to be specify using the `eigencell.correlation.th` parameter.
 * "clique": performs a clique clustering. The correlation coefficient cutoff needs to be specify using the `clique.correlation.th` parameter.

For instance, a classification of cell clusters based on their enrichment profiles can be performed using the following command: 

```r
#EnrichmentProfiles <- classifyEnrichmentProfiles(results, method.name = "kmeans", class.number = 10)
#print(EnrichmentProfiles)
```
The `classifyEnrichmentProfiles` function returns a `EnrichmentProfiles` object containing mainly the cluster associations.

# <a name="viewer_functions"/> 6. Vizualisation of statistical results

## <a name="abundant_clusters_viewer_function"/> 6.1 Abundant Clusters Viewer
The `abundantClustersViewer` function is used to visualize identified abundant clusters, that is to say results contained in an [AC objet](#object_structure_AC). This representation displays the p-value (shown as -log10(p-value) in X-axis) and the mean (Y-axis) of cells abundance in a two dimensional visualization. Each dot represents a cluster and abundant clusters are shown in red. The size of dots is proportional to the number of associated cells (--all samples are considered--). 

For instance, results contained in an `AC` objet can be shown using the following command:

```r
plot(resultsAC)
```

<img src="README.figures/AbundantClusters-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="volcano_viewer_function"/> 6.2 Volcano Viewer
The `volcanoViewer` function is used to visualize identifed differentially enriched clusters, that is to say results contained in a [DEC objet](#object_structure_DEC). The volcano plot [9] representation displays the p-value (shown as -log10(p-value) in Y-axis) and the foldchange (in X-axis) of abundance cells in a two dimensional visualization. By default, the foldchange is represented with a log2 transformation (which can be change using the `fc.log2` parameter). Each dot represents a cluster and differentially enriched clusters are shown in red. The size of dots is proportional to the number of associated cells (--all samples are considered--).

For instance, results contained in an `DEC` objet can be shown using the following command:

```r
plot(resultsDEC, fc.log2 = FALSE)
```

<img src="README.figures/VolcanoViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="correlated_clusters_viewer_function"/> 6.3 Correlated Clusters Viewer
The `correlatedClustersViewer` function is used to visualize identifed correlated clusters, that is to say results contained in a [CC objet](#object_structure_CC). This representation show the p-value (shown as -log10(p-value) in Y-axis) and the correlation coefficient (shown in X-axis) of correlated clusters in a two dimensional visualization. Each dot represents a cluster and correlated clusters are shown in red. The size of dots is proportional to the number of associated cells (--all samples are considered--).

For instance, results contained in an `CC` objet can be shown using the following command:

```r
plot(resultsCC)
```

<img src="README.figures/CorrelatedClusters-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="profile_viewer_function"/> 6.4 Profile Viewer
The `Profile Viewer` function is used to visualize classification results. This function can be used either on `EnrichmentProfiles` object or `PhenoProfiles` object and show classes of as a graph. In this graph, nodes represent clusters and edges connect clusters sharing the same class. The size of nodes is proportional to the number of associated cells.

-- how to display clique cluster of eigen cell clustering ? --

For instance, a `EnrichmentProfiles` object or `PhenoProfiles` object can be shown using the following command:

```r
#plot(EnrichmentProfiles)
```

# <a name="#viewer_functions"/> 7. Miscellaneous vizualisations

## <a name="cluster_viewer_function"/> 7.1 Cluster Viewer
The `clusterViewer` function returns a parallel coordinates respresentation for each cluster of each marker (clustering markers are shown in bold). This visualization allows to compare marker expressions between clusters and also with mean and quantiles expressions. The `show.mean` parameter specify the kind of information to display in the cluster viewer. The `both` value will show marker expressions for each samples together with mean expressions (black dashed line), `only` value will show only the mean expressions of samples and `none` value will mask mean expressions. In order to restrain vizualisation to specific markers and cluster, it is possible to provide a vector of marker names to `markers` parameter or in the same way a vector of clusters to `clusters` parameter.
<!-- 
By default all clusters will be computed, thus to avoid heavy processing, it is recommended to provide the vector of cluster to analyze in the "clusters" parameter.
We strongly advise to exclude "cell_length","FileNum","density","time" of this analysis because this data are pretty higher than other markers values.
-->  
For instance, the following command describe how to exlude "cell_length","FileNum","density" and "time" markers and then visualize cluster 1 and 8 :

```r
markers <- setdiff(sample(results@marker.names),c("cell_length","FileNum","density","time"))
gridExtra::grid.arrange(grobs = clusterViewer(results,clusters = c(1,8),show.mean = "both",markers = markers))
```

<img src="README.figures/ClusterViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="pheno_viewer_function"/> 7.2 Pheno Viewer
The `phenoViewer` function returns an heatmap of expression scores for each marker of each cluster.
<!-- Scores are cumputed for each makers between the lower bound and the upper bound of quantiles. -->
This function can only handle `SPADEResults` objects (but not a `Results` object). 


```r
phenoViewer(results)
## [START] - computing PhenoTable
## [END] - computing PhenoTable
```

<img src="README.figures/PhenoViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="tree_viewer_function"/> 7.3 Tree Viewer
The `treeViewer` function improve the SPADE tree showing the number of cells for each cluster (node of the tree) and if DEC, AC or CC object is provided (with parameter: `stat.object`), it highlights significant clusters. This function can only handle `SPADEResults` objects (but not a `Results` object). 


```r
treeViewer(results, stat.object = resultsDEC)
```

<img src="README.figures/TreeViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="kinetic_viewer_function"/> 7.4 Kinetic Viewer
The `kineticsViewer` function aims to represent the evolution of cells abundance for each cluster across the time for one or several individuals. It needs the `SPADEResult` object and to specify the `assignments` parameter. This parameter needs a dataframe with sample names in row names and 2 others column providing the timepoints and individuals associated with samples. 

<!-- By default all clusters will be computed, thus to avoid heavy processing, it is recommended to provide the vector of cluster ID to analyze in the "clusters" parameter. -->



```r
assignments <- data.frame(row.names = results@sample.names,
						  timepoints = c(1,1.5,2.1,2.4,3.2), 
						  individuals = c("HUMAN","MACAQUE","MACAQUE","HUMAN","HUMAN"))

gridExtra::grid.arrange(grobs = kineticsViewer(results, assignments = assignments, clusters=c(1,3)))
## These clusters will be cumpute :
## 1	3
```

<img src="README.figures/KineticViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="boxplot_viewer_function"/> 7.5 Boxplot Viewer
The `boxplotViewer` function aims to compare cell enrichment of a cluster across several biological conditions. Biological `conditions` are encoded in a named vector with samples in rownames providing correspondence with a biological condition (numeric or character are allowed). 


```r
conditions <- c(1,2,2,1,1)
names(conditions) <- results@sample.names
gridExtra::grid.arrange(grobs = boxplotViewer(results, label = TRUE, conditions = conditions, clusters=c(1,3)))
## These clusters will be cumpute :
## 1	3
```

<img src="README.figures/BoxplotViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="biplot_viewer_function"/> 7.6 Biplot Viewer
The `biplot` function can be use to visualize a two dimensions plot with an x-axis marker (`x.marker1`) and a y axis marker (`y.marker2`). It allows to filter cells point by clusters and samples with `clusters` and `samples` parameters. Moreover, cells can be shown for of all samples merged or with a facet for each sample.
When too cells dots are displayed, it can require heavy computation, to resolve this issue, you can resample your data to show less points with the `resample.ratio` parameter.

This function can only handle `SPADEResults` objects (but not a `Results` object). 

```r
#gridExtra::grid.arrange(biplot())
```

## <a name="distogram_viewer_function"/> 7.7 Distogram Viewer
The `distogramViewer` function displays correlations between all markers.


```r
distogramViewer(results)
```

<img src="README.figures/DistogramViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="streamgraph_viewer_function"/> 7.8 Streamgraph Viewer
The `streamgraphViewer` function aims to represent the evolution of cells abundance of all clusters selected using `clusters` parameter. To select and order sample, you should provide to `order` parameter, a named vector with samples in rownames providing correspondence with a position (integer) or NA to exclude a sample. 


```r
order <- c(1,3,NA,2,4)
names(order) <- results@sample.names
streamgraphViewer(results,order = order, clusters = c(1,2,3,14,5))
## These clusters will be cumpute :
## 1	2	3	14	5
```

<img src="README.figures/StreamgraphViewer-1.png" title="" alt="" style="display: block; margin: auto;" />

## <a name="MDS_viewer_function"/> 7.9 MDS Viewer (Multidimensional Scaling)
MDS methods aim to represent the similarities and differences among high dimensionality objects into a space of a lower dimensions, generally in two or three dimensions for visualization purposes [8]. In MDS representations, the Kruskal Stress (KS) indicates the amount of information lost during the dimensionality reduction process.

The `MDSViewer` function generates a MDS (Multidimensional Scaling) representation either of clusters or samples based on `space` parameter. If `space = "sample"`, `assignments` is required (as for [Kinetic Viewer](kinetic_viewer_function)).


```r
MDSViewer(results, space = "clusters", assignments = assignments,clusters = c(2,3,4,5,6,7,9,20,41,52,44,74,21,49,22))
## These clusters will be cumpute :
## 2	3	4	5	6	7	9	20	41	52	44	74	21	49	22
## MDS computation
## done
```

<img src="README.figures/MDSViewer_clusters-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
MDSViewer(results, space = "samples", assignments = assignments, clusters = c(2,3,4,5,6,7,9,20,41,52,44,74,21,49,22))
## These clusters will be cumpute :
## 2	3	4	5	6	7	9	20	41	52	44	74	21	49	22
## MDS computation
## done
```

<img src="README.figures/MDSViewer_samples-1.png" title="" alt="" style="display: block; margin: auto;" />

# <a name="license"/> 9. Export 
The `export()` function is available for all objects of this package. Use this function to export an object to a tab separated file able to be open with Microsoft Excel&copy; or with Libre Office Calc.


```r
export(AC,filename = "DESIRED_PATH.txt")
```

# <a name="license"/> 10. Generate report 
The `generateReport()` function allows to easily generate a PDF file containing all desired plots. In order to select the plots to include in the report, you could provided a character vector containing the names of desired plots among: 

 * "pheno" to include Pheno Viewer
 * "kinetic" to include Kinetic Viewer (need to provide the `assignement` parameter in the same way as [Kinetic Viewer](kinetic_viewer_function) function)
 * "cluster" to include Cluster Viewer
 * "kinetic_cluster" to include Kinetic Viewer juxtaposed with Cluster Viewer
 * "tree" to include Tree Viewer
 * "disto" to include Distogram Viewer
 * "stream" to include Streamgraph Viewer
 * "MDS" to include MDS Viewer (if `assignement` parameter is provide, it will plot the 2 spaces: "clusters" and "samples" of MDS represnetation, else only "clusters" space will be plotted)

The report will follows the order of plots names in the vector.

You can also provided a vector of stat objects (AC, DEC, CC) with the `stat.objects` parameter and a vector of profile.objects (`PhenoProfile` and `EnrichmentProfile`) with the `profile.objects` parameter


```r
generateReport(results,PDFfile = "DESIRED_PATH.pdf", assignments = assignments, report = c("pheno", "kinetic_cluster","tree", "disto","stream","MDS"))
```

*Generating a big report can take a minute or more.*

# <a name="license"/> 11. License
SPADEVizR is freely distributed under the GLP-3 license.

# <a name="references"/> 12. References 
[1] - Bendall, S. C., Simonds, E. F., Qiu, P., Amir, E. D., Krutzik, P. O., Finck, R., . Nolan, G. P. (2011). Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. Science (New York, N.Y.), 332(6030), 687-96.

[2] - Gregori, G., Rajwa, B., Patsekin, V., Jones, J., Furuki, M., Yamamoto, M., & Paul Robinson, J. (2014). Hyperspectral cytometry. Current Topics in Microbiology and Immunology, 377, 191-210.

[3] - http://www.flowjo.com/

[4] - https://www.cytobank.org

[5] - Qiu, P., Simonds, E. F., Bendall, S. C., Gibbs, K. D., Bruggner, R. V, Linderman, M. D., . Plevritis, S. K. (2011). Extracting a cellular hierarchy from high-dimensional cytometry data with SPADE. Nature Biotechnology, 29(10), 886-91.

[6] - Grammar of Graphics library http://ggplot2.org/

[7] - Ellis B, Haaland P, Hahne F, Meur NL, Gopalakrishnan N, Spidlen J and Jiang M. flowCore: flowCore: Basic structures for flow cytometry data. R package version 1.34.7.

[8] - Kruskal, J. B., and Wish, M. (1978), Multidimensional Scaling, Sage University Paper series on Quantitative Application in the Social Sciences, 07-011. Beverly Hills and London: Sage Publications.

[9] - Cui X1, Churchill GA. (2003). Statistical tests for differential expression in cDNA microarray experiments. Genome Biol.
