#' Geneset Enrichment (GSE) on scRNAseq Data
#' Calculates single cell GSE scores using "AUCell" method
#' @param seu1 Output of UMI count matrix normalized by Seurat workflows
#' @param geneset_list a List of gene/Transcript features.
#'      Preferably all or at least 1 feature in each geneset must be present in the rownames of the expression matrix.
#' @param norm_method "SCTransform" or "NormalizeData". Defaults to "SCTransform"
#' @param n_cores Integer of length 1. Specifies the number of parallel cores running the function.
#'      Defaults to 1
#' @return data.frame holding geneset enrichment values for each geneset (columns) and for each cell (rows).
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/tolgaturan-github/IBRIDGE}
#' @examples
#' enrichment_scores1<-sc_geneset_enrich(seu1, geneset_list1, "SCTransform", 100)
#' @import AUCell
#' @export
#'


sc_geneset_enrich<-function(seu1, geneset_list, norm_method="SCTransform",n_cores=1){
	library(AUCell)
	if (norm_method=="SCTransform"){
  		cells_rankings <- AUCell_buildRankings(seu1$SCT@data, BPPARAM=BiocParallel::MulticoreParam(n_cores), plotStats=FALSE, splitByBlocks=TRUE)}
	else if(norm_method=="NormalizeData"){
		cells_rankings <- AUCell_buildRankings(seu1$RNA@data, BPPARAM=BiocParallel::MulticoreParam(n_cores), plotStats=FALSE, splitByBlocks=TRUE)} 
else {
                stop("Invalid norm_method, please use 'SCTransform' or 'NormalizeData")
        }


	cells_AUC <- AUCell_calcAUC(geneset_list, cells_rankings, aucMaxRank=ceiling(0.10 * nrow(cells_rankings)))
	auc1 <- data.frame(t(cells_AUC@assays@data@listData$AUC))
	auc1}
