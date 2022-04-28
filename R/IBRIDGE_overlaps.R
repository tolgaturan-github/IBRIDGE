#' IBRIDGE_Overlaps
#' Integrates bulk with single-cell (or cell line) gene expression data and identifies features that can classfy cells into "inflamed" and "cold" bins.
#' @param seurat_object Output of UMI count matrix normalized by Seurat workflows.
#' @param type type of gene identifier the dataset has. "ens" for ensembl IDs and "gene" for HUGO gene exymbols. Defaults to "ens"
#' @param cohort one of 28 solid tumor cohorts from TCGA database. The cancer cohort from which the single-cell gene expression dataset is generated.
#' @param norm_method "SCTransform" or "NormalizeData". Defaults to "SCTransform"
#' @return list of 2 elements. One for "inflamed" and one for "cold" features.
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/tolgaturan-github/IBRIDGE}
#' @examples
#' IBRIDGE_features<-IBRIDGE_overlaps(seu1, "LUAD", "SCTransform")
#' @import AUCell
#' @import Matrix
#' @export
#'

IBRIDGE_overlaps<-function(seu1, type="ens", cohort=c("ACC", "BLCA", "CESC", "CHOL", "COAD", 
                                                               "ESCA", "GBM",  "HNSC", "KICH", "KIRC", "KIRP", 
                                                               "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", 
                                                               "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC", 
                                                               "UCS"), norm_method="SCTransform"){ 
    if (type=="ens"){
        bulk_features<-ens_list[[cohort]]
    }else if (type=="gene"){
        bulk_features<-genesymb_list[[cohort]]
    }
    lapply(bulk_features, head)
    variable1 <- VariableFeatures(seu1)
    bulk_features[[1]]<-intersect(bulk_features[[1]], variable1)
    bulk_features[[2]]<-intersect(bulk_features[[2]], variable1)
    
    data <- NULL
    if (norm_method == "SCTransform") {
        data <- seu1@assays$SCT@data
    } else if (norm_method == "NormalizeData") {
        data <- seu1@assays$RNA@data
    } else {
        stop("Normalization method is invalid! Please use 'SCTransform' or 'NormalizeData'") 
    }
    if (length(bulk_features[[1]]) == 0) {
	warn("No overlapped inflamed gene signature found")
 	inflamed_top100 <- c()
    } else {
    	inflamed.mean <- Matrix::rowMeans(data[bulk_features[[1]],])
    	inflamed_top100 <- head(bulk_features[[1]][order(inflamed.mean, decreasing = T)], n=100)
    }

    if (length(bulk_features[[2]]) == 0) {
	warn("No overlapped cold gene signature found")
	cold_top100 <- c()
    } else {
    	cold.mean <- Matrix::rowMeans(data[bulk_features[[2]],])
    	cold_top100 <- head(bulk_features[[2]][order(cold.mean, decreasing = T)][1:100], n=100)
    }
    IBRIDGE_features<-list(inflamed=inflamed_top100, cold=cold_top100)
    return(IBRIDGE_features)
}
