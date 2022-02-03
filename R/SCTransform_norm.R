#' SCTransform based scRNAseq Workflow 
#' This function accepts UMI count matrix and patient metadata as arguments 
#' and replaces NormalizeData(), ScaleData(), and FindVariableFeatures() functions from the Seurat package
#' The function also performs UMAP dimensio reduction and clustering
#' @param countmatrix Numeric matrix of UMI counts genes as rows and cell_ids as columns from malignant subset of the scRNAseq data
#'      This matrix is usually a subset of a scRNAseq dataset profiling the whole TIME.
#' @param metadata Character vector (1-dimension) of length matching the cell count for the input countmatrix
#'      Preferably all or at least 1 feature in each geneset must be present in the rownames of the expression matrix.
#' @return A seurat object with SCTransform normalized counts and PErson residuals as well as patient metadata and dimension reduction slots
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/tolgaturan-github/IBRIDGE}
#' @examples
#' seurat_object1<-SCTransform_normalization(countmatrix1, patient_metadata1)
#' @import Seurat
#' @export
#'


SCTransform_normalization<-function(countmatrix, metadata){
        library(Seurat)
       	seu1<-CreateSeuratObject(counts =countmatrix) 
        seu1@meta.data<-cbind(seu1@meta.data, Patient=metadata)
	seu1<-SCTransform(object = seu1, verbose = FALSE)
        seu1<-RunPCA(seu1)
        seu1<-RunUMAP(seu1, dims=1:20)
        seu1<-FindNeighbors(seu1, dims=1:20)
        seu1<-FindClusters(seu1)
        seu1}

