#' Classify single-cells or in vitro cell lines
#' Classifies the input single-cell or in vitro cell line data into "inflamed", "cold" or "unassigned" bins
#' @param seurat_object Output of UMI count matrix normalized by Seurat workflows
#' @param IBRIDGE_features Output of IBRIDGE_overlaps function.
#' @param norm_method "SCTransform" or "NormalizeData". Defaults to "SCTransform"
#' @param n_cores Integer of length 1. Specifies the number of parallel cores running the function.
#'      Defaults to 1
#' @return A seurat_object with metadata updated with IBRIDGE classification 
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/tolgaturan-github/IBRIDGE}
#' @examples
#' seu1<-Classify_cells(seu1, IBRIDGE_feautures, "SCTransform", 100)
#' @export



Classify_cells<-function(seu1, IBRIDGE_features,norm_method="SCTransform", n_cores=1){

	nes1<-sc_geneset_enrich(seu1, IBRIDGE_features, norm_method, n_cores)
	nes_df<-data.frame(nes1[,1:2], Inflamed_tertiles=assign_tertiles(nes1[,1], "Inflamed"), Cold_tertiles=assign_tertiles(nes1[,2], "Cold"))
	nes_df$Class<-ifelse(nes_df$Inflamed_tertiles=="Inflamed_high"&nes_df$Cold_tertiles!="Cold_high", "Inflamed", ifelse(nes_df$Cold_tertiles=="Cold_high"&nes_df$Inflamed_tertiles!="Inflamed_high", "Cold", "Unassigned"))
	seu1@meta.data$IBRIDGE_Class<-nes_df$Class
	return(seu1)
}

