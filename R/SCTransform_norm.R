#SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures() functions from the Seurat package#


SCTransform_normalization<-function(sce1){
        library(Seurat)
       	seu1<-CreateSeuratObject(counts =as.matrix(counts(sce1))) 
        seu1@meta.data<-cbind(seu1@meta.data, data.frame(colData(sce1)))
	seu1<-SCTransform(object = seu1, verbose = FALSE)
        seu1<-RunPCA(seu1)
        seu1<-RunUMAP(seu1, dims=1:20)
        seu1<-FindNeighbors(seu1, dims=1:20)
        seu1<-FindClusters(seu1)
        seu1}

