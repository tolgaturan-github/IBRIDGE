generate_sce<-function(dir1, meta_path=FALSE){
	library(SingleCellExperiment)
	library(GEOquery)
	if(paste0(dir1, "/barcodes.tsv.gz")==list.files(paste0(dir1), full.names=TRUE)[1]){gunzip(paste0(dir1, "/barcodes.tsv.gz"))}
	if(paste0(dir1, "/features.tsv.gz")==list.files(paste0(dir1), full.names=TRUE)[2]){gunzip(paste0(dir1, "/features.tsv.gz"))}
	if(paste0(dir1, "/matrix.mtx.gz")==list.files(paste0(dir1), full.names=TRUE)[3]){gunzip(paste0(dir1, "/matrix.mtx.gz"))}
	bar1<-read.table(paste0(dir1, "/barcodes.tsv"), sep="\t")
	feat1<-read.table(paste0(dir1, "/features.tsv"), sep="\t")
	matrix1<-read.table(paste0(dir1, "/matrix.mtx"), sep=" ", skip=2)
	mat_list1<-split(matrix1, as.factor(as.character(matrix1[,2])))
	names(mat_list1)<-as.character(bar1$V1)[as.integer(names(mat_list1))]
	matrix2<-sapply(mat_list1, function(x) {x[match(rownames(feat1), as.character(x[,1])),3]})
	rownames(matrix2)<-feat1[,1]
	matrix3<-matrix2[,as.character(bar1[bar1[,1] %in% colnames(matrix2),1])]
	matrix4<-apply(matrix3, 2, function(x) ifelse(is.na(x), 0, x))
	if (is.character(meta_path)){
		meta1<-read.table(meta_path, sep="\t", header=TRUE)
		rownames(meta1)<-as.character(bar1[,1])
		sce1 <- SingleCellExperiment(assays = list(counts = matrix4), colData = meta1)}
	else {
		sce1 <- SingleCellExperiment(assays = list(counts = matrix4))}
	sce1
	}
	

