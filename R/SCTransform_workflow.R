SCTransform_workflow<-function(dir1){
	sce1<-generate_sce(dir1)
	sce1<-generate_QC(sce1)
	seu1<-SCTransform_normalization(sce1)	
 	seu1}

	
