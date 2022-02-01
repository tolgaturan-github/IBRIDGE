generate_QC<-function(sce1){
	library(scater)
	sce1<-addPerFeatureQC(addPerCellQC(sce1))
	print(paste0("ColSums - Outlier cell count  lower than 5 NMADs is " , sum(isOutlier(colData(sce1)$sum, nmads = 5, type = "lower",log=TRUE)))) 
	print(paste0("ColSums - Outlier cell count  higher than 5 NMADs is " , sum(isOutlier(colData(sce1)$sum, nmads = 5, type = "higher",log=TRUE))))
	print(paste0("NonZero - Outlier cell count  lower than 5 NMADs is " , sum(isOutlier(colData(sce1)$detected, nmads = 5, type = "lower",log=TRUE))))
	print(paste0("Nonzero - Outlier cell count  higher than 5 NMADs is " , sum(isOutlier(colData(sce1)$detected, nmads = 5, type = "higher",log=TRUE))))

	rowData(sce1) <- cbind(rowData(sce1), norm_mean = calculateAverage(sce1), nexprs = nexprs(sce1, byrow = TRUE))
	sce1 <- sce1[rowData(sce1)$nexprs > 0, ]
	sce1}
