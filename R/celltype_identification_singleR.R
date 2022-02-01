celltype_id<-function(seu1, database=c("blueprint", "primarycellatlas","mousernaseq", "immgen"), type=c("ens", "gene"), cluster_only=TRUE){
	library(SingleCellExperiment)
	library(SingleR)
	library(celldex)
	library(scater)
	sce1<-logNormCounts(as.SingleCellExperiment(seu1))
	if("blueprint" %in% database){
		if (cluster_only){
                pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/hg38/hg38_ENS_to_Symb.txt", sep=",", header=TRUE)
		if (type=="ens"){
                        rownames(sce1)<-as.character(pl1$Gene.name)[match(rownames(sce1), pl1$Gene.stable.ID)]
			type<-"gene"
			}
                if (type=="gene"){}

		data1<-BlueprintEncodeData()
		main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
	 	fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
		seu1@meta.data$SingleR_bp_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
		seu1@meta.data$SingleR_bp_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]}
		else{
                pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/hg38/hg38_ENS_to_Symb.txt", sep=",", header=TRUE)
                if (type=="ens"){
                        rownames(sce1)<-as.character(pl1$Gene.name)[match(rownames(sce1), pl1$Gene.stable.ID)]
                        type<-"gene"
                        }
                if (type=="gene"){}

                data1<-BlueprintEncodeData()
                main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
                fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
                seu1@meta.data$SingleR_bp_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
                seu1@meta.data$SingleR_bp_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]
		main2<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main)
		fine2<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine)
		seu1@meta.data$SingleR_bp_main_cell<-data.frame(main2)$labels
		seu1@meta.data$SingleR_bp_fine_cell<-data.frame(fine2)$labels
		}
		}
	if("primarycellatlas" %in% database){
                if (cluster_only){
		pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/hg38/hg38_ENS_to_Symb.txt", sep=",", header=TRUE)
                if (type=="ens"){
                        rownames(sce1)<-as.character(pl1$Gene.name)[match(rownames(sce1), pl1$Gene.stable.ID)]
			type<-"gene"
		}
                if (type=="gene"){}
		data1<-HumanPrimaryCellAtlasData()
		main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
                fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
                seu1@meta.data$SingleR_prim_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
                seu1@meta.data$SingleR_prim_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]}
		else{
                pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/hg38/hg38_ENS_to_Symb.txt", sep=",", header=TRUE)
                if (type=="ens"){
                        rownames(sce1)<-as.character(pl1$Gene.name)[match(rownames(sce1), pl1$Gene.stable.ID)]
                        type<-"gene"
                        }
                if (type=="gene"){}

                data1<-HumanPrimaryCellAtlasData()
                main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
                fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
                seu1@meta.data$SingleR_prim_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
                seu1@meta.data$SingleR_prim_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]
                main2<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main)
                fine2<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine)
                seu1@meta.data$SingleR_prim_main_cell<-data.frame(main2)$labels
                seu1@meta.data$SingleR_prim_fine_cell<-data.frame(fine2)$labels
                }
                }

	if("mousernaseq" %in% database){
        	pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/mm11/Primary/GenCode/transcriptIDs_to_genes.txt", sep="\t", header=TRUE)
		if (type=="ens"){
                	rownames(sce1)<-as.character(pl1$Gene)[match(rownames(sce1), pl1$TranscriptId)]
			type<-"gene"
		}
        	if (type=="gene"){}
		data1<-MouseRNAseqData()
                main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
                fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
                seu1@meta.data$SingleR_mousernaseq_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
                seu1@meta.data$SingleR_mousernaseq_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]
		if(!cluster_only){
			seu1@meta.data$SingleR_mousernaseq_main_cell<-data.frame(SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main))$labels
			seu1@meta.data$SingleR_mousernaseq_fine_cell<-data.frame(SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine))$labels
			}}
	if("immgen" %in% database){
		pl1<-read.table("/bioinformatics/Users/turantx1/references/genomes/mm11/Primary/GenCode/transcriptIDs_to_genes.txt", sep="\t", header=TRUE)
		if (type=="ens"){
                        rownames(sce1)<-as.character(pl1$Gene)[match(rownames(sce1), pl1$TranscriptId)]
			type<-"gene"
		}
                if (type=="gene"){}
                data1<-ImmGenData()
                main1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main, clusters=seu1@meta.data$seurat_clusters)
                fine1<-SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine, clusters=seu1@meta.data$seurat_clusters)
                seu1@meta.data$SingleR_immgen_main<-data.frame(main1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(main1))), "labels"]
                seu1@meta.data$SingleR_immgen_fine<-data.frame(fine1)[match(seu1@meta.data$seurat_clusters, rownames(data.frame(fine1))), "labels"]
		if(!cluster_only){
                        seu1@meta.data$SingleR_immgen_main_cell<-data.frame(SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.main))$labels
                        seu1@meta.data$SingleR_immgen_fine_cell<-data.frame(SingleR(test = assays(sce1)$logcounts, ref = data1, labels = data1$label.fine))$labels
                        }}
	seu1}

