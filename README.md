IBRIDGE integrates bulk and single-cell datasets to profile granular heterogeneity of TIME
================
Tolga Turan
February 1, 2022

## INTRODUCTION:

Quantification of T cell infiltration levels in a tumor tissue, is an important component of Immuno-Oncology research. In bulk gene expression datasets, this can be approximated by using "Gene Set Enrichment" or "Cell Type Deconvolution" methods (with gene signatures such as ICR, TIS or IFNG, etc. ). However, reaching the same goal is not straigthforward in single-cell RNAseq datasets. Here, we report a novel method "IBRIDGE" that can integrate and *bridge* between bulk and single-cell gene expression data types to quantify T cell infiltration levels and classify the patients as "inflamed", "cold" or "unassigned/intermediate".

## INSTALLATION:

Install IBRIDGE and the packges it depends as follows:

``` {.r}
library(devtools)
install.packages("Seurat")
BiocManager::install("AUCell")
install_github("tolgaturan-github/IBRIDGE")
```

## USAGE:

As explained in the associated publication, malignant cancer cells are the only subpopulation that clusters by patient group. Therefore, we hypothesize that 'malignant' population can be used to assess T-cell infiltration levels from scRNASeq data by integrating the correlates of T-cell infiltration from the corresponding TCGA cohort. In this vignette, we utilized a Breast Cancer scRNAseq dataset (GSE176078) reported by Wu et al., 2021, to describe the IBRIDGE workflow.

Load necessary packages, the geneXcell count matrix and the patient metadata for the malignant subset of this dataset

``` {.r}
#Below test data are used to replicate this vignette but are excluded from the package due to space restrictions.
data("GSE176078_malignant_UMI_countmatrix", "GSE176078_malignant_metadata", package="IBRIDGE")
```

``` {.r}
library(IBRIDGE)
library(Seurat)
dim(GSE176078_malignant_UMI_countmatrix)
```

    ## [1] 29733 27915

``` {.r}
GSE176078_malignant_UMI_countmatrix[1:6, 1:6]
```

    ##        CID3586_TTCTCCTTCTGAAAGA CID3586_AACCATGGTTCAGGCC
    ## MALAT1                        0                       66
    ## MT-CO2                        0                       28
    ## MT-CO1                        0                       40
    ## MT-CO3                        0                       16
    ## RPLP1                        57                       86
    ## FTH1                         11                       42
    ##        CID3586_AACCGCGAGCGATGAC CID3586_AAGGAGCCATCTGGTA
    ## MALAT1                       88                       35
    ## MT-CO2                        0                        5
    ## MT-CO1                       10                        2
    ## MT-CO3                        3                        1
    ## RPLP1                        10                        7
    ## FTH1                          2                       10
    ##        CID3586_AAGGCAGAGATCCGAG CID3586_AAGTCTGAGCCAACAG
    ## MALAT1                       68                       89
    ## MT-CO2                       13                       27
    ## MT-CO1                       14                       22
    ## MT-CO3                        8                       10
    ## RPLP1                        53                       51
    ## FTH1                         26                       73

``` {.r}
table(GSE176078_malignant_metadata$Patient)
```

    ## 
    ##  CID3586  CID3921  CID3941  CID3948  CID3963  CID4066  CID4067 CID4290A 
    ##      698      442      197      261      221      786     2355     4082 
    ## CID44041  CID4461  CID4463  CID4465  CID4471  CID4495 CID44971 CID44991 
    ##      151      207      685      137     2178     1150     1622     4001 
    ##  CID4513  CID4515 CID45171  CID4523 CID4530N  CID4535 
    ##        8     2205      815     1170     2134     2410

Normalize count matrix using SCTransform based on Seurat workflow.

``` {.r}
seu1<-SCTransform_normalization(GSE176078_malignant_UMI_countmatrix, GSE176078_malignant_metadata$Patient)
```

To show that malignant population of single-cells clusters by patient, we can visualize patients vs. UMAP dimensions:

``` {.r}
library(ggplot2)
library(RColorBrewer)
DimPlot(seu1, group.by="Patient", size=0.5)+scale_colour_manual(values=c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(5, "Set2")))
```

![](README_files/figure-markdown_github/Visualize_clusters-1.png)

Now, let's identify IBRIDGE overlaps between single-cell and bulk expression data. Since this is a BRCA dataset, we pull TCGA features that correlate with inflamed and cold TCGA-BRCA samples. Then we can identify top highly variable genes form the single-cell data:

``` {.r}
IBRIDGE_features<-IBRIDGE_overlaps(seu1, "gene", "BRCA", "SCTransform")
lapply(IBRIDGE_features, head)
```

    ## $inflamed
    ## [1] "B2M"    "HLA-C"  "HLA-A"  "HLA-B"  "MUCL1"  "CALML5"
    ## 
    ## $cold
    ## [1] "NEAT1"   "COX6C"   "AZGP1"   "AGR2"    "PIP"     "SLC39A6"

After identifying the features we can apply geneset enrichment to see which of the cancer cells express inflamed/cold features. Then we can assign an inflamed/cold classification to each of the cells:

``` {.r}
seu1<-Classify_cells(seu1, IBRIDGE_features, "SCTransform", 100)
```

``` {.r}
table(seu1@meta.data$IBRIDGE_Class)
```

    ## 
    ##       Cold   Inflamed Unassigned 
    ##       8539       8540      10836

We can then visualize these classes on the UMAP dimensions:

``` {.r}
library(ggplot2)
DimPlot(seu1, group.by="IBRIDGE_Class", size=0.5)+scale_colour_manual(values=c("blue", "red", "gray"))
```

![](README_files/figure-markdown_github/Visualize_classes-1.png)

References:

Wu SZ, Al-Eryani G, et al. **A single-cell and spatially resolved atlas of human breast cancers.** *Nat Genet. 2021 Sep;53(9):1334-1347.*
