# Codes for single cell RNA sequencing
<h3>Install Packages</h3>

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("glmGamPoi")
BiocManager::install("DESeq2")
BiocManager::install("MAST")

install.packages(c("dplyr","Seurat","patchwork","data.table","glmGamPoi","sctransform"))
```

<h3>Call Library</h3>

```r
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(sctransform)
library(reticulate)
library(ggplot2)
library(data.table)
library(DESeq2)
library(MAST)
```
<h3>Call Database and check info</h3>

```r
E_all_combined <- readRDS("/media/data/naz/scmRNA_data/E_all_combined/E_all_combined_for_clustering.rds")
head(colnames(E_all_combined ))
table(E_all_combined $treatment)
```

<h3>Make cluster</h3>

```r
object <- FindNeighbors(object = object, dims = 1:40)
object <-FindClusters(object = object,graph.name = "RNA_snn")
object <- RunUMAP(object = object,graph = "RNA_snn")
DimPlot(object = object,reduction = "umap", group.by = "orig.ident")

E_all_combined_2 <-FindClusters(object = E_all_combined, graph.name = "RNA_snn", resolution = 0.2)
cluster <- DimPlot(object=E_all_combined_2, reduction="umap", label=TRUE, repel = T, pt.size = 0.5)
cluster

```
<h3>Modify cluster resulution on the basis of need</h3>

```r
E_11_13_combined_2 <-FindClusters(object = E_11_13_combined_1, graph.name = "RNA_snn", resolution = 0.18)
cluster <- DimPlot(object=E_11_13_combined_2, reduction="umap", label=TRUE)
cluster
```
<h3>Split / group cluster on the basis of treatment</h3>

```   r
cluster_overlapped_treatment <- DimPlot(E13_1, label = F, repel = T, pt.size = 0.5, group.by = "treatment", cols = c("#E04b41", "#41b0e0")) + ggtitle("Unsupervised clustering")
cluster_overlapped_treatment

cluster_by_treatment <- DimPlot(E_all_combined_2, reduction="umap", label=F, repel = T, pt.size = 0.1, split.by = "treatment", cols = c("#E04b41", "#41b0e0")) + ggtitle("Clustering by treatment")
cluster_by_treatment
```

<h3>Change the name of cluster created by seurat</h3>

``` r
Idents(E_11_13_combined_2) <- "RNA_snn_res.0.18"
E_11_13_combined_3 <- RenameIdents(E_11_13_combined_2, `0` = "CD14 Mono", `1` = "CD4 Naive T",
`2` = "CD4 Memory T", `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated",
`8` = "DC", `9` = "B Activated",`10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")
DimPlot(E_11_13_combined_3, label = TRUE)
``` 
<h3>Make dot plot on the basis of highly expressed gene from each cluster</h3>

```r
Idents(E_11_13_combined_3) <- factor(Idents(E_11_13_combined_3), levels = c("Mono/Mk Doublets",
    "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
    "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("Col1a1", "Col2a1", "Sox11", "Pax1", "Tbx3", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
#DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
DotPlot(E_11_13_combined_3, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8,
col.min=-2, col.max=2, scale.min = 0, scale.max = 100) +
    RotatedAxis()
```
<h3>Quick Check genes expression at each cluster, Cluster 4 is the myogenic</h3>

```
Idents(E_all_combined_1) <-"RNA_snn_res.0.15"
levels(E_all_combined_1)
p2 <- VlnPlot(E_all_combined_1, c("Scx","Msx","Runx2","Msx2","Dlx5","Zfp423","Pdgfra","Prrx1","Col1a1"), assay = "RNA", pt.size= 0.01, split.by = "treatment" , cols = c("#E04b41", "#41b0e0"))

#ggsave(plot = p2, 'E13_mRNA_CT_HFD_PC70_res0.2_vplot_fibro_adipo.tiff', dpi=300, width=8, height=6, path = "cluster_umap", limitsize = TRUE)

p2
```

<h3>Find conserved gene from targeted cluster</h3>

```
cluster5_conserved_markers <- FindConservedMarkers(E_all_combined_2,
                              ident.1 = 5,
                     	      grouping.var = "orig.ident",
                              only.pos = TRUE,
		              logfc.threshold = 0.25)
```





<h3>Using featurePLot to visulized the gene expression among all clusters</h3>

```r
p4<-FeaturePlot(E13_1, features = c("Zfp423","Pparg"), split.by = "treatment", min.cutoff = 0)

p3<-FeaturePlot(E13_1, features = c("Pdgfra","Tcf7l2"),split.by = "treatment")


ggsave(plot = p3, "E13_mRNA_CT_HFD_Norm_Scale_PC100_res0.2_cluster_split_featureplot_fibrogenic.tiff", dpi=600, width=10, height=10, path = "cluster_umap", limitsize = TRUE)

ggsave(plot = p4, "E13_mRNA_CT_HFD_Norm_Scale_PC100_res0.2_cluster_split_featureplot_adipogenic.tiff", dpi=600, width=10, height=10, path = "cluster_umap",

```

<h3>Subset Cluster</h3>

```r
Idents(E13_1) <-"RNA_snn_res.0.2"
adipogenic_cluster<-subset(x=E13_1, idents=c("8","13"))
saveRDS(adipogenic_cluster, file = "fibro_adipogenic_cluster/adipogenic_cluster.rds")
```


<h1>Pseudotiming</h1>
<h3>Call Packages</h3>

```r
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
set.seed(1234)
```

<h3>Make Pseudotiming</h3>

```r
myogenic_cluster<- readRDS("/media/data/Naz/data_analysis/E13_RNA_data/RDSdata/subcluster/fibrogenic_cluster.rds")
E11.5_myogenic_Diet <- DietSeurat(myofibro_cluster, counts=TRUE, data=TRUE, graph="umap", dimreducs = "umap")


E11.5_myogenic_Diet[["UMAP"]] <- E11.5_myogenic_Diet[["umap"]]  ## If THE UMAP from seurat needs to be kept, it only accepts "UMAP' or "PCA" and others. As we want to use the WNN umap, we need to name it as "UMAP" to be processed in Monocle3.



E11.5_myogenic_Diet.cds <- as.cell_data_set(E11.5_myogenic_Diet)
E11.5_myogenic_Diet.cds <- cluster_cells(cd=E11.5_myogenic_Diet.cds, reduction_method = "UMAP")
E11.5_myogenic_Diet.cds <- learn_graph(E11.5_myogenic_Diet.cds,use_partition = TRUE)

E11.5_myogenic_Diet.cds <- order_cells(E11.5_myogenic_Diet.cds, reduction_method = "UMAP")

plot_cells(cds = E11.5_myogenic_Diet.cds, label_cell_groups = FALSE, color_cells_by = "ident", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = FALSE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.65))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.65))+theme(axis.ticks.y = element_line(size = 0.65))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))

plot_cells(cds = E11.5_myogenic_Diet.cds, label_cell_groups = FALSE, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = TRUE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.65))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.65))+theme(axis.ticks.y = element_line(size = 0.65))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))
```
<h3>1st</h3>

```r
E11.5_myogenic_Diet.cds <- estimate_size_factors(E11.5_myogenic_Diet.cds)
E11.5_myogenic_Diet.cds@rowRanges@elementMetadata@listData[['gene_short_name']] <- rownames(E11.5_myogenic_Diet.cds[['RNA']])

RP_genes.module <- c("Prrx1","Pdgfra","Pdgfrb","Col1a1","Col3a1")

rowData(E11.5_myogenic_Diet.cds)$gene_name <- rownames(E11.5_myogenic_Diet.cds)
rowData(E11.5_myogenic_Diet.cds)$gene_short_name <- rowData(E11.5_myogenic_Diet.cds)$gene_name

RP_genes.module_cds <- E11.5_myogenic_Diet.cds[rowData(E11.5_myogenic_Diet.cds)$gene_short_name %in% RP_genes.module,label_by_short_name = FALSE]


plot_genes_in_pseudotime(RP_genes.module_cds,color_cells_by="ident")+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+ theme(strip.text = element_text(size=8, face="bold.italic"))

saveRDS(E11.5_myogenic_Diet.cds,"Monocle3.E11.5_myogenic_Diet.KEEPumap.rds")
```

<h3>3rd</h3>

```r


E11.5_myogenic_2.cds <- as.cell_data_set(E11.5_myogenic)

E11.5_myogenic_2.cds <- preprocess_cds(E11.5_myogenic_2.cds, num_dim = 96)
E11.5_myogenic_2.cds <- align_cds(E11.5_myogenic_2.cds)
E11.5_myogenic_2.cds <- reduce_dimension(E11.5_myogenic_2.cds, preprocess_method = "Aligned")

plot_cells(E11.5_myogenic_2.cds, label_groups_by_cluster=FALSE,  color_cells_by = "ident")


E11.5_myogenic_2.cds <- cluster_cells(cd=E11.5_myogenic_2.cds,reduction_method = "UMAP")
E11.5_myogenic_2.cds <- learn_graph(E11.5_myogenic_2.cds,use_partition = TRUE)

E11.5_myogenic_2.cds <- order_cells(E11.5_myogenic_2.cds, reduction_method = "UMAP")



plot_cells(cds = E11.5_myogenic_2.cds, label_cell_groups = FALSE, color_cells_by = "ident", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = FALSE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))
ggsave('F1.germ.ident.plot.tiff', dpi=300, width=3, height=2)



plot_cells(cds = E11.5_myogenic_2.cds, label_cell_groups = FALSE, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = FALSE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))
```
<h3>4th</h3>

```r
E11.5_myogenic_2.cds <- estimate_size_factors(E11.5_myogenic_2.cds)
E11.5_myogenic_2.cds@rowRanges@elementMetadata@listData[['gene_short_name']] <- rownames(E11.5_myogenic_2.cds[['RNA']])

RP_genes.module <- c("Pax7","Myog","Myf5","H19","Pax3")

rowData(E11.5_myogenic_2.cds)$gene_name <- rownames(E11.5_myogenic_2.cds)
rowData(E11.5_myogenic_2.cds)$gene_short_name <- rowData(E11.5_myogenic_2.cds)$gene_name

RP_genes.module_cds <- E11.5_myogenic_2.cds[rowData(E11.5_myogenic_2.cds)$gene_short_name %in% RP_genes.module,label_by_short_name = FALSE]


plot_genes_in_pseudotime(RP_genes.module_cds,color_cells_by="ident")+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+ theme(strip.text = element_text(size=8, face="bold.italic"))
```
<h3>5th</h3>

```r
RP_genes.module <- c("Mef2c")

rowData(E11.5_myogenic_2.cds)$gene_name <- rownames(E11.5_myogenic_2.cds)
rowData(E11.5_myogenic_2.cds)$gene_short_name <- rowData(E11.5_myogenic_2.cds)$gene_name

RP_genes.module_cds <- E11.5_myogenic_2.cds[rowData(E11.5_myogenic_2.cds)$gene_short_name %in% RP_genes.module,label_by_short_name = FALSE]


plot_genes_in_pseudotime(RP_genes.module_cds,color_cells_by="ident")+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+ theme(strip.text = element_text(size=8, face="bold.italic"))
```
