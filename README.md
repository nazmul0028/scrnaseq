# Codes for single cell RNA sequencing
<h3>Make cluster</h3>

```r
object <- FindNeighbors(object = object, dims = 1:40)
object <-FindClusters(object = object,graph.name = "RNA_snn")
object <- RunUMAP(object = object,graph = "RNA_snn")
DimPlot(object = object,reduction = "umap", group.by = "orig.ident")
```
<h3>Modify cluster resulution on the basis of need</h3>

```r
E_11_13_combined_2 <-FindClusters(object = E_11_13_combined_1, graph.name = "RNA_snn", resolution = 0.18)
cluster <- DimPlot(object=E_11_13_combined_2, reduction="umap", label=TRUE)
cluster
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
DotPlot(E_11_13_combined_3, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, col.min=-2, col.max=2, scale.min = 0, scale.max = 100) +
    RotatedAxis()
```
