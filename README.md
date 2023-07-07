# Codes for single cell rna-sequencing

</head>
<body>
<p>Change name of cluster</p>
<code>
Idents(E_11_13_combined_2) <- "RNA_snn_res.0.18"
E_11_13_combined_3 <- RenameIdents(E_11_13_combined_2, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")
DimPlot(E_11_13_combined_3, label = TRUE)
</code>
</body>
</html>
