# 02_annotation_DE_DA.R
#
# Cluster annotation, marker genes, differential expression (DE) and differential
# accessibility (DA) (Fig. 1c-e, Fig. 4, Fig. S2, Fig. S4).
# Input: integrated object from 01_preprocessing_integration.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
obj <- readRDS(file.path(project_dir, "processed_data", "integrated_multiome_seurat.rds"))
Idents(obj) <- "wsnn_clusters"

# ---------------------------------------------------------------------------
# Cluster annotation (Fig. 1b); clusters 1 = otic prosensory progenitors,
# 19 = early hair cells (EHC), 12 = late hair cells (LHC)
# ---------------------------------------------------------------------------
cluster_labels <- c(
  "0" = "supporting cells", "1" = "otic prosensory progenitors", "2" = "neural progenitors",
  "3" = "stroma cells", "4" = "mesenchymal cells", "5" = "epithelial cells", "6" = "neural progenitors",
  "7" = "Schwann cells", "8" = "stroma cells", "9" = "myogenic precursors", "10" = "epithelial cells",
  "11" = "early otic progenitors", "12" = "late hair cells", "13" = "stroma cells",
  "14" = "neural progenitor", "15" = "dark cells", "16" = "nerve", "17" = "melanocytes",
  "18" = "epithelial cells", "19" = "early hair cells"
)
obj$cluster_label <- unname(cluster_labels[as.character(obj$wsnn_clusters)])

# ---------------------------------------------------------------------------
# Cluster markers (Fig. 1c) and canonical markers (Fig. 1d)
# ---------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

dot_genes <- c("PAX2", "JAG1", "ATOH1", "POU4F3", "MYO7A", "FBXO2", "BRICD5", "OC90")
ggsave("marker_dotplot.pdf", DotPlot(obj, features = dot_genes, assay = "RNA") + RotatedAxis(),
       width = 8, height = 6)

# ---------------------------------------------------------------------------
# DE: EHC (19) vs otic prosensory progenitors (1), cluster 4 and cluster 8
# (Fig. 4b, Fig. S2). Wilcoxon test on the SCT assay; BH-adjusted P values.
# ---------------------------------------------------------------------------
obj <- PrepSCTFindMarkers(obj, assay = "SCT")
for (ref in c("1", "4", "8")) {
  de <- FindMarkers(obj, ident.1 = "19", ident.2 = ref, assay = "SCT",
                    min.pct = 0.1, logfc.threshold = 0.1)
  de$p_val_adj <- p.adjust(de$p_val, method = "BH")
  write.csv(de, sprintf("DE_EHC_vs_cluster%s.csv", ref))
}

# ---------------------------------------------------------------------------
# DA peaks: EHC (19) and LHC (12) vs otic prosensory progenitors (1).
# Logistic regression with sequencing depth (nCount_ATAC) as covariate.
# ---------------------------------------------------------------------------
DefaultAssay(obj) <- "ATAC"
for (grp in c("19", "12")) {
  da <- FindMarkers(obj, ident.1 = grp, ident.2 = "1", assay = "ATAC",
                    test.use = "LR", latent.vars = "nCount_ATAC",
                    min.pct = 0.1, logfc.threshold = 0.1, only.pos = TRUE)
  write.csv(da[da$p_val_adj <= 0.05, ], sprintf("DA_peaks_cluster%s_vs_cluster1.csv", grp))
}

# ---------------------------------------------------------------------------
# Chromatin accessibility around hair-cell regulators (Fig. 4d, Fig. S4)
# ---------------------------------------------------------------------------
for (g in c("FOXI1", "ATOH1", "POU4F3")) {
  p <- CoveragePlot(obj, region = g, assay = "ATAC", extend.upstream = 2500,
                    extend.downstream = 2500, peaks = TRUE, annotation = TRUE)
  ggsave(sprintf("coverage_%s.pdf", g), p, width = 9, height = 9)
}

saveRDS(obj, file.path(project_dir, "processed_data", "integrated_multiome_seurat_annotated.rds"))
