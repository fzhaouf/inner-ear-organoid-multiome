# 04_export_for_velocity.R
#
# Exports the integrated object for scVelo/CellRank (05_rna_velocity_cellrank.py):
# RNA counts, cell metadata (incl. WNN clusters and WNN UMAP) and PCA.
# Spliced/unspliced counts were generated per sample with velocyto (run10x) on the
# Cell Ranger ARC outputs with the GRCh38-2020-A gene annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
out_dir <- file.path(project_dir, "anndata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(file.path(project_dir, "processed_data", "integrated_multiome_seurat.rds"))

obj$barcode <- colnames(obj)
obj$UMAP_1 <- Embeddings(obj, "wnn.umap")[, 1]
obj$UMAP_2 <- Embeddings(obj, "wnn.umap")[, 2]
write.csv(obj@meta.data, file.path(out_dir, "metadata.csv"), quote = FALSE, row.names = FALSE)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
writeMM(counts, file.path(out_dir, "counts.mtx"))
write.table(data.frame(gene = rownames(counts)), file.path(out_dir, "gene_names.csv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(Embeddings(obj, "pca"), file.path(out_dir, "pca.csv"), quote = FALSE, row.names = FALSE)
