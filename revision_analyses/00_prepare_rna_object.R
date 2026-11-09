# 00_prepare_rna_object.R
#
# Builds a slim RNA-only Seurat object from the integrated multiome object so the
# revision scripts do not each reload the 2.7 GB file. Run once.
#
# Input : processed_data/integrated_multiome_seurat.rds from code/01 (13,551 cells; WNN clusters 0-19)
# Output: revision/intermediate/rna_slim_seurat.rds
#
# RNA normalisation and HVG selection follow celloracle/prepare_RNAlayer.R
# (LogNormalize, scale factor 1e4; vst, 3,000 HVGs computed on all cells).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source(file.path(Sys.getenv("TANG_PROJECT_DIR", "."),
                 "revision_analyses", "00_config.R"))

message("Loading integrated object ...")
obj <- readRDS(file.path(proj_dir, "processed_data", "integrated_multiome_seurat.rds"))  # from code/01
invisible(gc())

message("Assays: ", paste(Assays(obj), collapse = ", "),
        " | reductions: ", paste(Reductions(obj), collapse = ", "))
message("RNA assay class: ", class(obj[["RNA"]])[1])

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
keep_meta <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "experiment_day", "wsnn_clusters"),
                       colnames(obj@meta.data))
meta <- obj@meta.data[, keep_meta, drop = FALSE]
umap <- Embeddings(obj, reduction = "wnn.umap")
rm(obj)
invisible(gc())

slim <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
slim[["wnn.umap"]] <- CreateDimReducObject(embeddings = umap[colnames(slim), ], key = "wnnUMAP_", assay = "RNA")
slim$wsnn_clusters <- factor(as.character(slim$wsnn_clusters), levels = as.character(0:19))
slim$experiment_day <- factor(slim$experiment_day, levels = day_levels)
slim$cluster_label <- factor(
  paste0(levels(slim$wsnn_clusters), " ", cluster_labels[levels(slim$wsnn_clusters)])[slim$wsnn_clusters],
  levels = paste0(as.character(0:19), " ", cluster_labels[as.character(0:19)])
)

slim <- NormalizeData(slim, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
slim <- FindVariableFeatures(slim, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
hvg <- VariableFeatures(slim)
message("HVGs: ", length(hvg), " | FOXI1 in HVGs: ", "FOXI1" %in% hvg, " | ATOH1 in HVGs: ", "ATOH1" %in% hvg)

message("Cells: ", ncol(slim), " | genes: ", nrow(slim))
print(table(cluster = slim$wsnn_clusters, day = slim$experiment_day, useNA = "ifany"))

saveRDS(slim, rna_slim_rds)
message("Saved ", rna_slim_rds)
