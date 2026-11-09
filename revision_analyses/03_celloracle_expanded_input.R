# 03_celloracle_expanded_input.R
#
# Reviewer 1, comment 4b: "FOXI1 overexpression could also promote other
# differentiation pathways than EHC (Fig. 5C). To exclude this possibility, the
# CellOracle analysis should be performed for the clusters 2, 6, and 11.
# Additionally, knock-out perturbation will clarify the necessity of FOXI1."
#
# Step 1 of 2 (local, R): export the CellOracle input for an expanded subset.
#   - The original 1,972 cells used for Fig. 5 (cluster 1 "otic", 19 "EHC",
#     12 "LHC"; cluster-1 cells beyond the 95th percentile of distance from the
#     cluster-1 UMAP median removed, as in celloracle/01_prepare_RNAlayer.ipynb)
#   - plus all cells of clusters 2, 6 and 11.
#   - Raw RNA counts for the same 3,000 HVGs (vst on all 13,551 cells, as in
#     celloracle/prepare_RNAlayer.R), WNN UMAP coordinates and labels.
# Step 2 (CellOracle 0.20.0): 04_celloracle_expanded_FOXI1_OE_KO.py
#
# Output: revision/intermediate/celloracle_expanded_input/
#   counts_hvg.mtx (genes x cells), genes.csv, cells.csv, meta.csv, umap.csv,
#   cell_selection_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source(file.path(Sys.getenv("TANG_PROJECT_DIR", "."),
                 "revision_analyses", "00_config.R"))

out_dir <- file.path(inter_dir, "celloracle_expanded_input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- load_rna_slim()
umap <- Embeddings(obj, reduction = "wnn.umap")
cl <- as.character(obj$wsnn_clusters)

# CellOracle GRN units. Names of the original three units are kept unchanged.
unit_of <- c("1" = "otic", "19" = "early_hair_cell", "12" = "late_hair_cell",
             "2" = "neural_prog_C2", "6" = "neural_prog_C6", "11" = "early_otic_prog_C11")

# Original cluster-1 outlier filter (numpy median / quantile 0.95 / strict "<")
otic <- which(cl == "1")
centre <- apply(umap[otic, ], 2, median)
d <- sqrt(rowSums(sweep(umap[otic, ], 2, centre)^2))
thr <- quantile(d, 0.95, type = 7)   # same as numpy's default "linear" method
otic_keep <- otic[d < thr]
original_cells <- colnames(obj)[c(otic_keep, which(cl %in% c("19", "12")))]
message("Original subset: ", length(original_cells), " cells (", length(otic_keep), " cluster-1 cells kept, ",
        length(otic) - length(otic_keep), " removed)")
stopifnot(length(original_cells) == 1972)

# Cross-check against the 1,972 cells of the CellRank/velocity subset, if available
vel_h5ad <- file.path(proj_dir, "RNAvelocity", "anndata", "intergrated_anndata_velocity_subset3.h5ad")
if (file.exists(vel_h5ad) && requireNamespace("rhdf5", quietly = TRUE)) {
  vel_cells <- as.character(rhdf5::h5read(vel_h5ad, "/obs/barcode"))  # full Seurat barcodes
  message("Overlap with velocity subset cells: ", length(intersect(vel_cells, original_cells)),
          " / ", length(vel_cells))
}

added_cells <- colnames(obj)[cl %in% c("2", "6", "11")]
cells <- c(original_cells, added_cells)
message("Expanded subset: ", length(cells), " cells")

hvg <- VariableFeatures(obj)
stopifnot(length(hvg) == 3000, all(c("FOXI1", "ATOH1") %in% hvg))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[hvg, cells]

meta <- data.frame(
  cell = cells,
  cluster_id = cl[match(cells, colnames(obj))],
  experiment_day = as.character(obj$experiment_day[match(cells, colnames(obj))]),
  nCount_RNA = obj$nCount_RNA[match(cells, colnames(obj))],
  nFeature_RNA = obj$nFeature_RNA[match(cells, colnames(obj))],
  in_original_fig5_subset = cells %in% original_cells
)
meta$cell_type <- unname(unit_of[meta$cluster_id])
meta$color <- unname(cluster_colors[meta$cluster_id])

Matrix::writeMM(counts, file.path(out_dir, "counts_hvg.mtx"))
write.csv(data.frame(gene = rownames(counts)), file.path(out_dir, "genes.csv"), row.names = FALSE)
write.csv(data.frame(cell = colnames(counts)), file.path(out_dir, "cells.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "meta.csv"), row.names = FALSE)
write.csv(data.frame(cell = cells, UMAP_1 = umap[cells, 1], UMAP_2 = umap[cells, 2]),
          file.path(out_dir, "umap.csv"), row.names = FALSE)

sel <- as.data.frame(table(cluster_id = meta$cluster_id, cell_type = meta$cell_type,
                           in_original_fig5_subset = meta$in_original_fig5_subset))
sel <- sel[sel$Freq > 0, ]
write.csv(sel, file.path(out_dir, "cell_selection_summary.csv"), row.names = FALSE)
print(sel)
message("Wrote ", out_dir)
