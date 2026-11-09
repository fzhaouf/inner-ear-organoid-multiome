# 06_celloracle_input.R
#
# Inputs for CellOracle (Fig. 5):
#   1. RNA layer: raw counts of the hair-cell lineage (clusters 1, 19, 12), 3,000 HVGs
#      (vst on all cells), WNN UMAP -> exported as matrix/CSV files
#   2. ATAC layer: Cicero peak co-accessibility computed on all cells (base GRN)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(monocle3)
  library(cicero)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
out_dir <- file.path(project_dir, "celloracle")
dir.create(file.path(out_dir, "rna_export"), showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(file.path(project_dir, "processed_data", "integrated_multiome_seurat.rds"))

# ---------------------------------------------------------------------------
# 1. RNA layer
# ---------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
unit <- c("1" = "otic", "19" = "early_hair_cell", "12" = "late_hair_cell")
obj$cell_type <- unname(unit[as.character(obj$wsnn_clusters)])
lin <- subset(obj, cells = colnames(obj)[!is.na(obj$cell_type)])

exp_dir <- file.path(out_dir, "rna_export")
writeMM(GetAssayData(lin, assay = "RNA", layer = "counts"), file.path(exp_dir, "counts.mtx"))
write.csv(data.frame(gene = rownames(lin)), file.path(exp_dir, "genes.csv"), row.names = FALSE)
write.csv(data.frame(gene = VariableFeatures(obj)), file.path(exp_dir, "hvg.csv"), row.names = FALSE)
write.csv(data.frame(cell = colnames(lin), cell_type = lin$cell_type, experiment_day = lin$experiment_day,
                     UMAP_1 = Embeddings(lin, "wnn.umap")[, 1], UMAP_2 = Embeddings(lin, "wnn.umap")[, 2]),
          file.path(exp_dir, "meta.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. ATAC layer: Cicero co-accessibility on all cells
# ---------------------------------------------------------------------------
atac <- GetAssayData(obj, assay = "ATAC", layer = "counts")
atac@x[atac@x > 0] <- 1
rownames(atac) <- gsub("-", "_", rownames(atac))
parts <- do.call(rbind, strsplit(rownames(atac), "_"))
peakinfo <- data.frame(site_name = rownames(atac), chr = parts[, 1],
                       bp1 = as.integer(parts[, 2]), bp2 = as.integer(parts[, 3]), row.names = rownames(atac))
cellinfo <- data.frame(cells = colnames(atac), row.names = colnames(atac))

cds <- new_cell_data_set(atac, cell_metadata = cellinfo, gene_metadata = peakinfo)
cds <- cds[Matrix::rowSums(exprs(cds)) != 0, ]
set.seed(2017)
cds <- detect_genes(cds)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, method = "LSI")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "LSI")
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)

genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[!grepl("_|random", names(genome))]
conns <- run_cicero(cicero_cds, genomic_coords = data.frame(V1 = names(genome), V2 = as.numeric(genome)))
write.csv(conns, file.path(out_dir, "cicero_connections.csv"))
write.csv(data.frame(x = rownames(atac)), file.path(out_dir, "all_peaks.csv"))
