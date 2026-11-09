# 03_chromatin_motif_chromVAR.R
#
# Chromatin landscape (Fig. 2):
#   a  genomic annotation of cell-type DARs (ChIPseeker)
#   b  TF motif enrichment in DARs (JASPAR2020, Signac FindMotifs)
#   c-d per-cell TF motif activity (chromVAR)
# Input: integrated object from 01_preprocessing_integration.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2020)
  library(TFBSTools)
  library(chromVAR)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(dplyr)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
obj <- readRDS(file.path(project_dir, "processed_data", "integrated_multiome_seurat.rds"))

# Cell-type grouping used for the Fig. 2 DAR / motif analysis
fig2_types <- c(
  "0" = "Supporting cells", "1" = "Otic prosensory progenitors", "2" = "Neural progenitor",
  "3" = "Fibrocytes", "4" = "Mesenchymal cells", "5" = "Endothelial cells", "6" = "Glia",
  "7" = "Schwann cells", "8" = "Fibrocytes", "9" = "Cycling cells", "10" = "Glia",
  "11" = "Early otic progenitors", "12" = "Late hair cells", "13" = "Fibroblasts",
  "14" = "Endolymphatic membrane", "15" = "Dark cells", "16" = "Nerve", "17" = "Melanocytes",
  "18" = "Epithelial cells", "19" = "Early hair cells"
)
obj$cell_type <- unname(fig2_types[as.character(obj$wsnn_clusters)])
DefaultAssay(obj) <- "ATAC"
Idents(obj) <- "cell_type"

# ---------------------------------------------------------------------------
# Cell-type DARs (one vs rest; Wilcoxon; adjusted P < 0.05); types with >= 20 cells
# ---------------------------------------------------------------------------
keep <- names(which(table(obj$cell_type) >= 20))
sub <- subset(obj, cell_type %in% keep)
dars <- lapply(setNames(keep, keep), function(ct) {
  res <- FindMarkers(sub, ident.1 = ct, assay = "ATAC", test.use = "wilcox",
                     min.pct = 0.05, logfc.threshold = 0.1, only.pos = TRUE)
  res$peak <- rownames(res)
  res[res$p_val_adj < 0.05, ]
})

# a. Genomic annotation of DARs (promoter = TSS +/- 3 kb)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annot <- bind_rows(lapply(names(dars), function(ct) {
  if (nrow(dars[[ct]]) < 10) return(NULL)
  ann <- annotatePeak(StringToGRanges(dars[[ct]]$peak), tssRegion = c(-3000, 3000),
                      TxDb = txdb, level = "gene", verbose = FALSE)
  data.frame(ann@annoStat, cell_type = ct)
}))
write.csv(annot, "DAR_genomic_annotation.csv", row.names = FALSE)

# b. Motif enrichment in DARs (JASPAR2020 human CORE; hypergeometric test; BH)
pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
sub <- AddMotifs(sub, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motifs <- bind_rows(lapply(names(dars), function(ct) {
  if (nrow(dars[[ct]]) < 30) return(NULL)
  data.frame(FindMotifs(sub, features = dars[[ct]]$peak), cell_type = ct)
}))
motifs$padj <- p.adjust(motifs$pvalue, method = "BH")
write.csv(motifs, "DAR_motif_enrichment.csv", row.names = FALSE)

# c-d. chromVAR motif activity on all cells, plotted on the WNN UMAP
obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
obj <- RunChromVAR(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(obj) <- "chromvar"
ids <- ConvertMotifID(obj, name = c("ATOH1(var.2)", "POU4F3"), assay = "ATAC")
p <- FeaturePlot(obj, features = ids, reduction = "wnn.umap", min.cutoff = "q5", max.cutoff = "q95")
ggplot2::ggsave("chromVAR_ATOH1_POU4F3.pdf", p, width = 10, height = 5)

saveRDS(obj, file.path(project_dir, "processed_data", "integrated_multiome_seurat_chromvar.rds"))
