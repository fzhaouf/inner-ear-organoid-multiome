# Shared settings for the revision R scripts (sourced by 00, 01, 02, 04a).

proj_dir <- Sys.getenv(
  "TANG_PROJECT_DIR",
  unset = "."
)
rev_dir   <- file.path(proj_dir, "revision")
fig_dir   <- file.path(rev_dir, "figures")
res_dir   <- file.path(rev_dir, "results")
inter_dir <- file.path(rev_dir, "intermediate")
for (d in c(fig_dir, res_dir, inter_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Slim RNA-only object written by 00_prepare_rna_object.R
rna_slim_rds <- file.path(inter_dir, "rna_slim_seurat.rds")

# Cluster names as shown in manuscript Fig. 1B (WNN clusters, resolution 0.1)
cluster_labels <- c(
  "0"  = "Supporting cells",
  "1"  = "Otic prosensory progenitors",
  "2"  = "Neural progenitors",
  "3"  = "Stroma cells",
  "4"  = "Mesenchymal cells",
  "5"  = "Epithelial cells",
  "6"  = "Neural progenitors",
  "7"  = "Schwann cells",
  "8"  = "Stroma cells",
  "9"  = "Myogenic precursors",
  "10" = "Epithelial cells",
  "11" = "Early otic progenitors",
  "12" = "Late hair cells (LHC)",
  "13" = "Stroma cells",
  "14" = "Neural progenitors",
  "15" = "Dark cells",
  "16" = "Nerve",
  "17" = "Melanocytes",
  "18" = "Epithelial cells",
  "19" = "Early hair cells (EHC)"
)

# Same colours as Fig. 1B (scanpy palette stored in
# RNAvelocity/anndata/intergrated_anndata_velocity.h5ad: uns/wsnn_clusters_colors)
cluster_colors <- setNames(
  c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6"),
  as.character(0:19)
)

day_levels <- c("day45", "day50", "day55")

load_rna_slim <- function() {
  if (!file.exists(rna_slim_rds)) stop("Run 00_prepare_rna_object.R first: ", rna_slim_rds)
  readRDS(rna_slim_rds)
}

# Save a ggplot/patchwork object as both PDF and PNG
save_plot <- function(p, name, width, height, dpi = 300) {
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white")
}

# ---------------------------------------------------------------------------
# UMAP plotting helpers (plain ggplot2).
# Seurat 5.1's FeaturePlot() calls patchwork's `&`, which fails with the
# installed ggplot2 4.0.3 / patchwork 1.3.0.9000, so feature plots are drawn here.
# Colours follow Seurat defaults (lightgrey -> blue).
# ---------------------------------------------------------------------------
umap_df <- function(obj, genes = NULL) {
  emb <- Seurat::Embeddings(obj, reduction = "wnn.umap")
  df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                   cluster = obj$wsnn_clusters, day = obj$experiment_day,
                   row.names = colnames(obj))
  if (length(genes)) {
    ex <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
    for (g in genes) df[[g]] <- as.numeric(ex[g, rownames(df)])
  }
  df
}

umap_base_theme <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10, face = "bold"))
}

cluster_centres <- function(df) {
  stats::aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = df, FUN = stats::median)
}

# Expression on UMAP; expressing cells drawn on top (Seurat order = TRUE).
# pt_size_pos: optional larger point size for expressing cells (for rare transcripts).
umap_feature_plot <- function(df, gene, title = gene, pt_size = 0.3, label_clusters = FALSE,
                              facet = NULL, cols = c("lightgrey", "blue"), pt_size_pos = NULL) {
  df <- df[order(df[[gene]]), ]
  df$.size <- if (is.null(pt_size_pos)) pt_size else ifelse(df[[gene]] > 0, pt_size_pos, pt_size)
  p <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2, colour = .data[[gene]])) +
    ggplot2::geom_point(ggplot2::aes(size = .size), stroke = 0) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_colour_gradient(low = cols[1], high = cols[2], name = "log-norm.\nexpression") +
    ggplot2::labs(title = title, x = "wnnUMAP_1", y = "wnnUMAP_2") +
    umap_base_theme()
  if (label_clusters) {
    cen <- cluster_centres(df)
    p <- p + ggrepel::geom_text_repel(data = cen, ggplot2::aes(UMAP_1, UMAP_2, label = cluster),
                                      inherit.aes = FALSE, size = 3, min.segment.length = Inf, seed = 1)
  }
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), nrow = 1)
  p
}

# Clusters on UMAP with Fig. 1B colours and cluster numbers.
umap_cluster_plot <- function(df, title = "WNN clusters", pt_size = 0.2) {
  cen <- cluster_centres(df)
  ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2, colour = cluster)) +
    ggplot2::geom_point(size = pt_size, stroke = 0) +
    ggplot2::scale_colour_manual(values = cluster_colors, guide = "none") +
    ggrepel::geom_text_repel(data = cen, ggplot2::aes(UMAP_1, UMAP_2, label = cluster),
                             inherit.aes = FALSE, size = 3.2, min.segment.length = Inf, seed = 1) +
    ggplot2::labs(title = title, x = "wnnUMAP_1", y = "wnnUMAP_2") +
    umap_base_theme()
}
