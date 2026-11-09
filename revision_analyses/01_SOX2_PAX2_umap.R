# 01_SOX2_PAX2_umap.R
#
# Reviewer 1, comment 2a: cells were FACS-sorted as SOX2-mGFP+, PAX2-mCherry+ or
# double-positive, but no UMAP showed these cohorts.
# Tang's planned response: "feature plots of expressions of SOX2 and PAX2 now are provided."
#
# The reporter transcripts (mEGFP, mCherry) are not in the GRCh38 reference, so the
# endogenous SOX2 and PAX2 transcripts are shown on the integrated WNN UMAP
# (all 13,551 cells; Fig. 1B embedding). RNA assay, LogNormalize (scale factor 1e4).
#
# Outputs
#   figures/R1-2a_SOX2_PAX2_featureplots.{pdf,png}   cluster map + SOX2 + PAX2
#   figures/R1-2a_SOX2_PAX2_codetection.{pdf,png}    SOX2+ only / PAX2+ only / double+
#   results/R1-2a_SOX2_PAX2_detection_by_cluster.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source(file.path(Sys.getenv("TANG_PROJECT_DIR", "."),
                 "revision_analyses", "00_config.R"))

obj <- load_rna_slim()
genes <- c("SOX2", "PAX2")
df <- umap_df(obj, genes)

# 1. Cluster map + SOX2 + PAX2 expression
p1 <- umap_cluster_plot(df, title = "WNN clusters (Fig. 1B)") |
  umap_feature_plot(df, "SOX2", label_clusters = TRUE) |
  umap_feature_plot(df, "PAX2", label_clusters = TRUE)
save_plot(p1, "R1-2a_SOX2_PAX2_featureplots", width = 15, height = 4.6)

# 2. Co-detection (>= 1 UMI). Single-nucleus RNA has dropouts, so these calls
#    under-estimate the reporter-positive fractions.
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes, rownames(df)]
df$sox2_pos <- as.numeric(cnt["SOX2", ]) > 0
df$pax2_pos <- as.numeric(cnt["PAX2", ]) > 0
lev <- c("Neither detected", "SOX2+ only", "PAX2+ only", "SOX2+ PAX2+")
df$codetect <- factor(ifelse(df$sox2_pos & df$pax2_pos, lev[4],
                      ifelse(df$sox2_pos, lev[2], ifelse(df$pax2_pos, lev[3], lev[1]))), levels = lev)
df_ord <- df[order(df$codetect), ]
cen <- cluster_centres(df)
p2 <- ggplot(df_ord, aes(UMAP_1, UMAP_2, colour = codetect)) +
  geom_point(size = 0.3, stroke = 0) +
  scale_colour_manual(values = c("grey88", "#e41a1c", "#4daf4a", "#984ea3"), name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggrepel::geom_text_repel(data = cen, aes(UMAP_1, UMAP_2, label = cluster), inherit.aes = FALSE,
                           size = 3, min.segment.length = Inf, seed = 1) +
  labs(title = "SOX2 / PAX2 transcript detection", x = "wnnUMAP_1", y = "wnnUMAP_2") +
  umap_base_theme()
save_plot(p2, "R1-2a_SOX2_PAX2_codetection", width = 7, height = 5.2)

# Combined supplementary figure (panels A-D), styled like manuscript Fig. 1B:
# no axes or panel titles, cluster numbers on the map, legend "number + cell type".
fig1b_names <- c("supporting cells", "otic prosensory progenitors (Otic)", "neural progenitors", "stroma cells",
                 "mesenchymal cells", "epithelial cells", "neural progenitors", "Schwann cells", "stroma cells",
                 "myogenic precursors", "epithelial cells", "early otic progenitors", "late hair cells (LHC)",
                 "stroma cells", "neural progenitor", "dark cells", "nerve", "melanocytes", "epithelial cells",
                 "early hair cells (EHC)")
legend_lab <- setNames(paste(0:19, fig1b_names), as.character(0:19))
txt <- 16
fig_theme <- theme_void(base_size = txt) +
  theme(legend.text = element_text(size = txt - 2), legend.title = element_text(size = txt - 1),
        legend.key.height = unit(0.55, "cm"), plot.margin = margin(5, 5, 5, 5))
cen_lab <- cluster_centres(df)
num_layer <- ggrepel::geom_text_repel(data = cen_lab, aes(UMAP_1, UMAP_2, label = cluster), inherit.aes = FALSE,
                                      size = 5, min.segment.length = Inf, seed = 1)
pt <- 0.6
xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
x0 <- xr[1] - 0.02 * diff(xr); y0 <- yr[1] - 0.02 * diff(yr)
len <- 0.12 * min(diff(xr), diff(yr))
umap_arrows <- list(
  annotate("segment", x = x0, y = y0, xend = x0 + len, yend = y0,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.6),
  annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + len,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.6),
  annotate("text", x = x0 + len / 2, y = y0 - 0.035 * diff(yr), label = "UMAP_1", size = 4),
  annotate("text", x = x0 - 0.035 * diff(xr), y = y0 + len / 2, label = "UMAP_2", size = 4, angle = 90)
)

fA <- ggplot(df, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = pt, stroke = 0) + num_layer + umap_arrows +
  scale_colour_manual(values = cluster_colors, labels = legend_lab, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1)) + fig_theme

expr_panel <- function(gene) {
  d <- df[order(df[[gene]]), ]
  ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[gene]])) +
    geom_point(size = pt, stroke = 0) + num_layer + umap_arrows +
    scale_colour_gradient(low = "lightgrey", high = "blue", name = bquote(italic(.(gene)))) +
    guides(colour = guide_colourbar(barheight = unit(4, "cm"))) + fig_theme
}
fB <- expr_panel("SOX2")
fC <- expr_panel("PAX2")
fD <- ggplot(df_ord, aes(UMAP_1, UMAP_2, colour = codetect)) +
  geom_point(size = pt, stroke = 0) + num_layer + umap_arrows +
  scale_colour_manual(values = c("grey88", "#e41a1c", "#4daf4a", "#984ea3"), drop = FALSE, name = NULL,
                      labels = c("neither detected", expression(italic(SOX2)^"+"*" only"), expression(italic(PAX2)^"+"*" only"),
                                 expression(italic(SOX2)^"+"*italic(PAX2)^"+"))) +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) + fig_theme

tag_theme <- theme(plot.tag = element_text(size = 22, face = "bold"))
fA <- fA + tag_theme; fB <- fB + tag_theme; fC <- fC + tag_theme; fD <- fD + tag_theme
p_comb <- (fA | fB) / (fC | fD) + plot_annotation(tag_levels = "a")
save_plot(p_comb, "R1-2a_FigSX_SOX2_PAX2_combined", width = 17, height = 13)

# 3. Detection table by cluster
summarise_detection <- function(d) {
  d %>% summarise(
    n_cells = n(),
    pct_SOX2_pos = round(100 * mean(sox2_pos), 1),
    pct_PAX2_pos = round(100 * mean(pax2_pos), 1),
    pct_SOX2_only = round(100 * mean(sox2_pos & !pax2_pos), 1),
    pct_PAX2_only = round(100 * mean(!sox2_pos & pax2_pos), 1),
    pct_double_pos = round(100 * mean(sox2_pos & pax2_pos), 1),
    pct_either_pos = round(100 * mean(sox2_pos | pax2_pos), 1),
    .groups = "drop")
}
tab <- df %>% group_by(cluster) %>% summarise_detection() %>%
  mutate(cluster = as.character(cluster), cluster_name = cluster_labels[cluster], .after = cluster)
tab_all <- summarise_detection(df) %>% mutate(cluster = "All", cluster_name = "All cells", .before = 1)
tab <- bind_rows(tab, tab_all)
write.csv(tab, file.path(res_dir, "R1-2a_SOX2_PAX2_detection_by_cluster.csv"), row.names = FALSE)
print(as.data.frame(tab))
