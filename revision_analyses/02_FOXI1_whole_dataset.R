# 02_FOXI1_whole_dataset.R
#
# Reviewer 1, comment 3a: "The authors should show FOXI1 expression in a UMAP plot of
# the entire dataset and discuss temporal changes in FOXI1 expression."
# (Fig. 4C only showed FOXI1 in the otic -> EHC -> LHC subset; Reviewer 1 comment 34
# asks for the same whole-dataset view.)
#
# Uses all 13,551 cells on the integrated WNN UMAP (Fig. 1B embedding),
# RNA assay, LogNormalize (scale factor 1e4).
#
# Outputs
#   figures/R1-3a_FOXI1_featureplot_all_cells.{pdf,png}
#   figures/R1-3a_FOXI1_featureplot_by_day.{pdf,png}
#   figures/R1-3a_FOXI1_dotplot_by_cluster.{pdf,png}
#   results/R1-3a_FOXI1_by_cluster.csv
#   results/R1-3a_FOXI1_by_day.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source(file.path(Sys.getenv("TANG_PROJECT_DIR", "."),
                 "revision_analyses", "00_config.R"))

obj <- load_rna_slim()
Idents(obj) <- "wsnn_clusters"
gene <- "FOXI1"

df_u <- umap_df(obj, gene)

# 1. FOXI1 on the whole dataset (cluster numbers as in Fig. 1B).
#    FOXI1 is detected in very few cells, so expressing cells are drawn larger
#    (pt_size_pos) and on top; the right panel zooms into the otic -> hair-cell branch.
p_all <- umap_feature_plot(df_u, gene, title = "FOXI1, all cells (days 45-55)",
                           label_clusters = TRUE, pt_size_pos = 1.6)
hc <- df_u[df_u$cluster %in% c("12", "19"), ]
zx <- range(hc$UMAP_1) + c(-2, 2)
zy <- c(min(hc$UMAP_2) - 5, max(hc$UMAP_2) + 1)
cen_in <- subset(cluster_centres(df_u), UMAP_1 > zx[1] & UMAP_1 < zx[2] & UMAP_2 > zy[1] & UMAP_2 < zy[2])
p_zoom <- umap_feature_plot(df_u, gene, title = "Zoom: hair-cell branch (clusters 1, 19, 12)",
                            pt_size = 0.8, pt_size_pos = 2.5) +
  ggplot2::geom_text(data = cen_in, ggplot2::aes(UMAP_1, UMAP_2, label = cluster), inherit.aes = FALSE,
                     size = 3.5, nudge_x = 0.6) +
  ggplot2::coord_cartesian(xlim = zx, ylim = zy)
p_all2 <- (p_all + ggplot2::annotate("rect", xmin = zx[1], xmax = zx[2], ymin = zy[1], ymax = zy[2],
                                     fill = NA, colour = "black", linetype = "dashed", linewidth = 0.3)) |
  p_zoom
save_plot(p_all2, "R1-3a_FOXI1_featureplot_all_cells", width = 12, height = 5.2)

# 2. Same plot split by sampling day (shared colour scale)
day_n <- table(df_u$day)
df_u$day_lab <- factor(paste0(sub("day", "Day ", df_u$day), " (n = ", format(day_n[as.character(df_u$day)], big.mark = ","), ")"))
df_u$day_lab <- factor(df_u$day_lab, levels = unique(df_u$day_lab[order(df_u$day)]))
p_day <- umap_feature_plot(df_u, gene, title = "FOXI1 by sampling day", facet = "day_lab", pt_size_pos = 1.6)
save_plot(p_day, "R1-3a_FOXI1_featureplot_by_day", width = 14, height = 4.8)

# 3. Dot plot across all 20 clusters (size = % cells with >= 1 UMI; colour = mean log-norm. expression)
cnt_g <- GetAssayData(obj, assay = "RNA", layer = "counts")[gene, rownames(df_u)]
dot <- data.frame(cluster = df_u$cluster, pos = as.numeric(cnt_g) > 0, expr = df_u[[gene]]) %>%
  group_by(cluster) %>%
  summarise(pct = 100 * mean(pos), avg = mean(expr), .groups = "drop") %>%
  mutate(label = factor(paste0(cluster, " ", cluster_labels[as.character(cluster)]),
                        levels = rev(paste0(0:19, " ", cluster_labels[as.character(0:19)]))))
p_dot <- ggplot(dot, aes(x = "FOXI1", y = label)) +
  geom_point(aes(size = pct, colour = avg)) +
  scale_size_area(max_size = 8, name = "% expressing") +
  scale_colour_gradient(low = "lightgrey", high = "blue", name = "Mean log-norm.\nexpression") +
  labs(x = NULL, y = NULL, title = "FOXI1 by cluster") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic"), plot.title = element_text(face = "bold"))
save_plot(p_dot, "R1-3a_FOXI1_dotplot_by_cluster", width = 4.8, height = 5.8)

# Combined supplementary figure, styled like manuscript Fig. 1B / Fig. S(R1-2a):
# a = all cells (cluster numbers), b-d = day 45 / 50 / 55; shared colour scale;
# FOXI1-expressing cells drawn larger and on top (FOXI1 is detected in very few cells).
txt <- 16
fig_theme <- theme_void(base_size = txt) +
  theme(legend.text = element_text(size = txt - 2), legend.title = element_text(size = txt - 1),
        plot.margin = margin(5, 5, 5, 5), plot.tag = element_text(size = 22, face = "bold"))
xr <- range(df_u$UMAP_1); yr <- range(df_u$UMAP_2)
x0 <- xr[1] - 0.02 * diff(xr); y0 <- yr[1] - 0.02 * diff(yr); len <- 0.12 * min(diff(xr), diff(yr))
umap_arrows <- list(
  annotate("segment", x = x0, y = y0, xend = x0 + len, yend = y0,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.6),
  annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + len,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.6),
  annotate("text", x = x0 + len / 2, y = y0 - 0.035 * diff(yr), label = "UMAP_1", size = 4),
  annotate("text", x = x0 - 0.035 * diff(xr), y = y0 + len / 2, label = "UMAP_2", size = 4, angle = 90))
lims <- c(0, max(df_u[[gene]]))
cen_f <- cluster_centres(df_u)
foxi1_panel <- function(d, numbers = FALSE, day_label = NULL) {
  d <- d[order(d[[gene]]), ]
  d$.sz <- ifelse(d[[gene]] > 0, 2.2, 0.6)
  p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[gene]])) +
    geom_point(aes(size = .sz), stroke = 0) + scale_size_identity() +
    scale_colour_gradient(low = "lightgrey", high = "blue", limits = lims, name = bquote(italic(.(gene)))) +
    guides(colour = guide_colourbar(barheight = unit(4, "cm"))) +
    coord_cartesian(xlim = xr + c(-0.06, 0.02) * diff(xr), ylim = yr + c(-0.06, 0.02) * diff(yr)) +
    umap_arrows + fig_theme
  if (numbers) p <- p + ggrepel::geom_text_repel(data = cen_f, aes(UMAP_1, UMAP_2, label = cluster),
                                                 inherit.aes = FALSE, size = 5, min.segment.length = Inf, seed = 1)
  if (!is.null(day_label)) p <- p + annotate("text", x = xr[2], y = yr[2], label = day_label,
                                             hjust = 1, vjust = 1, size = 6)
  p
}
gA <- foxi1_panel(df_u, numbers = TRUE)
gB <- foxi1_panel(df_u[df_u$day == "day45", ], day_label = "Day 45")
gC <- foxi1_panel(df_u[df_u$day == "day50", ], day_label = "Day 50")
gD <- foxi1_panel(df_u[df_u$day == "day55", ], day_label = "Day 55")
p_fig <- (gA | gB) / (gC | gD) + plot_annotation(tag_levels = "a")
save_plot(p_fig, "R1-3a_FigSX_FOXI1_combined", width = 15, height = 12)

# 4. Tables. Positive = >= 1 FOXI1 UMI.
cnt  <- GetAssayData(obj, assay = "RNA", layer = "counts")
expr <- GetAssayData(obj, assay = "RNA", layer = "data")[gene, ]
df <- data.frame(cluster = as.character(obj$wsnn_clusters),
                 day = as.character(obj$experiment_day),
                 umi = as.numeric(cnt[gene, ]),
                 lib = as.numeric(Matrix::colSums(cnt)),
                 expr = as.numeric(expr))
n_pos_total <- sum(df$umi > 0)

summ <- function(d) {
  d %>% summarise(
    n_cells = n(),
    n_FOXI1_pos = sum(umi > 0),
    pct_FOXI1_pos = round(100 * mean(umi > 0), 2),
    mean_lognorm_expr = round(mean(expr), 4),
    pseudobulk_CPM = round(1e6 * sum(umi) / sum(lib), 2),
    .groups = "drop")
}

by_cluster <- df %>% group_by(cluster) %>% summ() %>%
  mutate(pct_of_all_FOXI1_pos_cells = round(100 * n_FOXI1_pos / n_pos_total, 1),
         cluster_name = cluster_labels[cluster], .after = cluster) %>%
  arrange(as.integer(cluster))
write.csv(by_cluster, file.path(res_dir, "R1-3a_FOXI1_by_cluster.csv"), row.names = FALSE)

group_of <- function(cl) {
  dplyr::case_when(cl == "19" ~ "EHC (cluster 19)",
                   cl == "12" ~ "LHC (cluster 12)",
                   cl %in% c("1", "11") ~ "Otic progenitors (clusters 1 and 11)",
                   TRUE ~ "All other clusters")
}
by_day <- bind_rows(
  df %>% group_by(day) %>% summ() %>% mutate(cell_group = "All cells", .before = 1),
  df %>% mutate(cell_group = group_of(cluster)) %>% group_by(cell_group, day) %>% summ()
) %>%
  mutate(cell_group = factor(cell_group, levels = c("All cells", "EHC (cluster 19)", "LHC (cluster 12)",
                                                    "Otic progenitors (clusters 1 and 11)", "All other clusters"))) %>%
  arrange(cell_group, day)
write.csv(by_day, file.path(res_dir, "R1-3a_FOXI1_by_day.csv"), row.names = FALSE)

print(by_cluster, n = Inf)
print(by_day, n = Inf)
