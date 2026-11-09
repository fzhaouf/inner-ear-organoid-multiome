"""05_rna_velocity_cellrank.py

RNA velocity (scVelo) and fate mapping (CellRank) (Fig. 3, Fig. S3N).
Python 3.9.18, scVelo 0.3.3, CellRank 2.0.7, Scanpy 1.9.8.

Inputs: anndata/ export from 04_export_for_velocity.R and velocyto loom files
        (D45.loom, D50.loom, D55.loom).
"""
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scvelo as scv
import cellrank as cr
from scipy import io

# ---------------------------------------------------------------------------
# 1. AnnData from the Seurat export, merged with spliced/unspliced counts
# ---------------------------------------------------------------------------
adata = anndata.AnnData(X=io.mmread("anndata/counts.mtx").transpose().tocsr())
adata.obs = pd.read_csv("anndata/metadata.csv")
adata.obs.index = adata.obs["barcode"]
adata.var.index = open("anndata/gene_names.csv").read().splitlines()
adata.obsm["X_pca"] = pd.read_csv("anndata/pca.csv").to_numpy()
adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy()   # WNN UMAP

looms = []
for day, f in [("day45", "D45.loom"), ("day50", "D50.loom"), ("day55", "D55.loom")]:
    ld = scv.read(f, cache=True)
    ld.obs.index = [f"{day}_{bc.split(':')[1].replace('x', '')}-1" for bc in ld.obs.index]
    ld.var_names_make_unique()
    looms.append(ld)
ldata = looms[0].concatenate(looms[1:])
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
adata = scv.utils.merge(adata, ldata)
adata.obs["wsnn_clusters"] = adata.obs["wsnn_clusters"].astype("category")

# ---------------------------------------------------------------------------
# 2. Velocity on all cells (stochastic model; Fig. 3a)
# ---------------------------------------------------------------------------
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)
labels = {0: "Support", 1: "Otic", 12: "LHC", 19: "EHC"}
adata.obs["cell_type"] = adata.obs["wsnn_clusters"].map(lambda x: labels.get(x, str(x)))
scv.pl.velocity_embedding_stream(adata, basis="umap", color="cell_type", density=3, smooth=0.5)

# ---------------------------------------------------------------------------
# 3. Hair-cell lineage subset (clusters 1, 19, 12); dynamical model
#    Cluster-1 cells beyond the 95th percentile of distance from the cluster-1
#    UMAP median are removed (scattered cells).
# ---------------------------------------------------------------------------
sub = adata[adata.obs["cell_type"].isin(["Otic", "EHC", "LHC"])].copy()
otic = (sub.obs["cell_type"] == "Otic").values
d = np.linalg.norm(sub.obsm["X_umap"][otic] - np.median(sub.obsm["X_umap"][otic], axis=0), axis=1)
keep = ~otic
keep[np.where(otic)[0]] = d < np.quantile(d, 0.95)
sub = sub[keep].copy()                                    # 1,972 cells

sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca")
scv.pp.filter_and_normalize(sub)
scv.pp.moments(sub)
scv.tl.recover_dynamics(sub)
scv.tl.velocity(sub, mode="dynamical")
scv.tl.velocity_graph(sub)

# ---------------------------------------------------------------------------
# 4. CellRank: kernels, GPCCA macrostates, fate probabilities, lineage drivers (Fig. 3b-e)
# ---------------------------------------------------------------------------
vk = cr.kernels.VelocityKernel(sub).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(sub).compute_transition_matrix()
K = 0.8 * vk + 0.2 * ck
K.compute_transition_matrix()

g = cr.estimators.GPCCA(K)
g.fit(cluster_key="cell_type", n_states=7)           # resolves to 8 macrostates
g.set_terminal_states(states=["EHC", "LHC_2", "Otic_5"])
g.set_initial_states(states=["Otic_1"])
g.compute_fate_probabilities()
g.compute_lineage_drivers(lineages=["EHC"], return_drivers=True).to_csv("lineage_drivers_EHC.csv")
g.compute_lineage_drivers(lineages=["LHC_2"], return_drivers=True).to_csv("lineage_drivers_LHC.csv")

# ---------------------------------------------------------------------------
# 5. Diffusion pseudotime (root: otic cell with the lowest first diffusion
#    component) and GAM gene trends weighted by EHC fate probability
#    (Fig. 3e heatmap; Fig. S3N)
# ---------------------------------------------------------------------------
sc.pp.pca(sub, n_comps=50)
sc.pp.neighbors(sub, n_neighbors=30, n_pcs=30)
sc.tl.diffmap(sub)
otic = (sub.obs["cell_type"] == "Otic").values
sub.uns["iroot"] = int(np.where(otic)[0][np.argmin(sub.obsm["X_diffmap"][otic, 0])])
sc.tl.dpt(sub)
sub.obs["pseudotime"] = sub.obs["dpt_pseudotime"]

model = cr.models.GAM(sub, n_knots=6)
cr.pl.gene_trends(sub, model=model, lineages=["EHC"], genes=["MYO7A", "FOXI1"], time_key="pseudotime",
                  same_plot=True, ncols=1, hide_cells=True, weight_threshold=(1e-4, 1e-4),
                  save="FigS3N_gene_trends.pdf")
drivers = g.compute_lineage_drivers(lineages="EHC")
cr.pl.heatmap(sub, model=model, lineages=["EHC"], genes=drivers.head(78).index, cluster_key="cell_type",
              show_fate_probabilities=True, time_key="pseudotime", show_all_genes=False,
              weight_threshold=(1e-3, 1e-3), save="Fig3e_EHC_driver_heatmap.pdf")
