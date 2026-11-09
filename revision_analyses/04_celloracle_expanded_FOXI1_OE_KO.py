#!/usr/bin/env python
"""04_celloracle_expanded_FOXI1_OE_KO.py

Reviewer 1, comment 4b: "FOXI1 overexpression could also promote other
differentiation pathways than EHC (Fig. 5C). To exclude this possibility, the
CellOracle analysis should be performed for the clusters 2, 6, and 11.
Additionally, knock-out perturbation will clarify the necessity of FOXI1."

Step 2 of 2 (CellOracle 0.20.0, as in the original analysis).
Input from 03_celloracle_expanded_input.R: the original Fig. 5 cells
(cluster 1 otic, 19 EHC, 12 LHC; 1,972 cells) plus all cells of clusters 2, 6
and 11 (4,429 cells in total).

CellOracle steps and parameters follow celloracle/04_GRN_celloracle.ipynb and
05_perturbation_simulations.ipynb: same base GRN; PCA with the same elbow rule
(capped at 50 PCs); KNN imputation with k = 49 (value used for Fig. 5);
cluster-specific ridge GRNs (alpha = 10); simulation GRNs from
links.filter_links() defaults; FOXI1 OE = 0.6 and KO = 0; n_propagation = 3;
transition probabilities over 200 UMAP neighbours (knn_random, sampled_fraction = 1);
sigma_corr = 0.05; grid smooth = 0.8, n_grid = 40.

Quantification added for the reviewer's question (where do cells move?).
For each cell i, CellOracle's transition probabilities P_ij over its UMAP
neighbours are summed per destination cluster c and compared with uniform
neighbour weights (the baseline CellOracle itself subtracts when computing the
embedding shift):
    shift_i(c) = sum_{j in c} P_ij  -  (# neighbours of i in c) / (# neighbours of i)
shift_i(c) > 0 means the simulated perturbation moves cell i toward cluster c.
The same quantity from CellOracle's randomized simulation (permuted delta_X) is
the negative control. Per source cluster we report the mean shift (percentage
points) and a paired Wilcoxon test of observed vs randomized (BH-adjusted).

Usage:
    module load celloracle/0.20.0
    python 04_celloracle_expanded_FOXI1_OE_KO.py [--base .]
"""

import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import scipy.io
from scipy.stats import wilcoxon

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import scanpy as sc
import celloracle as co

parser = argparse.ArgumentParser()
parser.add_argument("--base", default=".")
parser.add_argument("--k", type=int, default=49, help="KNN imputation k (49 in the original Fig. 5 analysis)")
parser.add_argument("--oe_value", type=float, default=0.6, help="FOXI1 OE value (0.6 in the manuscript)")
parser.add_argument("--scale_oe", type=float, default=10, help="arrow scale for OE plots (original: 10)")
parser.add_argument("--scale_ko", type=float, default=25, help="arrow scale for KO plots (original: 25)")
parser.add_argument("--min_mass", type=float, default=None,
                    help="grid mass filter; default rescales the original 5.8 to the new grid size")
parser.add_argument("--n_jobs", type=int, default=8)
args = parser.parse_args()

IN_DIR = f"{args.base}/revision/intermediate/celloracle_expanded_input"
BASE_GRN = f"{args.base}/celloracle/base_GRN_dataframe.parquet"   # from code/07_celloracle.py
RES = f"{args.base}/revision/results/R1-4b_celloracle"
FIG = f"{args.base}/revision/figures/R1-4b_celloracle"
OBJ = f"{args.base}/revision/intermediate/celloracle_expanded_objects"
for d in (RES, FIG, OBJ):
    os.makedirs(d, exist_ok=True)

GOI = "FOXI1"
ORDER = ["1", "19", "12", "2", "6", "11"]   # source / destination order in tables and heatmaps
NAME = {"1": "C1 otic prosensory", "19": "C19 EHC", "12": "C12 LHC",
        "2": "C2 neural prog.", "6": "C6 neural prog.", "11": "C11 early otic prog."}
GROUPS = {"EHC+LHC (19+12)": ["19", "12"], "Other branches (2+6+11)": ["2", "6", "11"]}
plt.rcParams["savefig.dpi"] = 300


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def savefig(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def grid_sigma(emb, n_grid=40, smooth=0.8):
    """Gaussian kernel width CellOracle uses in calculate_p_mass (velocyto calculate_grid_arrows)."""
    steps = []
    for d in range(emb.shape[1]):
        m, M = emb[:, d].min(), emb[:, d].max()
        m = m - 0.025 * abs(M - m)
        M = M + 0.025 * abs(M - m)
        steps.append((M - m) / (n_grid - 1))
    return smooth * float(np.mean(steps))


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    o = np.argsort(q)
    r = q[o] * len(q) / (np.arange(len(q)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    tmp = np.empty_like(q)
    tmp[o] = np.minimum(r, 1)
    out[ok] = tmp
    return out


def cluster_legend(ax):
    handles = [Line2D([0], [0], marker="o", ls="", color=COLOR[c], label=NAME[c], markersize=6) for c in ORDER]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8, frameon=False)


# ---------------------------------------------------------------------------
# 1. AnnData from the R export
# ---------------------------------------------------------------------------
log("Loading input")
counts = scipy.io.mmread(f"{IN_DIR}/counts_hvg.mtx").T.tocsr().astype(np.float32)   # cells x genes
genes = pd.read_csv(f"{IN_DIR}/genes.csv")["gene"].astype(str).values
meta = pd.read_csv(f"{IN_DIR}/meta.csv", dtype={"cluster_id": str}).set_index("cell")
umap = pd.read_csv(f"{IN_DIR}/umap.csv").set_index("cell").loc[meta.index]

adata = sc.AnnData(X=counts, obs=meta.copy(), var=pd.DataFrame(index=genes))
unit = meta.groupby("cluster_id")["cell_type"].first().to_dict()
COLOR = meta.groupby("cluster_id")["color"].first().to_dict()
adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"], categories=[unit[c] for c in ORDER])
adata.uns["cell_type_colors"] = [COLOR[c] for c in ORDER]          # Fig. 1B colours
adata.obsm["X_umap"] = umap[["UMAP_1", "UMAP_2"]].values
adata.layers["raw_count"] = adata.X.copy()
n_by_cluster = adata.obs["cluster_id"].value_counts().reindex(ORDER).to_dict()
log(f"{adata.n_obs} cells x {adata.n_vars} genes; cells per cluster: {n_by_cluster}")

# ---------------------------------------------------------------------------
# 2. Oracle object, PCA, KNN imputation (as in 04_GRN_celloracle.ipynb)
# ---------------------------------------------------------------------------
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="cell_type", embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=pd.read_parquet(BASE_GRN))
oracle.perform_PCA()
n_comps = int(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0])
n_comps = min(n_comps, 50)
k = args.k
log(f"PCA dims = {n_comps}; KNN imputation k = {k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=args.n_jobs)

imp = sc.get.obs_df(oracle.adata, keys=[GOI], layer="imputed_count")[GOI]
cl = oracle.adata.obs["cluster_id"].astype(str).values
foxi1_tab = pd.DataFrame({
    "cluster_id": ORDER,
    "cluster": [NAME[c] for c in ORDER],
    "n_cells": [int((cl == c).sum()) for c in ORDER],
    "pct_cells_imputed_FOXI1_gt0": [round(100 * float((imp[cl == c] > 0).mean()), 2) for c in ORDER],
    "mean_imputed_FOXI1": [float(imp[cl == c].mean()) for c in ORDER],
    "max_imputed_FOXI1": [float(imp[cl == c].max()) for c in ORDER],
})
foxi1_tab.to_csv(f"{RES}/R1-4b_FOXI1_imputed_expression_by_cluster.csv", index=False)

# ---------------------------------------------------------------------------
# 3. Cluster-specific GRNs (as in 04/05 notebooks)
# ---------------------------------------------------------------------------
log("Inferring cluster-specific GRNs")
links = oracle.get_links(cluster_name_for_GRN_unit="cell_type", alpha=10, verbose_level=10, n_jobs=args.n_jobs)
links.to_hdf5(file_path=f"{OBJ}/expanded.celloracle.links")
links.filter_links()   # defaults (p = 0.001, coef_abs, top 10,000), as used for simulation in 05 notebook
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
oracle.to_hdf5(f"{OBJ}/expanded.celloracle.oracle")

# FOXI1 out-degree per cluster GRN (context for the simulations)
deg = []
for unit_name, df in links.filtered_links.items():
    deg.append({"cell_type": unit_name, "FOXI1_out_edges": int((df["source"] == GOI).sum()),
                "FOXI1_in_edges": int((df["target"] == GOI).sum())})
pd.DataFrame(deg).to_csv(f"{RES}/R1-4b_FOXI1_edges_per_cluster_GRN.csv", index=False)

# ---------------------------------------------------------------------------
# 4. Perturbation settings
# ---------------------------------------------------------------------------
# CellOracle accepts values up to 2*max - min of the imputed expression
upper = 2 * float(imp.max()) - float(imp.min())
oe_value = args.oe_value if args.oe_value <= upper else upper
if oe_value != args.oe_value:
    log(f"WARNING: OE value {args.oe_value} exceeds CellOracle's allowed range; using {oe_value:.3f}")

emb = oracle.adata.obsm["X_umap"]
in_orig = oracle.adata.obs["in_original_fig5_subset"].astype(str).str.upper().eq("TRUE").values
min_mass = args.min_mass if args.min_mass is not None else 5.8 * grid_sigma(emb) / grid_sigma(emb[in_orig])

# ---------------------------------------------------------------------------
# 5. Reference panels: clusters and imputed FOXI1
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=[13, 6])
oracle.plot_cluster_whole(ax=ax[0], s=5)
cluster_legend(ax[0])
ax[0].set_title("Expanded CellOracle subset (clusters 1, 19, 12 + 2, 6, 11)")
o = np.argsort(imp.values)
sca = ax[1].scatter(emb[o, 0], emb[o, 1], c=imp.values[o], cmap="viridis", s=5, rasterized=True)
plt.colorbar(sca, ax=ax[1], shrink=0.6, label="imputed FOXI1")
ax[1].set_title("FOXI1 (CellOracle imputed expression)")
ax[1].axis("off")
savefig(fig, "R1-4b_expanded_subset_clusters_and_FOXI1")

# ---------------------------------------------------------------------------
# 6. Simulations
# ---------------------------------------------------------------------------
def cluster_shift():
    """Per-cell net transition probability toward each destination cluster (observed, randomized)."""
    M = np.zeros((len(cl), len(ORDER)))
    for j, c in enumerate(ORDER):
        M[cl == c, j] = 1
    K = oracle.embedding_knn.toarray()
    K = K / K.sum(1, keepdims=True)
    base = K @ M
    obs = oracle.transition_prob @ M - base
    rnd = oracle.transition_prob_random @ M - base
    for g, members in GROUPS.items():
        idx = [ORDER.index(c) for c in members]
        obs = np.column_stack([obs, obs[:, idx].sum(1)])
        rnd = np.column_stack([rnd, rnd[:, idx].sum(1)])
    return obs, rnd


DEST = ORDER + list(GROUPS.keys())
summary_rows, per_cell = [], {}
mass_plot_done = False

for label, value, scale in [("OE", oe_value, args.scale_oe), ("KO", 0.0, args.scale_ko)]:
    log(f"Simulating FOXI1 {label} (value = {value:.3f})")
    oracle.simulate_shift(perturb_condition={GOI: value}, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
    if not mass_plot_done:
        oracle.suggest_mass_thresholds(n_suggestion=12)
        savefig(plt.gcf(), "QC_mass_threshold_suggestions")
        mass_plot_done = True
    oracle.calculate_mass_filter(min_mass=min_mass, plot=False)

    # vector field on clusters (as Fig. 5C) + randomized control
    fig, ax = plt.subplots(1, 2, figsize=[13, 6])
    for a, rnd_ctrl in zip(ax, [False, True]):
        oracle.plot_cluster_whole(ax=a, s=5)
        if rnd_ctrl:
            oracle.plot_simulation_flow_random_on_grid(scale=scale, ax=a, show_background=False)
        else:
            oracle.plot_simulation_flow_on_grid(scale=scale, ax=a, show_background=False)
    cluster_legend(ax[0])
    ax[0].set_title(f"FOXI1 {label}: simulated cell-state shift")
    ax[1].set_title("Randomized simulation (negative control)")
    savefig(fig, f"R1-4b_FOXI1_{label}_vector_field")

    # cluster-level shift
    obs, rnd = cluster_shift()
    delta = np.asarray(oracle.adata.layers["delta_X"])
    perturbed = np.abs(delta).sum(1) > 0
    per_cell[label] = (obs, rnd, perturbed, oracle.delta_embedding.copy())
    for s in ORDER:
        i = cl == s
        for j, d in enumerate(DEST):
            diff = obs[i, j] - rnd[i, j]
            try:
                p = wilcoxon(obs[i, j], rnd[i, j]).pvalue if np.any(diff != 0) else np.nan
            except ValueError:
                p = np.nan
            summary_rows.append({
                "perturbation": f"FOXI1_{label}", "perturb_value": value,
                "source_cluster": NAME[s], "destination": NAME.get(d, d),
                "n_source_cells": int(i.sum()), "n_source_cells_perturbed": int(perturbed[i].sum()),
                "mean_shift_observed_pp": 100 * float(obs[i, j].mean()),
                "mean_shift_randomized_pp": 100 * float(rnd[i, j].mean()),
                "mean_shift_obs_minus_random_pp": 100 * float(diff.mean()),
                "wilcoxon_p": p,
            })

summary = pd.DataFrame(summary_rows)
summary["wilcoxon_padj_BH"] = np.nan
for pert in summary["perturbation"].unique():
    m = summary["perturbation"] == pert
    summary.loc[m, "wilcoxon_padj_BH"] = bh(summary.loc[m, "wilcoxon_p"].values)
summary.to_csv(f"{RES}/R1-4b_cluster_transition_shift.csv", index=False)

# per-cell table (for custom plots)
pc = oracle.adata.obs[["cluster_id", "cell_type", "experiment_day"]].copy()
pc["UMAP_1"], pc["UMAP_2"], pc["FOXI1_imputed"] = emb[:, 0], emb[:, 1], imp.values
for label, (obs, rnd, perturbed, de) in per_cell.items():
    pc[f"{label}_perturbed"] = perturbed
    pc[f"{label}_dx"], pc[f"{label}_dy"] = de[:, 0], de[:, 1]
    for j, d in enumerate(DEST):
        tag = d if d in NAME else d.split(" ")[0]
        pc[f"{label}_shift_to_{tag}"] = obs[:, j]
        pc[f"{label}_shift_to_{tag}_random"] = rnd[:, j]
pc.to_csv(f"{RES}/R1-4b_per_cell_shift.csv.gz")

# heatmaps: observed minus randomized, source x destination cluster
fig, axes = plt.subplots(1, 2, figsize=[17, 5.6])
fig.subplots_adjust(wspace=0.25)
for ax, pert in zip(axes, ["FOXI1_OE", "FOXI1_KO"]):
    sub = summary[summary["perturbation"] == pert]
    mat = sub.pivot(index="source_cluster", columns="destination", values="mean_shift_obs_minus_random_pp")
    padj = sub.pivot(index="source_cluster", columns="destination", values="wilcoxon_padj_BH")
    rows = [NAME[c] for c in ORDER]
    cols = [NAME[c] for c in ORDER] + list(GROUPS.keys())
    mat, padj = mat.loc[rows, cols], padj.loc[rows, cols]
    vmax = max(float(np.nanmax(np.abs(mat.values))), 1e-6)
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for r in range(mat.shape[0]):
        for c in range(mat.shape[1]):
            v, q = mat.values[r, c], padj.values[r, c]
            star = "*" if (not np.isnan(q) and q < 0.05) else ""
            ax.text(c, r, f"{v:.2f}{star}", ha="center", va="center", fontsize=7,
                    color="white" if abs(v) > 0.6 * vmax else "black")
    ax.axvline(len(ORDER) - 0.5, color="black", lw=0.8)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows if pert == "FOXI1_OE" else [], fontsize=8)
    ax.set_xlabel("Destination cluster")
    if pert == "FOXI1_OE":
        ax.set_ylabel("Source cluster")
    ax.set_title(f"{pert.replace('_', ' ')}: net transition shift vs randomized (pp)")
    plt.colorbar(im, ax=ax, shrink=0.8, label="percentage points")
fig.text(0.01, -0.16, "* BH-adjusted paired Wilcoxon p < 0.05 (observed vs randomized simulation)", fontsize=8)
savefig(fig, "R1-4b_cluster_transition_shift_heatmap")

params = {
    "n_cells": int(adata.n_obs), "cells_per_cluster": {NAME[c]: int(n) for c, n in n_by_cluster.items()},
    "n_hvg": int(adata.n_vars), "n_pca_dims": n_comps, "knn_k": k,
    "grn_alpha": 10, "simulation_links_filter": "filter_links() defaults: p=0.001, coef_abs, top 10000",
    "FOXI1_OE_value": oe_value, "FOXI1_KO_value": 0.0, "n_propagation": 3,
    "transition_n_neighbors": 200, "sigma_corr": 0.05, "grid_smooth": 0.8, "n_grid": 40,
    "min_mass": min_mass, "arrow_scale_OE": args.scale_oe, "arrow_scale_KO": args.scale_ko,
    "max_imputed_FOXI1": float(imp.max()), "celloracle_version": co.__version__,
}
with open(f"{RES}/R1-4b_parameters.json", "w") as fh:
    json.dump(params, fh, indent=2)
log("Done")
