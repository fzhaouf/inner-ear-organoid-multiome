#!/usr/bin/env python
"""05_celloracle_shift_relative_to_neighbourhood.py

Follow-up to 04b (R1-4b). Cluster 19 (EHC) is tiny, so few cluster-1 cells border it,
while clusters 2 and 11 share long borders with cluster 1. Raw mean shifts over all
cluster-1 cells therefore favour large neighbouring clusters. This script re-expresses the
FOXI1 OE / KO shifts of cluster-1 cells relative to each destination's share of their
UMAP neighbourhood.

The neighbourhood is CellOracle's: 201 nearest UMAP neighbours, self excluded
(estimate_transition_prob with n_neighbors=200, sampled_fraction=1).

Input : results/R1-4b_celloracle/R1-4b_per_cell_shift.csv.gz (from 04b)
Output: results/R1-4b_celloracle/R1-4b_C1_shift_relative_to_neighbourhood.csv
Run   : python 05_celloracle_shift_relative_to_neighbourhood.py [--rev <revision dir>]
"""
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

parser = argparse.ArgumentParser()
parser.add_argument("--rev", default="./revision")
args = parser.parse_args()
res = os.path.join(args.rev, "results", "R1-4b_celloracle")

pc = pd.read_csv(os.path.join(res, "R1-4b_per_cell_shift.csv.gz"), index_col=0, dtype={"cluster_id": str})
emb = pc[["UMAP_1", "UMAP_2"]].values
K = NearestNeighbors(n_neighbors=201).fit(emb).kneighbors_graph(mode="connectivity")
cl = pc["cluster_id"].values
order = ["1", "19", "12", "2", "6", "11"]
name = {"1": "C1 otic prosensory", "19": "C19 EHC", "12": "C12 LHC",
        "2": "C2 neural prog.", "6": "C6 neural prog.", "11": "C11 early otic prog."}
M = np.column_stack([(cl == c).astype(float) for c in order])
base = (K @ M) / 201.0           # each cell's neighbour share per cluster

rows = []
src = cl == "1"
for j, c in enumerate(order):
    if c == "1":
        continue
    b = base[src, j]
    exposed = b > 0
    for pert in ["OE", "KO"]:
        d = pc.loc[src, f"{pert}_shift_to_{c}"].values - pc.loc[src, f"{pert}_shift_to_{c}_random"].values
        rows.append({
            "perturbation": f"FOXI1_{pert}", "source": name["1"], "destination": name[c],
            "pct_C1_cells_with_neighbours_in_destination": round(100 * exposed.mean(), 1),
            "mean_baseline_neighbour_share_pct": round(100 * b.mean(), 2),
            "mean_shift_obs_minus_random_pp_all_C1": round(100 * d.mean(), 3),
            "mean_shift_obs_minus_random_pp_exposed_C1": round(100 * d[exposed].mean(), 2) if exposed.any() else np.nan,
            "relative_shift_pct_of_baseline": round(100 * d.mean() / b.mean(), 1) if b.mean() > 0 else np.nan,
        })
out = pd.DataFrame(rows)
out.to_csv(os.path.join(res, "R1-4b_C1_shift_relative_to_neighbourhood.csv"), index=False)
print(out.to_string(index=False))
