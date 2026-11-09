"""07_celloracle.py

CellOracle GRN inference and FOXI1 in silico perturbation (Fig. 5).
Python 3.10.19, CellOracle 0.20.0 (requires bedtools and the hg38 genome via genomepy).
Inputs from 06_celloracle_input.R (celloracle/).
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.applications import Pseudotime_calculator, Gradient_calculator, Oracle_development_module

D = "celloracle"

# ---------------------------------------------------------------------------
# 1. Base GRN: Cicero-linked peaks (co-accessibility >= 0.8) annotated to TSSs,
#    scanned for TF motifs (gimme.vertebrate.v5.0, FPR 0.02, score >= 10)
# ---------------------------------------------------------------------------
peaks = pd.read_csv(f"{D}/all_peaks.csv", index_col=0).x.values
conns = pd.read_csv(f"{D}/cicero_connections.csv", index_col=0)
tss = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss, cicero_connections=conns)
peak = integrated[integrated.coaccess >= 0.8][["peak_id", "gene_short_name"]].reset_index(drop=True)
peak = ma.check_peak_format(peak, "hg38")

tfi = ma.TFinfo(peak_data_frame=peak, ref_genome="hg38")
tfi.scan(fpr=0.02, motifs=None)                    # default motifs: gimme.vertebrate.v5.0
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
base_GRN = tfi.to_dataframe()
base_GRN.to_parquet(f"{D}/base_GRN_dataframe.parquet")   # also used by revision_analyses/04_*

# ---------------------------------------------------------------------------
# 2. RNA data: hair-cell lineage (clusters 1, 19, 12); scattered cluster-1 cells
#    (> 95th percentile of distance from the cluster-1 UMAP median) removed
# ---------------------------------------------------------------------------
E = f"{D}/rna_export"
adata = sc.AnnData(X=scipy.io.mmread(f"{E}/counts.mtx").T.tocsr().astype(np.float32),
                   obs=pd.read_csv(f"{E}/meta.csv").set_index("cell"),
                   var=pd.DataFrame(index=pd.read_csv(f"{E}/genes.csv")["gene"].values))
adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].values

otic = (adata.obs["cell_type"] == "otic").values
d = np.linalg.norm(adata.obsm["X_umap"][otic] - np.median(adata.obsm["X_umap"][otic], axis=0), axis=1)
keep = ~otic
keep[np.where(otic)[0]] = d < np.quantile(d, 0.95)
adata = adata[keep].copy()                           # 1,972 cells

# Pseudotime for the development flow (root: otic progenitor cell), computed on
# log-normalized expression of all genes (Seurat LogNormalize equivalent)
norm = adata.copy()
sc.pp.normalize_total(norm, target_sum=1e4)
sc.pp.log1p(norm)
sc.pp.pca(norm, n_comps=50)
sc.pp.neighbors(norm, n_neighbors=30, use_rep="X_pca")
sc.tl.diffmap(norm)
pt = Pseudotime_calculator(adata=norm, obsm_key="X_umap", cluster_column_name="cell_type")
pt.set_lineage(lineage_dictionary={"Otic_to_Hair": ["otic", "early_hair_cell", "late_hair_cell"]})
pt.set_root_cells(root_cells={"Otic_to_Hair": "day45_GTAATAGCAGCAACCT-1"})
pt.get_pseudotime_per_each_lineage()
adata.obs["Pseudotime"] = pt.adata.obs["Pseudotime"]
adata = adata[:, adata.var_names.isin(pd.read_csv(f"{E}/hvg.csv")["gene"])].copy()   # 3,000 HVGs for CellOracle

# ---------------------------------------------------------------------------
# 3. Oracle object, KNN imputation, cluster-specific GRNs, network scores
# ---------------------------------------------------------------------------
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="cell_type", embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()
n_comps = min(int(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]), 50)
k = int(0.025 * oracle.adata.shape[0])               # 49
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=4)

links = oracle.get_links(cluster_name_for_GRN_unit="cell_type", alpha=10)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
links.get_network_score()                            # degree / betweenness / eigenvector centrality (Fig. 5d-f)
links.merged_score.to_csv("celloracle_network_scores.csv")

# ---------------------------------------------------------------------------
# 4. FOXI1 overexpression (value 0.6) and development flow (Fig. 5b-c)
# ---------------------------------------------------------------------------
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

oracle.simulate_shift(perturb_condition={"FOXI1": 0.6}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)
oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
oracle.calculate_mass_filter(min_mass=5.8)

gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=5.8)
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 3})
gradient.calculate_gradient()

dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)
dev.visualize_development_module_layout_0(s=5, scale_for_simulation=10, s_grid=50, scale_for_pseudotime=45, vm=0.02)
