# Single-cell multiome analysis of hair cell development in human inner ear organoids

Analysis code for:
*Single-cell multiome atlas unveils a driver gene during the development of human vestibular hair cells*
(Zhao F, Dykxhoorn D, Liu X-Z, Chen L, Tang P-C).

The scripts document the main analyses and parameters used in the study. File paths are relative to a
project directory (`PROJECT_DIR` / `TANG_PROJECT_DIR`, default `.`) and may need to be adapted.

## Data

Raw and processed single-cell multiome data (days 45, 50 and 55) are available from GEO under accession
**GSE325916**. Cell Ranger ARC (v2.0.2) with the GRCh38-2020-A-2.0.0 reference was used for alignment
and counting; velocyto (v0.17.17, `run10x`) was used to obtain spliced/unspliced counts.

## Code

| Script | Analysis | Figures |
|---|---|---|
| `code/01_preprocessing_integration.R` | QC, unified peak set, SCTransform, rPCA integration, ATAC LSI integration, WNN, UMAP, clustering | Fig. 1 |
| `code/02_annotation_DE_DA.R` | Cluster annotation, marker genes, DE (EHC vs clusters 1, 4, 8), DA peaks, coverage plots | Fig. 1, Fig. 4, Fig. S2, Fig. S4 |
| `code/03_chromatin_motif_chromVAR.R` | DAR genomic annotation, JASPAR2020 motif enrichment, chromVAR | Fig. 2 |
| `code/04_export_for_velocity.R` | Export of the integrated object for Python | – |
| `code/05_rna_velocity_cellrank.py` | scVelo RNA velocity, CellRank fate probabilities, lineage drivers, diffusion pseudotime, gene trends | Fig. 3, Fig. S3N |
| `code/06_celloracle_input.R` | CellOracle RNA input and Cicero co-accessibility | Fig. 5 |
| `code/07_celloracle.py` | Base GRN, cluster-specific GRNs, network centrality, FOXI1 overexpression simulation | Fig. 5 |

Analyses added during revision (`revision_analyses/`):

| Script | Analysis |
|---|---|
| `00_config.R`, `00_prepare_rna_object.R` | Shared settings; RNA object used by the scripts below |
| `01_SOX2_PAX2_umap.R` | *SOX2* and *PAX2* expression on the integrated UMAP |
| `02_FOXI1_whole_dataset.R` | *FOXI1* expression across all cells and by sampling day |
| `03_celloracle_expanded_input.R` | CellOracle input including clusters 2, 6 and 11 |
| `04_celloracle_expanded_FOXI1_OE_KO.py` | CellOracle *FOXI1* overexpression and knock-out with clusters 2, 6 and 11 |
| `05_celloracle_shift_relative_to_neighbourhood.py` | Cluster-level quantification of simulated cell-state shifts |

Run the scripts in numeric order within each folder.

## Software

- R 4.2.2: Seurat 5.1.0, Signac 1.13.0, chromVAR 1.20.2, JASPAR2020 0.99.10, biomaRt 2.54.0,
  EnsDb.Hsapiens.v86 2.99.0, ChIPseeker 1.34.1, monocle3 1.3.1, cicero 1.3.9
- Python 3.9.18: scVelo 0.3.3, CellRank 2.0.7, Scanpy 1.9.8
- Python 3.10.19: CellOracle 0.20.0, GimmeMotifs 0.17.2, BEDTools 2.30.0

Full list: `environment/software_versions.csv`.
