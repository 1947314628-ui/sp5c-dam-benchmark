# SP5C DAM-state benchmarking — companion code repository

Code and processed data for the manuscript:
**"Orofacial inflammatory pain is accompanied by a disease-associated microglial state transition in the trigeminal spinal nucleus caudalis"** (manuscript in preparation).

## Availability

- **Project name**: SP5C DAM-state benchmark — analysis code and processed data
- **Project home page**: https://github.com/1947314628-ui/sp5c-dam-benchmark
- **Archived version**: Zenodo DOI, minted on the `v1.0.0` release
- **Operating system(s)**: platform independent (analysis run on Linux; figures generated on Windows 11)
- **Programming language**: Python 3.11.2; R 4.6.0
- **Other requirements**: `environment/python_requirements.txt`, `environment/R_sessionInfo.txt`
- **License**: MIT (see `LICENSE`); data remain subject to the original GEO terms
- **Any restrictions to use by non-academics**: none

This repository reproduces every figure and table of the three benchmarks (scoring methods, cell–cell communication tools, integration methods) and the case-study analyses shown in the manuscript. All processed data files shipped here are plain CSV; raw data are public on GEO.

## Repository layout

```
data_r/                     Processed outputs of the cloud pipeline (all manuscript CSVs)
tables/                     Supplementary tables (communication scores, composition, etc.)
demo_case/                  R scripts that draw manuscript Figs 1–5 and Figs S2–S3
benchmarks/
  01_scoring/               Fig B1: DAM-scoring methods compared across 3 datasets
  02_communication/         Fig B2: CellChat canonical vs LIANA (CellPhoneDB / CellChat / consensus)
  03_integration/           Fig B3: harmony / uncorrected baseline / BBKNN / scVI
environment/
  R_sessionInfo.txt         Exact R session (Windows, R 4.6.0)
  python_requirements.txt   Python analysis pins (scanpy 1.11.5, liana 1.10.0, ...)
```

## Data provenance

- **GSE285766** — snRNA-seq (+ATAC) of mouse trigeminal Sp5C, CFA-induced orofacial pain model, 6 samples (CON/CFA/TB × 2). Download from GEO; the `data_r/` CSVs are the outputs of the cloud processing pipeline (QC → integration with Harmony → Leiden 38 clusters → annotation → marker/DAM-state export).
- **GSE172167 / GSE162807** — validation datasets (sciatic nerve injury; sham/SNI), used only through the cell-level DAM-score CSVs in `data_r/`.
- h5ad files are **not** shipped (~3.4 GB); they are regenerable from GEO with the pipeline described in `benchmarks/03_integration/README_methods.md`, or available from the authors.

## Environment

- **R**: 4.6.0 (Windows x64), full `sessionInfo()` in `environment/R_sessionInfo.txt`.
  Key: Seurat 5, CellChat, presto (compiled with Rtools45), ggplot2, patchwork.
- **Python (local)**: venv with the pins in `environment/python_requirements.txt`.
- **Python (cloud, benchmark 03)**: 192-core CPU server, Python 3.11.2, scanpy + harmonypy 2.0.0, bbknn 1.6.0, scvi-tools 1.4.2 (torch 2.13.0). No GPU.

## Benchmark 1 — DAM-scoring methods (Fig B1)

```bash
cd benchmarks/01_scoring
Rscript figB1_scoring_benchmark.R     # reads ../../data_r/*.csv
```

Compares the 21-gene net-score, the Keren-Shaul 11-gene subset score, and k-means thresholding across GSE285766 (870 microglia), GSE172167 and GSE162807.
Anchor findings: net-score ≡ cloud `dam_score` (ρ = 1); threshold-based classification is robust on 2 of 3 datasets; the net-score is **not** correlated with the KS-11 score (ρ = −0.03) because the steady-state module dominates its variance; state separation is driven by homeostatic-gene downregulation.

## Benchmark 2 — communication tools (Fig B2)

```bash
cd benchmarks/02_communication
Rscript install_presto.R                          # one-time, compiles presto
# 1) build the CFA-only Seurat object from raw GEO files (set SP5C_GEX to the
#    directory holding per-sample {features,matrix,barcodes}_gex files):
SP5C_GEX=/path/to/gex Rscript build_cfa_seurat.R
# 2) CellChat canonical (needs ~8 GB RAM; run steps 2 and 3 sequentially!)
Rscript run_cellchat.R
# 3) LIANA-py (cellphonedb / cellchat / rank_aggregate flavours, mouseconsensus):
python run_liana_py.py
# 4) figure + comparison tables:
Rscript figB2_comm_tools.R
```

Anchor findings: LIANA's CellChat flavour reproduces canonical CellChat exactly on all 243 shared keys (ρ = 1, log-log slope = 1, R² = 1 — implementation fidelity, not a bug); CellChat vs CellPhoneDB on the same 243 keys ρ = 0.96; same resource, different scoring function ρ = 0.579; top-30 overlap Jaccard = 0.25 across tools vs 0.88 within LIANA — **resource choice matters far more than the scoring function**. Only one microglia-involving pair (Apoe→Trem2, Astro→Microglia, LIANA-consensus) reaches any top-30.

## Benchmark 3 — integration methods (Fig B3)

```bash
# Cloud (G26, 192 cores, ~40 min; needs the annotated h5ad from the pipeline):
python run_integration_benchmark.py
# Local figure from the returned CSVs (already shipped in this directory):
Rscript figB3_integration.R                  # reads cells_umap_cluster.csv.gz + methods_summary.csv
```

| method | clusters | aRI cluster vs sample | weighted entropy | downstream runtime (s) |
|---|---|---|---|---|
| harmony (authoritative) | 38 | 0.0037 | 0.924 | — (reused 2026-08-31 embedding) |
| baseline (PCA, uncorrected) | 39 | 0.0033 | 0.922 | 122.4 |
| BBKNN | 33 | 0.0054 | 0.923 | 60.2 |
| scVI | 36 | 0.0124 | 0.923 | 61.2 |

Unified downstream: kNN k = 15, UMAP min_dist = 0.3, Leiden resolution = 0.8, seed 42.
Anchor finding: aRI vs sample ≈ 0 for every method — batch effect is negligible even uncorrected, the integration method does not change downstream conclusions, and the Harmony result is robust. DAM-like microglia colocalise in all four embedding spaces (Fig B3A).

## Demo case — manuscript figures

```bash
cd demo_case
Rscript fig1_overview.R      # Fig 1 (incl. panel G marker dot-plot)
Rscript fig2_dam.R           # Fig 2
Rscript fig3_communication.R # Fig 3 (A chord diagram + B/C)
Rscript fig4_pain_genes.R    # Fig 4
Rscript fig5_validation.R    # Fig 5
Rscript figS2_gse162807.R    # Fig S2
Rscript figS3_deg_shap_heatmap.R  # Fig S3
```

All scripts read `../data_r/` and `../tables/` only; figures are written into the working directory. `export_for_r.py` is the cloud-side CSV exporter kept for reference.

## Known caveats

- The cloud scanpy pipeline scripts (raw → QC → Harmony clustering → annotation) will be added under `pipeline/` in a subsequent commit; the method summary is already documented in `benchmarks/03_integration/README_methods.md`.
- The Harmony row in Benchmark 3 reuses the embedding of the authoritative 2026-08-31 run rather than recomputing it (documented in `README_methods.md`); its metrics are therefore listed without downstream runtime.
- Microglial positive-threshold choice and sex information for GSE285766 are noted as limitations in the manuscript.

## License

MIT. Data remain subject to the original GEO terms.
