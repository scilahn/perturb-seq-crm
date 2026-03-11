# CLAUDE.md — perturb-seq-crm

> **GitHub:** https://github.com/scilahn/perturb-seq-crm  
> **S3 bucket:** `s3://learn-perturb-seq`  
> **Owner:** Richard Ahn, PhD — Translational Bioinformatics / Computational Biology

---

## 1. Project Purpose

This repository supports a structured self-directed learning program for
Perturb-seq / pooled perturbation single-cell analysis, with a focus on
**cardiovascular (primary), renal (secondary), and metabolic/hepatic (tertiary)**
biology. The learning strategy uses the **scPerturb portal** as the primary entry
point before transitioning to raw GEO datasets for novel hypothesis testing.

The two parallel goals are:
1. Develop end-to-end fluency in Perturb-seq analysis (pertpy, scvi-tools, scanpy)
2. Test novel biological hypotheses not yet explored by the original dataset authors

Full learning plan: `docs/scPerturb_Learning_Outline.docx`

---

## 2. Repository Structure

```
perturb-seq-crm/
│
├── CLAUDE.md                          ← this file
├── README.md
├── .gitignore                         ← excludes data/ entirely; use S3
│
├── envs/
│   ├── environment.yml                ← conda environment spec
│   └── requirements.txt               ← pip requirements
│
├── data/                              ← GITIGNORED; all data lives in S3
│   ├── raw/                           ← symlinks or download scripts only
│   │   ├── scperturb/                 ← harmonized .h5ad from scperturb.org
│   │   ├── GSE181897/                 ← Replogle 2022 (raw GEO)
│   │   ├── GSE222980/                 ← Jin 2023 cardiac iPSC-CM
│   │   ├── GSE173559/                 ← Tian 2021 cardiac iPSC-CM
│   │   ├── GSE217246/                 ← Xu 2023 kidney organoids
│   │   ├── GSE223448/                 ← Dhindsa 2023 cardiometabolic
│   │   ├── GSE264667/                 ← Replogle 2024 essential genes
│   │   └── zenodo_8187772/            ← Tran 2023 primary hepatocytes
│   └── processed/                     ← QC-filtered .h5ad objects
│       ├── scperturb/
│       └── geo/
│
├── src/                               ← shared library; importable as `crm`
│   ├── __init__.py
│   ├── preprocessing/
│   │   ├── __init__.py
│   │   ├── qc.py                      ← cell/gene QC filters
│   │   ├── guide_assignment.py        ← guide RNA assignment utilities
│   │   └── normalization.py
│   ├── perturbation_analysis/
│   │   ├── __init__.py
│   │   ├── differential.py            ← pseudobulk DE wrappers
│   │   ├── signatures.py              ← signature scoring, GSEA wrappers
│   │   ├── embeddings.py              ← SCVI / CPA model training
│   │   ├── edistance.py               ← E-distance / E-test utilities
│   │   └── interactions.py            ← gene interaction inference
│   ├── gwas_crossref/
│   │   ├── __init__.py
│   │   ├── opentargets.py             ← Open Targets API queries
│   │   └── gwas_catalog.py            ← GWAS catalog locus mapping
│   └── visualization/
│       ├── __init__.py
│       ├── umap_plots.py
│       ├── volcano.py
│       └── heatmaps.py
│
├── notebooks/
│   ├── 00_pipeline_dev/               ← scratch, methods learning, prototyping
│   │   ├── 00a_anndata_basics.ipynb
│   │   ├── 00b_edistance_tutorial.ipynb
│   │   └── 00c_pertpy_quickstart.ipynb
│   ├── 01_replogle_2022/              ← Phase 2 primary learning dataset
│   │   ├── 01a_qc_preprocessing.ipynb
│   │   ├── 01b_dimred_visualization.ipynb
│   │   ├── 01c_differential_expression.ipynb
│   │   ├── 01d_signature_scoring.ipynb
│   │   └── 01e_cpa_gene_programs.ipynb
│   ├── 02_jin_2023_cardiac/           ← Phase 3 cardiac iPSC-CM TFs
│   │   ├── 02a_qc_de.ipynb
│   │   ├── 02b_tf_signatures.ipynb
│   │   └── 02c_gwas_crossref.ipynb
│   ├── 03_tian_2021_cardiac/          ← Phase 3 large-scale cardiac CRISPRi
│   │   ├── 03a_qc_de.ipynb
│   │   └── 03b_replication_jin.ipynb
│   ├── 04_xu_2023_renal/              ← Phase 3 kidney organoid CRISPRi
│   │   ├── 04a_qc_de.ipynb
│   │   └── 04b_ckd_gwas_crossref.ipynb
│   ├── 05_dhindsa_2023_cardiometabolic/
│   ├── 06_tran_2023_hepatocytes/
│   ├── 07_replogle_2024_essential/
│   └── cross_dataset/                 ← Phase 3 meta-analyses and comparisons
│       ├── cd_01_harmonized_de.ipynb
│       ├── cd_02_conserved_effects.ipynb
│       ├── cd_03_gwas_mapping.ipynb
│       └── cd_04_hypothesis_catalogue.ipynb
│
├── results/
│   ├── figures/                       ← publication-quality plots (.pdf, .svg)
│   ├── tables/                        ← DE results, hypothesis catalogue (.csv)
│   └── reports/                       ← HTML notebook exports
│
├── batch/                             ← AWS Batch job definitions
│   ├── job_definitions/
│   │   ├── scvi_training.json
│   │   ├── genome_scale_de.json
│   │   └── cross_dataset_meta.json
│   └── submit_jobs.py
│
└── docs/
    ├── scPerturb_Learning_Outline.docx ← full learning plan
    └── dataset_inventory.xlsx          ← prioritized dataset table
```

---

## 3. S3 Data Architecture

All data lives in `s3://learn-perturb-seq`. The `data/` directory in the repo
is **gitignored** and should never contain actual data files — only download
scripts and symlinks.

```
s3://learn-perturb-seq/
├── raw/
│   ├── scperturb/                     ← .h5ad files from scperturb.org / Zenodo
│   │   ├── ReplogleWeissman2022_rpe1.h5ad
│   │   ├── ReplogleWeissman2022_k562_essential.h5ad
│   │   ├── ReplogleWeissman2022_k562_gwps.h5ad
│   │   ├── TianKampmann2019.h5ad
│   │   ├── TianKampmann2021_CRISPRi.h5ad
│   │   ├── TianKampmann2021_CRISPRa.h5ad
│   │   ├── NormanWeissman2019.h5ad
│   │   ├── SrivatsanTrapnell2020_sciplex3.h5ad
│   │   └── SchiebingerLander2019.h5ad
│   └── geo/
│       ├── GSE181897/                 ← Replogle 2022 raw counts
│       ├── GSE222980/                 ← Jin 2023 cardiac iPSC-CM
│       ├── GSE173559/                 ← Tian 2021 cardiac iPSC-CM
│       ├── GSE217246/                 ← Xu 2023 kidney organoids
│       ├── GSE223448/                 ← Dhindsa 2023 cardiometabolic
│       └── GSE264667/                 ← Replogle 2024 essential genes
├── processed/
│   ├── scperturb/                     ← post-QC .h5ad ready for analysis
│   └── geo/                           ← post-guide-assignment .h5ad
└── results/
    ├── figures/
    ├── tables/
    └── model_checkpoints/             ← scVI / CPA trained model weights
```

### Syncing data

```bash
# Pull a specific processed dataset to local data/ for interactive work
aws s3 cp s3://learn-perturb-seq/processed/scperturb/ReplogleWeissman2022_rpe1.h5ad \
    data/processed/scperturb/

# Push results back to S3
aws s3 sync results/ s3://learn-perturb-seq/results/

# List all processed datasets
aws s3 ls s3://learn-perturb-seq/processed/ --recursive
```

---

## 4. Compute Environment

### Interactive analysis (EC2)

| Parameter | Value |
|-----------|-------|
| Instance type | `r7i.4xlarge` (128GB RAM, 16 vCPU) |
| OS | Ubuntu 22.04 LTS |
| Storage | 2TB gp3 EBS (local .h5ad cache) |
| Access | SSH → JupyterLab on port 8888, or Claude Code via SSH |
| Cost control | CloudWatch alarm: auto-stop after 30 min CPU < 5% |

Scale up to `r7i.8xlarge` (256GB) for genome-scale Replogle datasets.  
Use `g4dn.2xlarge` (NVIDIA T4, 32GB GPU RAM) for scVI / CPA model training.

### Batch jobs (heavy compute)

Heavy jobs (genome-scale DE, scVI training, cross-dataset meta-analysis) should
be submitted via AWS Batch rather than run interactively. See `batch/` for job
definitions. Use Spot Instance fleets for 60–80% cost reduction.

### Environment setup

```bash
# Clone repo
git clone https://github.com/scilahn/perturb-seq-crm.git
cd perturb-seq-crm

# Create conda environment
conda env create -f envs/environment.yml
conda activate perturb-seq-crm

# Install src as editable package
pip install -e .

# Verify key packages
python -c "import pertpy; import scvi; import scanpy; print('Environment OK')"
```

### Core package versions (pinned in environment.yml)

```
pertpy >= 0.7.0
scvi-tools >= 1.1.0
scanpy >= 1.10.0
anndata >= 0.10.0
decoupler >= 1.6.0
pydeseq2 >= 0.4.0
boto3 >= 1.34.0
```

---

## 5. Learning Phases & Current Status

The project follows five sequential phases derived from `docs/scPerturb_Learning_Outline.docx`.
Update the status markers below as phases are completed.

### Phase Overview

| Phase | Focus | Primary Datasets | Duration | Status |
|-------|-------|-----------------|----------|--------|
| **Phase 1** | Environment + AnnData fundamentals | scPerturb portal | 1–2 wks | ⬜ Not started |
| **Phase 2** | Core analysis methods (QC, DE, signatures, CPA) | ReplogleWeissman2022 | 3–4 wks | ⬜ Not started |
| **Phase 3** | Advanced methods + cross-dataset analysis | scPerturb: Jin, Tian, Xu, Replogle | 3–4 wks | ⬜ Not started |
| **Phase 4** | Transition to raw GEO; novel hypothesis testing | GEO: GSE222980, GSE173559, GSE217246 | 2–3 wks | ⬜ Not started |
| **Phase 5** | Expand to Tier 2–3; multi-omics; outputs | All Tier 1–3 GEO datasets | Ongoing | ⬜ Not started |

### Phase 1 — Environment & Data Foundations
**Goal:** Working compute environment, data access confirmed, AnnData fluency.

Key tasks:
- [ ] Launch EC2 `r7i.4xlarge`, configure JupyterLab + Claude Code SSH access
- [ ] Set up S3 bucket `s3://learn-perturb-seq` with folder structure above
- [ ] Initialize this repo with directory structure; push first commit
- [ ] Download 2–3 small scPerturb .h5ad datasets (< 5GB each) to validate pipeline
- [ ] Complete `notebooks/00_pipeline_dev/00a_anndata_basics.ipynb`
- [ ] Verify: can load `.h5ad`, inspect `adata.obs['perturbation']`, identify control cells

**Checkpoint:** EC2 running → JupyterLab accessible → scPerturb files in S3 → AnnData loaded correctly.

### Phase 2 — Core Perturb-seq Analysis Methods
**Primary dataset:** `ReplogleWeissman2022` (scPerturb-harmonized, start with `_rpe1` sub-experiment)

Key tasks:
- [ ] Cell + gene QC; assess cells-per-perturbation distribution; flag ambiguous guide assignments
- [ ] HVG selection, log-normalization, PCA (50 components), UMAP
- [ ] SCVI latent space + SCVI-corrected UMAP
- [ ] Pseudobulk DE: each perturbation vs. non-targeting controls (PyDESeq2)
- [ ] `pertpy.tl.differential_response()` — compare with pseudobulk results
- [ ] GSEA / ORA on DE gene lists: MSigDB hallmarks, cardiac gene sets
- [ ] Perturbation similarity matrix (cosine distance on DE profiles) + leiden clustering
- [ ] Train `scvi.external.CPA` model; inspect perturbation embeddings
- [ ] Apply CINEMA-OT (`pertpy.tl.cinemaot`); compare program membership with CPA
- [ ] Validate: pick 5–10 known essential genes; confirm DE results match expected biology

**Checkpoint:** Full pipeline notebook complete → CPA model trained → ≥5 cardiac/metabolic perturbations manually validated.

### Phase 3 — Advanced Methods & Cross-Dataset Analysis
**Datasets:** scPerturb versions of Jin 2023, Tian 2021 (cardiac), Xu 2023 (renal), plus Replogle

Key tasks:
- [ ] Apply Phase 2 pipeline to `TianKampmann2021` (cardiac-relevant lysosomal/ferroptosis screen)
- [ ] Apply to Jin 2023 scPerturb version (if available) or raw GEO (see Phase 4 note below)
- [ ] Apply to Xu 2023 (kidney organoids) — note: may require direct GEO download
- [ ] Cross-dataset harmonized DE: restrict to shared gene universe, compare effect directions
- [ ] Identify conserved vs. context-specific perturbation effects across cardiac datasets
- [ ] GWAS locus mapping: pull cardiac + renal + metabolic GWAS catalog entries; map to perturbed genes
- [ ] Gene interaction inference (Norman 2019 approach) on cardiac TF pairs
- [ ] In silico perturbation prediction with trained CPA model
- [ ] Assemble hypothesis catalogue: ≥10 testable novel hypotheses, ranked by evidence strength

**Checkpoint:** Cross-dataset notebook complete → GWAS cross-reference table built → hypothesis catalogue drafted → top 3 hypotheses selected for Phase 4.

### Phase 4 — Transition to Raw GEO Data
**Datasets:** GSE222980 (Jin), GSE173559 (Tian cardiac), GSE217246 (Xu renal)

Key tasks:
- [ ] Download raw count matrices from GEO; run guide RNA assignment (MAGeCK or Cell Ranger multi)
- [ ] Apply QC: min 20 cells per perturbation, doublet removal (scrublet)
- [ ] Re-run DE pipeline on raw data; assess concordance with scPerturb-harmonized results
- [ ] Port heavy compute jobs to AWS Batch (Spot fleet: `r7i` + `g4dn`)
- [ ] Execute top 3 hypotheses from Phase 3 catalogue with full statistical rigor
- [ ] Validate findings with orthogonal evidence (Open Targets, STRING, proteomics literature)

**Checkpoint:** Raw GEO data processed → concordance with scPerturb results documented → ≥3 novel hypotheses tested → results in `results/`.

### Phase 5 — Novel Hypothesis Testing & Output (Ongoing)
Key tasks:
- [ ] Extend to Tier 2–3 datasets: Dhindsa 2023 (GSE223448), Tran 2023 (Zenodo), Replogle 2024 (GSE264667)
- [ ] Multi-omics integration: cross-reference perturbation gene programs with Olink/SomaScan signatures
- [ ] Finalize `src/` library with full docstrings and unit tests
- [ ] Draft preprint / conference abstract from top findings

---

## 6. Key Datasets

### Tier 1 — Start Here (scPerturb-harmonized available)

| Rank | scPerturb Name | GEO Accession | Biology | Cell Type | Notes |
|------|---------------|---------------|---------|-----------|-------|
| 1 | `ReplogleWeissman2022_rpe1` | GSE181897 | Genome-scale; metabolic/mitochondrial | RPE1 | Best learning dataset; start here |
| 2 | `ReplogleWeissman2022_k562_gwps` | GSE181897 | Genome-scale; cardiac/metabolic genes | K562 | 2.5M cells; use Batch for full analysis |
| 3 | `TianKampmann2021_CRISPRi` | GSE124703 | Lysosomal/ferroptosis → cardiac stress | iPSC neurons | Pathway overlap with cardiac ischemia |
| 4 | `TianKampmann2019` | GSE152988 | CRISPRi; mitochondrial + stress | iPSC neurons | Smaller; good for pipeline testing |

### Tier 1 — GEO only (not yet in scPerturb)

| Rank | GEO Accession | Biology | Cell Type | Notes |
|------|---------------|---------|-----------|-------|
| 2 | GSE222980 | Cardiac TF networks, sarcomere assembly | iPSC-CMs | Primary cardiac Perturb-seq; top hypothesis target |
| 3 | GSE173559 | Cardiac hypertrophy, heart development | iPSC-CMs | Large-scale; high statistical power |
| 4 | GSE217246 | Renal tubular biology, podocyte function | iPSC-kidney organoids | Best available renal perturbation dataset |

### Tier 2 — High Value

| Rank | Accession | Biology | Access |
|------|-----------|---------|--------|
| 5 | GSE159721 | Cardiac enhancers (iPSC-CMs) | Open |
| 6 | GSE223448 | Cardiometabolic: LDL, TG, T2D genes | Open |
| 7 | Zenodo 8187772 | Primary hepatocytes; NAFLD-relevant | Open |
| 8 | GSE264667 | Metabolic essentials; OXPHOS; lipids | Open |
| 9 | Zenodo 7041849 | scPerturb compendium (harmonized) | Open |

---

## 7. Analysis Conventions

### AnnData field names (scPerturb-harmonized)

Always use these field names when writing portable code:

```python
adata.obs['perturbation']        # perturbed gene name; 'control' for non-targeting
adata.obs['perturbation_type']   # 'CRISPRi', 'CRISPRa', 'CRISPR-cas9', 'drug'
adata.obs['nperts']              # int: number of perturbations per cell
adata.obs['cell_line']           # cell line or cell type
adata.obs['cancer']              # bool: is this a cancer-derived cell?

# Control cell selection
controls = adata[adata.obs['perturbation'] == 'control']
# Alternative control labels across datasets
CONTROL_LABELS = ['control', 'non-targeting', 'safe-harbor', 'CTRL', 'NT']
```

### QC thresholds (defaults — override per dataset as needed)

```python
MIN_CELLS_PER_PERT  = 20      # hard minimum; 200 preferred for E-test reliability
MIN_UMIS_PER_CELL   = 1000    # scPerturb recommendation
MIN_GENES_PER_CELL  = 500
MAX_PCT_MITO        = 20      # (%) for iPSC-CM datasets; adjust for other cell types
N_HVG               = 2000    # highly variable genes for PCA/SCVI
N_PCA_COMPONENTS    = 50
CONTROL_FRAC_MIN    = 0.10    # flag dataset if control cells < 10%
CONTROL_FRAC_MAX    = 0.30
```

### Differential expression

Use **pseudobulk** (PyDESeq2) as the primary method. Use `pertpy.tl.differential_response()`
as a secondary confirmation. Never use naive Wilcoxon on Perturb-seq data without
acknowledging its limitations.

```python
import pertpy as pt
import scanpy as sc

# Pseudobulk DE (preferred)
# Aggregate by perturbation, then run PyDESeq2
# See src/perturbation_analysis/differential.py

# pertpy method (secondary)
pert_data = pt.data.norman_2019()  # example
de = pt.tl.differential_response(pert_data, target_col='perturbation')
```

### E-distance / E-test (perturbation effect quantification)

Use E-distance as the primary measure of perturbation effect size (per scPerturb paper).
Minimum thresholds: ≥200 cells per perturbation, ≥1000 UMIs per cell.

```python
from scperturb import edist
estats = edist(adata, obs_key='perturbation')
# Returns pairwise E-distances; compare each perturbation to 'control'
```

### Gene program inference models

| Model | Use case | Package |
|-------|----------|---------|
| SCVI | Batch correction, latent space | `scvi.model.SCVI` |
| CPA | Perturbation effect decomposition, in silico prediction | `scvi.external.CPA` |
| CINEMA-OT | Causal gene program inference (non-parametric) | `pertpy.tl.cinemaot` |
| GEARS | Gene interaction prediction with prior knowledge | `gears` |

### GWAS cross-reference workflow

1. Pull cardiac/renal/metabolic GWAS loci from GWAS Catalog API (`src/gwas_crossref/gwas_catalog.py`)
2. Map loci to genes via Open Targets Genetics (`src/gwas_crossref/opentargets.py`)
3. Intersect with perturbed gene list per dataset
4. Flag GWAS genes with **no existing perturbation data** — these are the highest-priority hypothesis targets

---

## 8. Notebook Conventions

- **Naming:** `NN_descriptive_name.ipynb` where `NN` is zero-padded (01, 02, …)
- **Header cell:** Every notebook must begin with a markdown cell containing: dataset name,
  accession, phase, date, and a 1–2 sentence objective statement
- **Final cell:** Always save a timestamped HTML export to `results/reports/`
- **Figures:** Save all figures to `results/figures/` with descriptive names; commit figure
  metadata (not the figure files) to git
- **No data in notebooks:** Never store raw or processed expression matrices in the notebook
  itself; always load from S3 / local `data/` cache
- **Checkpoints:** At phase checkpoints, tag the repo: `git tag phase1-complete`, etc.

---

## 9. Git Workflow

```bash
# Feature branches for each notebook / analysis unit
git checkout -b analysis/replogle-2022-de

# Commit frequently with descriptive messages
git commit -m "feat(replogle): add pseudobulk DE pipeline with PyDESeq2"
git commit -m "fix(qc): correct min_cells_per_pert threshold for genome-scale data"
git commit -m "docs(hypothesis): add Phase 3 hypothesis catalogue v0.1"

# Never commit data files — check .gitignore before staging
git status  # verify no .h5ad, .csv, or large files staged

# Tag phase completions
git tag phase1-complete -m "Phase 1 checkpoint: environment + AnnData basics complete"
git push origin --tags
```

### .gitignore essentials

```
data/
*.h5ad
*.h5
*.loom
*.zarr
*.pkl
*.pt               # PyTorch model weights
__pycache__/
.ipynb_checkpoints/
.env
.venv
*.egg-info/
```

---

## 10. Claude Code Usage Notes

When working with Claude Code in this repo:

- **Data access:** Always pull data from S3 to `data/` before analysis; never hardcode
  S3 paths inside notebooks (use `src/` utility functions that read from a config)
- **AWS CLI:** Ensure `aws configure` is complete with appropriate IAM credentials before
  starting any session that touches S3 or Batch
- **Memory-aware operations:** For Replogle genome-scale data (2.5M cells), load only the
  sub-experiment needed; use `h5py` to extract specific perturbations without loading the
  full object if memory is constrained
- **Batch job submission:** Heavy jobs (scVI training, genome-scale DE) should be submitted
  via `batch/submit_jobs.py` rather than run interactively on EC2
- **Hypothesis tracking:** When a new hypothesis emerges during analysis, immediately add it
  to `notebooks/cross_dataset/cd_04_hypothesis_catalogue.ipynb` with supporting evidence

---

## 11. Key References & Resources

| Resource | URL | Use |
|----------|-----|-----|
| scPerturb portal | scperturb.org | Download harmonized .h5ad datasets |
| scPerturb paper | doi:10.1038/s41592-023-02144-y | E-distance methods; dataset table |
| pertpy docs | pertpy.readthedocs.io | Primary analysis toolkit |
| scvi-tools docs | docs.scvi-tools.org | SCVI, CPA, CINEMA-OT |
| Norman lab GitHub | github.com/normanlabfiles | Published Replogle 2022 analysis code |
| GWAS Catalog | ebi.ac.uk/gwas | Cardiac/renal/metabolic loci |
| Open Targets Genetics | genetics.opentargets.org | Gene-disease prioritization |
| STRING | string-db.org | Protein interaction networks |
| MSigDB | gsea-msigdb.org/msigdb | Gene set collections for enrichment |
| Replogle 2022 data portal | gwps.wi.mit.edu | Raw processed data for genome-scale screen |

---

## 12. Contact

Richard Ahn, PhD  
sungho.richard@gmail.com | (858) 336-4846  
https://scholar.google.com/citations?user=zvqPwUcAAAAJ&hl=en
