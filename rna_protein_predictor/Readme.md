# Cross-assay RNA-protein disagreement in cancer

This repository contains a reproducible analysis of RNA expression, Human Protein Atlas (HPA) tumor immunohistochemistry, subcellular localization, and matched Clinical Proteomic Tumor Analysis Consortium (CPTAC) mass-spectrometry proteomics.

The project began as a model of protein detectability from bulk RNA expression. After rebuilding the evaluation and adding matched CPTAC data, the principal question became:

> Why are many proteins reported as absent or mostly absent by tumor IHC while being reproducibly quantified by mass spectrometry in the same cancer type?

The current results support systematic, compartment-dependent **cross-assay disagreement**. They do not establish that an IHC-negative protein is biologically absent, nor that every RNA-protein difference reflects post-transcriptional regulation.

## Project status

The download-to-results pipeline, seven-cancer HPA analysis, four-cancer CPTAC atlas, threshold sensitivity analysis, and gene-clustered localization models are implemented. The study is being prepared for a genetics conference and subsequent publication.

The chronological decision record is in [`docs/PROJECT_LOG.md`](docs/PROJECT_LOG.md). Planned analyses and reporting rules are in [`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md).

## Validated findings

### 1. Localization adds modest information beyond RNA

The audited master dataset contains 106,484 unique Ensembl gene-cancer rows, 15,212 genes, seven cancer groups, 40 localization features, and zero duplicated Ensembl-cancer keys.

Matched logistic-regression models were evaluated using seven complete leave-one-cancer-out (LOCO) folds.

| Feature set | Mean ROC-AUC | Mean PR-AUC | Brier score |
|---|---:|---:|---:|
| RNA | 0.7121 | 0.7185 | 0.2178 |
| RNA + gene features | 0.7128 | 0.7192 | 0.2177 |
| RNA + localization | 0.7272 | 0.7430 | 0.2129 |
| Full | 0.7276 | 0.7434 | 0.2128 |

In 2,000 gene-clustered paired bootstrap replicates, adding localization to RNA increased:

- ROC-AUC by 0.01527; 95% interval 0.01190-0.01873
- PR-AUC by 0.02486; 95% interval 0.02041-0.02940
- Empirical P < 0.0005; no sampled difference was non-positive

The increment persisted when restricted to localization-annotated genes. The effect is real but substantially smaller than the original proof-of-concept estimate.

### 2. Naive single-gene discordance lists do not survive rigorous testing

In CPTAC COAD, 96 tumors and 8,757 shared RNA/protein genes were analyzed. Patient-level high-RNA/low-protein events were assessed across BCM, Broad, and WashU RNA processing, with sensitivity adjustment for tumor purity, immune/stromal scores, and xCell composition.

No individual gene survived FDR correction. The number of raw P < 0.05 results was lower than the number expected by chance:

| RNA processing | Genes tested | Expected at P < 0.05 | Observed |
|---|---:|---:|---:|
| BCM | 6,956 | 348 | 177 |
| Broad | 6,858 | 343 | 155 |
| WashU | 6,648 | 332 | 238 |

These are processing robustness analyses of mostly the same patients, not independent biological replication cohorts. The result does not support a claim of widespread recurrent single-gene discordance in COAD under the tested definition.

### 3. HPA-CPTAC detection disagreement is widespread and recurrent

The cross-assay atlas contains 37,101 gene-cancer pairs, 11,267 unique genes, and 419 CPTAC tumors across four cancers.

| Cancer | Tumors | HPA-CPTAC overlap | Primary disagreement rate | Among HPA-low genes |
|---|---:|---:|---:|---:|
| COAD | 96 | 7,768 | 15.1% | 58.5% |
| BRCA | 119 | 9,923 | 20.3% | 64.3% |
| GBM | 99 | 9,924 | 33.6% | 74.9% |
| PAAD | 105 | 9,486 | 20.0% | 61.6% |

The primary event is defined as:

```text
HPA detection fraction <= 0.49
AND
CPTAC protein coverage >= 0.80
```

At the primary threshold, 2,135 genes recur in at least two cancers, 1,303 in at least three cancers, and 604 in all four cancers.

At the strictest threshold - HPA detection fraction 0 and CPTAC coverage 1 - 40.5%-59.1% of HPA-zero genes are nevertheless quantified in every CPTAC tumor. Ninety-eight genes meet this strict definition in all four cancers.

### 4. Cross-assay disagreement is compartment-dependent

Localization effects were validated using gene-clustered generalized estimating equations (GEE), cancer adjustment, gene length, HPA patient count, relevant coverage covariates, FDR correction, and restriction to genes with a localization annotation.

Among CPTAC-quantified proteins, HPA negativity was associated with:

| Localization | Odds ratio | 95% CI | Direction |
|---|---:|---:|---|
| Predicted secreted | 4.02 | 3.33-4.85 | More HPA-negative |
| Cell junctions | 1.66 | 1.31-2.11 | More HPA-negative |
| Actin filaments | 1.65 | 1.24-2.21 | More HPA-negative |
| Plasma membrane | 1.64 | 1.48-1.83 | More HPA-negative |
| Vesicles | 1.34 | 1.21-1.49 | More HPA-negative |
| Nucleoplasm | 0.72 | 0.66-0.79 | Less HPA-negative |
| Mitochondria | 0.64 | 0.55-0.74 | Less HPA-negative |
| Nuclear speckles | 0.48 | 0.38-0.63 | Less HPA-negative |

These are adjusted marginal associations from separate localization models; they are not independent causal effects.

## Data sources

- [UCSC Xena Toil recompute](https://xenabrowser.net/datapages/): TCGA, TARGET, and GTEx RNA expression and phenotype data
- [Human Protein Atlas](https://v23.proteinatlas.org/about/download): pinned HPA v23 tumor pathology IHC and subcellular localization
- [GENCODE](https://www.gencodegenes.org/human/release_23.html): human release 23 gene annotation
- [CPTAC](https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer): matched cancer transcriptomics and mass-spectrometry proteomics

Downloaded source URLs are declared in `config/downloads.tsv`. File sizes, timestamps, and SHA-256 checksums are recorded in `data/raw/download_inventory.tsv`.

## Reproduce the HPA/Xena analysis

Python 3.9 or newer and `curl` are required. A virtual environment is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
./run_all.sh
```

The primary Xena expression archive is approximately 1.32 GB. The pipeline is restartable:

```bash
# Show commands without running them
./run_all.sh --dry-run

# Quick bootstrap test
RNA_PROTEIN_SKIP_INSTALL=1 ./run_all.sh --bootstrap 20

# Continue completed stages after interruption
RNA_PROTEIN_SKIP_INSTALL=1 ./run_all.sh --resume

# Run only the final models
RNA_PROTEIN_SKIP_INSTALL=1 ./run_all.sh \
  --from-stage loco_logistic \
  --bootstrap 2000
```

Per-stage logs are written to `outputs/logs/`; run manifests and successful stage markers are written to `outputs/logs/` and `outputs/state/`.

## Reproduce the CPTAC extension

Download the four overlapping cohorts:

```bash
python3 scripts/cptac_pilot.py --cancer coad --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer brca --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer gbm  --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer pdac --rna-source washu --protein-source bcm
```

Run the COAD patient-level analyses:

```bash
python3 -m src.analyze_cptac_discordance --cancer coad
python3 -m src.analyze_cptac_robustness
```

Build and validate the four-cancer atlas:

```bash
python3 -m src.build_cross_cancer_disagreement_atlas
python3 -m src.validate_cross_assay_disagreement
```

## Important outputs

| Output | Purpose |
|---|---|
| `outputs/reports/dataset_audit.md` | Dataset dimensions, duplication, missingness, and coverage audit |
| `outputs/tables/loco_metrics_logistic_all.csv` | Seven-cancer LOCO performance |
| `outputs/tables/loco_bootstrap_summary_logistic_all.csv` | Paired bootstrap effect estimates |
| `outputs/tables/cptac_coad_robustness_classes.csv` | COAD processing/composition sensitivity |
| `outputs/tables/cross_cancer_hpa_cptac_atlas.csv` | Gene-cancer cross-assay atlas |
| `outputs/tables/cross_cancer_disagreement_recurrence.csv` | Cross-cancer recurrence |
| `outputs/tables/cross_assay_threshold_sensitivity.csv` | Sixteen threshold combinations |
| `outputs/tables/cross_assay_localization_gee.csv` | Gene-clustered localization effects |
| `outputs/reports/cross_assay_validation_overview.json` | Validated primary summary |

Large downloaded matrices and patient-level prediction tables are excluded from Git. Small result summaries and manifests can be versioned.

## Repository structure

```text
rna_protein_predictor/
├── config/                 # Download declarations
├── data/
│   ├── raw/                # Downloaded inputs
│   ├── processed/          # HPA/Xena analysis tables
│   └── cptac/              # Downloaded CPTAC cohort matrices
├── docs/                   # Analysis plan, decisions, reproducibility
├── notebooks/              # Early exploratory scripts
├── outputs/
│   ├── logs/               # Per-stage execution logs
│   ├── reports/            # JSON and Markdown summaries
│   └── tables/             # Audited analytical results
├── scripts/                # Download and orchestration entry points
├── src/                    # Reusable analysis modules
├── Makefile
├── requirements.txt
└── run_all.sh
```

## Interpretation boundaries

- HPA IHC non-detection is not equivalent to biological protein absence.
- CPTAC coverage indicates successful MS quantification, not absolute abundance.
- HPA and CPTAC use different patients, tissue sections, and assay technologies.
- Bulk-tissue composition can influence RNA and protein measurements.
- HPA-derived localization may share ascertainment processes with HPA pathology; independent GO/UniProt validation is still required.
- The current GEE models test localization features separately.
- Empirical bootstrap P-values are bounded by the replicate count; a zero count is reported as P < 1 / number of replicates, not P = 0.

## Remaining publication milestones

1. Freeze and characterize the 98 strict recurrent genes.
2. Repeat localization analyses using independent GO/UniProt annotations.
3. Validate protein quantification using an independent UMich/PDC workflow.
4. Test HPA release-version sensitivity.
5. Analyze the available PDAC tumor/normal samples separately.
6. Generate the primary multi-panel figure and supplementary robustness figures.
7. Expand biological interpretation to pathways, druggability, and protein class only after independent validation.

## Citation and data-use statement

Data used in this project were generated by TCGA, TARGET, GTEx, the Human Protein Atlas, and the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). Publications should cite the specific source datasets and CPTAC cohort papers used in the final frozen analysis.

## Historical proof of concept

Earlier repository versions reported a single train/test split, an approximately 0.10 ROC-AUC localization improvement, and threshold-derived discordant-gene lists. Those numbers are superseded by the audited LOCO, bootstrap, matched proteomics, FDR, and gene-clustered analyses documented above.
