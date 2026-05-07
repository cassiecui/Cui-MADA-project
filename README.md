# Cui MADA Project — HMP 16S Microbiome Analysis Across Five Major Body Sites

This repository contains the full data-analysis project for the MADA course
(Modern Applied Data Analysis,
[course content overview](https://andreashandel.github.io/MADAcourse/content/content-overview.html)).
The project uses 16S rRNA V3–V5 data from the Human Microbiome Project (HMP)
and tests whether alpha diversity, community composition, microbial
co-occurrence network structure, and supervised classifiability differ
across five major body sites: Airways, Gut, Oral, Skin and Urogenital.

The pipeline is designed to demonstrate, end-to-end, the major MADA topics:
data wrangling and cleaning, exploratory data analysis, statistical
hypothesis testing, machine learning model fitting, regularization,
ensemble methods, cross-validation, model comparison and selection,
performance evaluation (confusion matrix + ROC), permutation variable
importance, and bootstrap-based uncertainty quantification.

Main research questions:

- **H1** — Alpha diversity differs across the five body sites
  (Kruskal-Wallis + BH-pairwise Wilcoxon).
- **H2** — Samples cluster by body site in beta-diversity ordination space
  (PERMANOVA + BETADISPER on Bray-Curtis and Jaccard, PCoA, hierarchical
  clustering).
- **H3** — Co-occurrence networks differ in structure across body sites
  (per-site Spearman networks, Louvain community detection, topology
  metrics).
- **H4** — Body site can be predicted from taxon abundances; tree-based
  ensembles outperform a regularized GLM and a single decision tree
  (`tidymodels` workflow with 5-fold CV and bootstrap CIs).

## Repository structure

```
Cui-MADA-project/
├── assets/                       # Bibliography, CSL, schematics
├── code/
│   ├── processing-code/          # HMP data download + cleaning (step 1)
│   ├── eda-code/                 # Exploratory data analysis (step 2)
│   └── analysis-code/            # Formal stats + networks (step 3)
├── data/
│   ├── raw-data/                 # Inputs (HMP is downloaded at runtime)
│   └── processed-data/           # ps_filt.rds, ps_rel.rds
├── products/
│   ├── manuscript/               # Manuscript.qmd + supplement
│   ├── report/                   # HTML report
│   ├── presentation/             # Slides
│   └── poster/                   # Poster
└── results/
    ├── figures/                  # All figures used in the manuscript
    ├── tables/                   # All tables used in the manuscript (RDS)
    └── output/                   # Larger intermediate objects (networks)
```

## Prerequisites

- R ≥ 4.3 and Quarto
- The following R packages (CRAN + Bioconductor):

  ```r
  install.packages(c(
    "tidyverse", "here", "broom", "vegan", "igraph", "ggraph",
    "tidygraph", "Hmisc", "RColorBrewer", "skimr", "readxl",
    "knitr", "rmarkdown",
    # H4: tidymodels workflow + extra learners
    "tidymodels", "glmnet", "ranger", "rpart", "vip",
    "xgboost"     # optional; H4 skips XGBoost if not installed
  ))

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c(
    "phyloseq", "HMP16SData", "SummarizedExperiment", "microbiome"
  ))
  ```

## Reproduction instructions

Run the following three steps in order. Each step writes its outputs to
`data/processed-data/` or `results/`. The working directory is set
automatically by the `here` package so commands can be run from anywhere in
the project.

1. **Processing** — downloads the HMP V3–V5 data via `HMP16SData::V35()`,
   converts it to a `phyloseq` object, removes zero-count taxa, filters
   rare taxa (present in ≤ 5 samples), transforms to relative abundance,
   and writes `data/processed-data/ps_filt.rds` and
   `data/processed-data/ps_rel.rds`:

   ```r
   quarto::quarto_render(
     here::here("code", "processing-code", "processingfile-v1.qmd"))
   ```

2. **Exploratory data analysis** — renders sequencing-depth, top-phyla,
   alpha-diversity and PCoA-preview figures, and writes the
   `exploratory_summary_by_subsite.rds` table:

   ```r
   quarto::quarto_render(here::here("code", "eda-code", "eda.qmd"))
   ```

3. **Formal statistical analysis (H1–H4)** — runs Kruskal-Wallis + pairwise
   Wilcoxon for alpha diversity (H1), PERMANOVA + BETADISPER for beta
   diversity (H2), per-site co-occurrence networks (H3), and the
   `tidymodels` body-site classifier (H4) with five candidate models
   (null, elastic-net GLM, decision tree, random forest, XGBoost), 5-fold
   cross-validation, model comparison, ROC curves, permutation variable
   importance, and bootstrap 95% CIs on the held-out test performance.
   Writes every table and figure used in the manuscript:

   ```r
   source(here::here("code", "analysis-code", "statistical-analysis.R"))
   ```

4. **Manuscript and supplement** — rendering these files pulls in the
   saved tables and figures:

   ```r
   quarto::quarto_render(
     here::here("products", "manuscript", "Manuscript.qmd"))
   quarto::quarto_render(
     here::here("products", "manuscript", "supplement",
                "Supplementary-Material.qmd"))
   ```

## File conventions

- Lower-case names with `-` as the word separator.
- R scripts end in `.R`; Quarto files in `.qmd`.
- Tables are saved as `.rds`; figures as `.png`.
- All randomness uses `set.seed(123)` for reproducibility.

## Data

The primary input is the HMP 16S V3–V5 dataset obtained at runtime via
`HMP16SData::V35()`. No raw files need to be stored in the repository; the
download is reproducible and checked against the package version.

## MADA topics covered

| Topic (course content overview) | Where it appears in the project |
|---------------------------------|---------------------------------|
| Data wrangling / cleaning | `code/processing-code/processingfile-v1.qmd` |
| Exploratory data analysis | `code/eda-code/eda.qmd` |
| Visualization and tables | All figures in `results/figures/` and tables in `results/tables/` |
| Statistical hypothesis testing | H1 (Kruskal-Wallis + Wilcoxon), H2 (PERMANOVA + BETADISPER), H3 (network topology) in `statistical-analysis.R` |
| Model fitting and tuning | H4 — `tune_grid()` over five candidate workflows |
| Cross-validation | H4 — 5-fold stratified CV (`vfold_cv`) on the training split |
| Regularization | H4 — elastic-net `multinom_reg(penalty, mixture)` |
| Tree-based / ensemble methods | H4 — single decision tree, random forest, XGBoost |
| Model comparison and selection | H4 — `ml_cv_metrics.rds`, `ml_cv_comparison.png` |
| Performance metrics | Accuracy, macro-AUC, confusion matrix, one-vs-rest ROC |
| Variable importance | H4 — `vip::vi`, `ml_variable_importance.png` |
| Uncertainty quantification | H4 — bootstrap 95% CIs (B = 1000) on test accuracy and macro-AUC; PERMANOVA permutations and BETADISPER for H2 |
| Reproducibility | `here::here()`, fixed `set.seed(123)`, programmatic data download |

## Citations

Example reproducible projects referenced here include McKay et al.
(2020a, 2020b); see `assets/dataanalysis-references.bib`.
