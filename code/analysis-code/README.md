# analysis-code

Stage 3 of the pipeline. `statistical-analysis.R` is one self-contained R
script that runs every formal analysis reported in the manuscript:

- **H1 — Alpha diversity (Observed, Shannon, Simpson, Chao1):**
  Kruskal-Wallis tests with Benjamini-Hochberg-adjusted pairwise Wilcoxon
  follow-ups.
- **H2 — Beta diversity:** Bray-Curtis and Jaccard distances, PERMANOVA
  (`vegan::adonis2`) with BETADISPER, PCoA (`cmdscale`) and Ward-D2
  hierarchical clustering.
- **H3 — Co-occurrence networks:** per-body-site Spearman networks
  (`Hmisc::rcorr`), thresholded at |r| ≥ 0.5 / p < 0.01, with Louvain
  community detection (`igraph::cluster_louvain`) and standard topology
  summaries.
- **H4 — Supervised classification (tidymodels):** four candidate models
  (null baseline, elastic-net multinomial GLM, single decision tree,
  random forest, XGBoost), tuned by 5-fold cross-validation on a
  stratified 75/25 train/test split. The script saves CV metrics for all
  models, the test-set confusion matrix and ROC curves for the winning
  model, permutation-based variable importance, and bootstrap 95%
  confidence intervals on the held-out accuracy and macro-AUC.

All randomness uses `set.seed(123)`. Outputs are written to
`results/figures/`, `results/tables/`, and `results/output/`.
Run after `processingfile-v1.qmd` has produced `ps_filt.rds`.
