# tables

RDS files containing every table reported in the manuscript and supplement.

| File | Source | Description |
|------|--------|-------------|
| `exploratory_summary_by_subsite.rds` | eda.qmd | Per-site sample size, library size, mean alpha diversity. |
| `alpha_kruskal.rds` | statistical-analysis.R (H1) | Kruskal-Wallis chi-sq, df, p-values for the four alpha-diversity indices. |
| `alpha_pairwise.rds` | statistical-analysis.R (H1) | BH-adjusted pairwise Wilcoxon p-values. |
| `permanova.rds` | statistical-analysis.R (H2) | PERMANOVA + BETADISPER for Bray-Curtis and Jaccard. |
| `network_stats.rds` | statistical-analysis.R (H3) | Per-site network topology (nodes, edges, density, modularity, ...). |
| `ml_cv_metrics.rds` | statistical-analysis.R (H4) | Mean ± SD CV accuracy and macro-AUC for every candidate model. |
| `ml_test_metrics.rds` | statistical-analysis.R (H4) | Held-out test-set accuracy and AUC for every candidate model. |
| `ml_test_uncertainty.rds` | statistical-analysis.R (H4) | Bootstrap 95% CIs for the winning model's test accuracy and AUC. |
| `ml_best_params.rds` | statistical-analysis.R (H4) | Tuned hyperparameters for each model. |
| `ml_confusion_matrix.rds` | statistical-analysis.R (H4) | `yardstick::conf_mat` object for the winning model. |
| `ml_variable_importance.rds` | statistical-analysis.R (H4) | Top-20 permutation-based variable importance scores. |
