# code

All analysis code for the project. The pipeline is split into three stages,
each in its own subfolder. Run them in order.

```
code/
├── processing-code/   # Stage 1: download HMP data, build phyloseq objects
├── eda-code/          # Stage 2: exploratory data analysis
└── analysis-code/     # Stage 3: formal hypothesis tests + ML modeling
```

## Run order

1. `processing-code/processingfile-v1.qmd` — downloads `HMP16SData::V35()`,
   filters rare taxa, and writes
   `data/processed-data/ps_filt.rds` and `ps_rel.rds`.
2. `eda-code/eda.qmd` — sequencing-depth, top-phyla, alpha-diversity boxplots,
   and the `exploratory_summary_by_subsite.rds` table.
3. `analysis-code/statistical-analysis.R` — Kruskal-Wallis + pairwise Wilcoxon
   for alpha diversity (H1), PERMANOVA + BETADISPER for beta diversity (H2),
   per-site co-occurrence networks (H3), and the supervised body-site
   classifier with cross-validation, model comparison, and bootstrap-based
   uncertainty (H4).

All randomness uses `set.seed(123)`. Outputs are written to
`results/figures/`, `results/tables/`, and `results/output/`.
