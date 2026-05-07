# eda-code

Exploratory data analysis on the processed phyloseq objects.

`eda.qmd` reads `data/processed-data/ps_filt.rds` and `ps_rel.rds`, then
produces:

- `results/figures/seq_depth.png` — sequencing-depth histogram (motivates
  rarefaction in `analysis-code/`).
- `results/figures/top_phyla.png` — top-10 phyla relative abundance.
- `results/figures/alpha_diversity_5sites.png` and `shannon_boxplot.png` —
  per-body-site alpha diversity (preview of H1).
- `results/figures/pcoa_preview.png` — PCoA preview (Bray-Curtis, before
  formal PERMANOVA).
- `results/tables/exploratory_summary_by_subsite.rds` — sample size, mean
  library size, and mean Observed/Shannon/Simpson by body site.

The formal hypothesis tests are run from `code/analysis-code/`, not here.
