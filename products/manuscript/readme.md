# manuscript

`Manuscript.qmd` is the main paper:
*Comparing Microbial Diversity and Co-Occurrence Networks Across Five Major
Body Sites in the Human Microbiome Project.* It renders to both `html` and
`docx` and pulls every figure and table from `results/` via `here::here()`.

`supplement/Supplementary-Material.qmd` is the methods/results supplement
(extra method detail, pairwise alpha-diversity contrasts, full network
statistics, ML model-comparison tables and ROC curves with bootstrap CIs).
It renders to `html` and `pdf`.

Both files use:

- bibliography: `../../assets/dataanalysis-references.bib`
- citation style: `../../assets/american-journal-of-epidemiology.csl`
  (main) and `../../../assets/vancouver-author-date.csl` (supplement).
