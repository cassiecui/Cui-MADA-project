# products

Public-facing outputs of the project.

```
products/
├── manuscript/        # Main paper (Manuscript.qmd → docx + html) and supplement
├── presentation/      # revealjs slides (presentation.qmd) for the class talk
├── poster/            # placeholder for a conference poster
└── report/            # not used in this project (kept empty)
```

The manuscript and the supplement read figures and tables directly from
`results/figures/` and `results/tables/` via `here::here()`, so re-running
the analysis pipeline automatically updates every figure and number in the
paper without any manual copy/paste.

Both `assets/dataanalysis-references.bib` and the CSL files in `assets/`
are referenced from each Quarto file in this folder.
