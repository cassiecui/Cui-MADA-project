# data

```
data/
├── raw-data/          # see raw-data/README.md
└── processed-data/    # see processed-data/readme.md
```

The primary input for this project is the Human Microbiome Project 16S V3–V5
dataset, which is downloaded programmatically at runtime via
`HMP16SData::V35()` inside `code/processing-code/processingfile-v1.qmd`. As a
result there is no large raw FASTQ/biom file checked into the repository — a
small placeholder is kept under `raw-data/` for documentation only.

`processed-data/` holds the two `phyloseq` snapshots produced by the
processing step:

- `ps_filt.rds` — count-level data after dropping zero-count and rare taxa.
- `ps_rel.rds` — relative-abundance copy of the same object.

Both files are recreated every time the processing Quarto file is rendered.
