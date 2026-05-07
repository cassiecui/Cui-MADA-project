# processing-code

Stage 1 of the pipeline: pull the Human Microbiome Project (HMP) 16S V3–V5
dataset, build a `phyloseq` object, and apply basic filtering.

`processingfile-v1.qmd` is the working file. It:

1. Loads `HMP16SData::V35()` and converts the `SummarizedExperiment` to a
   `phyloseq` object via `microbiome::as_phyloseq()`.
2. Removes zero-count taxa and rare taxa (present in ≤ 5 samples).
3. Builds a relative-abundance copy.
4. Writes `data/processed-data/ps_filt.rds` and
   `data/processed-data/ps_rel.rds`.

The HMP data are downloaded at runtime — no large raw file is committed to
the repo. Subsampling, rarefaction and the body-site simplification used in
the manuscript are applied later, inside
`code/analysis-code/statistical-analysis.R`, so `ps_filt.rds` reflects only
the count-level filtering described above.
