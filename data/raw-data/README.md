# raw-data

The "raw" data for this project is the Human Microbiome Project 16S rRNA
V3–V5 dataset, fetched programmatically by
`code/processing-code/processingfile-v1.qmd` via:

```r
se <- HMP16SData::V35()       # SummarizedExperiment from Bioconductor
ps <- microbiome::as_phyloseq(se)
```

Because the download is fully reproducible (the package version pins the
release), no FASTQ/biom file is committed here. The only file in this folder
is `exampledata.xlsx`, which is a tiny example dataset left over from the
project template; it is **not used** in any of the analyses and is kept only
because the supplement references the template's structure.

The codebook for the HMP variables we actually use (sample metadata column
`HMP_BODY_SITE`, simplified to Airways / Gut / Oral / Skin / Urogenital) is
documented in the manuscript and in `processingfile-v1.qmd`.
