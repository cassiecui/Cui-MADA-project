# processed-data

Created by `code/processing-code/processingfile-v1.qmd` from the live HMP
download. Both files are `phyloseq` objects.

| File | Description |
|------|-------------|
| `ps_filt.rds` | Filtered count-level object: zero-count taxa removed, taxa present in ≤ 5 samples removed. |
| `ps_rel.rds`  | Same samples and taxa, transformed to within-sample relative abundance (`x / sum(x)`). |

The simplified `SampleType` factor (Airways / Gut / Oral / Skin / Urogenital)
is added downstream in `eda.qmd` and `statistical-analysis.R`. Re-running the
processing Quarto file recreates these RDS files identically thanks to fixed
seeds.
