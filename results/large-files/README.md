# large-files

`.gitignore`d folder for output files larger than ~20 MB (e.g. tuning
results, bootstrap resample objects). The folder itself and this README are
tracked; nothing else is.

In this project the largest committed object is the per-site network list
in `results/output/networks.rds` (~60 KB), so no file is currently being
parked here. The folder is kept so that re-running the pipeline with a
larger `N_PER_SITE` or finer hyperparameter grid has somewhere to write
intermediate objects without polluting git.
