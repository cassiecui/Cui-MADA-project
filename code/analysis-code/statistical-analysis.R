###############################################################################
# statistical-analysis.R
#
# Formal statistical analysis for the 5-body-site HMP microbiome project.
# This script loads the processed phyloseq data and performs:
#   - H1: Alpha diversity comparison across body sites (Kruskal-Wallis,
#         pairwise Wilcoxon with BH correction)
#   - H2: Beta diversity / clustering (Bray-Curtis + Jaccard; PERMANOVA with
#         BETADISPER; PCoA; hierarchical clustering)
#   - H3: Microbial co-occurrence networks per body site (Spearman-based
#         correlation networks with Louvain community detection)
#
# All results (tables, figures, model objects) are written to
# results/tables and results/figures.
#
# Usage:  Run after processingcode.R / processingfile-v2.qmd have produced
#         data/processed-data/ps_filt.rds
###############################################################################

## ---- packages --------
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(phyloseq)
  library(vegan)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(Hmisc)
  library(broom)
  library(RColorBrewer)
  library(tidymodels)
})
tidymodels_prefer()

options(lifecycle_verbosity = "quiet")
theme_set(theme_bw(base_size = 12))

colors_body <- c(
  "Airways"    = "#87CEEB",
  "Gut"        = "#8B4513",
  "Oral"       = "#DC143C",
  "Skin"       = "#FFB6C1",
  "Urogenital" = "#9370DB"
)

fig_dir   <- here("results", "figures")
tab_dir   <- here("results", "tables")
out_dir   <- here("results", "output")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## ---- loaddata --------
# IMPORTANT: ps_filt.rds from the processing step may contain the full
# filtered HMP dataset (~4500 samples x ~34000 taxa). Running PERMANOVA /
# network analysis on that is impractical (the Bray distance alone has
# ~10M pairs). We therefore apply the stratified subsampling + prevalence
# filter + rarefaction here, mirroring the group's original R script. The
# size parameters are tuned so the full pipeline completes in a few minutes.

# --- Tunable knobs (keep small so the pipeline finishes reliably) -----------
N_PER_SITE       <- 200    # max samples per body site
MIN_TAXA_PREV    <- 10     # keep taxa present in >= this many samples
MIN_TAXA_ABUND   <- 50     # keep taxa with total count >= this
RAREFY_PCTILE    <- 0.05   # rarefy to this percentile of sample depths
N_PERM           <- 99     # permutations for PERMANOVA / BETADISPER
# ---------------------------------------------------------------------------

ps <- readRDS(here("data", "processed-data", "ps_filt.rds"))

# Ensure the simplified body-site column exists
if (!"SampleType" %in% colnames(sample_data(ps))) {
  sample_data(ps)$SampleType <- sample_data(ps)$HMP_BODY_SITE |>
    gsub("Gastrointestinal Tract", "Gut", x = _) |>
    gsub("Urogenital Tract", "Urogenital", x = _)
}

# Keep only the five major body sites
keep_sites <- names(colors_body)
ps <- prune_samples(
  sample_data(ps)$SampleType %in% keep_sites, ps
)
cat(sprintf("[setup] %d samples across %d sites before subsampling\n",
            nsamples(ps), length(unique(sample_data(ps)$SampleType))))

# Stratified subsample: up to N_PER_SITE per site (or all, if fewer)
# NOTE: avoid phyloseq::subset_samples() here — it uses NSE and won't see
# `site` from a closure reliably.
set.seed(123)
site_vec <- as.character(sample_data(ps)$SampleType)
names(site_vec) <- sample_names(ps)
selected <- unlist(lapply(keep_sites, function(s) {
  site_names <- names(site_vec)[site_vec == s]
  if (length(site_names) == 0) return(character(0))
  sample(site_names, min(N_PER_SITE, length(site_names)))
}))
ps <- prune_samples(selected, ps)

# Drop rare taxa (present in too few samples OR too low total abundance)
cat(sprintf("[setup] %d samples, %d taxa before prevalence/abundance filter\n",
            nsamples(ps), ntaxa(ps)))
ps <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= MIN_TAXA_PREV & sum(x) >= MIN_TAXA_ABUND,
  prune = TRUE
)
ps <- prune_samples(sample_sums(ps) > 0, ps)
cat(sprintf("[setup] %d samples, %d taxa after filter\n",
            nsamples(ps), ntaxa(ps)))

# Rarefy to an even depth (RAREFY_PCTILE-th percentile) so that
# alpha/beta diversity are not confounded by sequencing depth.
set.seed(123)
min_depth <- sort(sample_sums(ps))[max(1, floor(nsamples(ps) * RAREFY_PCTILE))]
cat(sprintf("[setup] rarefying to %d reads/sample\n", min_depth))
ps_rare <- rarefy_even_depth(ps, sample.size = min_depth, rngseed = 123,
                             verbose = FALSE)
cat(sprintf("[setup] final dataset: %d samples x %d taxa\n",
            nsamples(ps_rare), ntaxa(ps_rare)))


###############################################################################
## H1: Alpha diversity differs by body site
###############################################################################

## ---- alpha_diversity --------
alpha_div <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Shannon", "Simpson", "Chao1")
)
alpha_div$SampleType <- sample_data(ps_rare)$SampleType

metrics <- c("Observed", "Shannon", "Simpson", "Chao1")

# Overall Kruskal-Wallis tests
kw_table <- purrr::map_dfr(metrics, function(m) {
  kw <- kruskal.test(as.formula(paste(m, "~ SampleType")), data = alpha_div)
  tibble(metric = m,
         chi_sq = unname(kw$statistic),
         df     = unname(kw$parameter),
         p      = kw$p.value)
})
saveRDS(kw_table, file.path(tab_dir, "alpha_kruskal.rds"))

# Pairwise Wilcoxon with Benjamini-Hochberg correction (long format)
pairwise_list <- lapply(metrics, function(m) {
  pw <- pairwise.wilcox.test(alpha_div[[m]], alpha_div$SampleType,
                             p.adjust.method = "BH")
  pw_df <- as.data.frame(as.table(pw$p.value))
  names(pw_df) <- c("site_1", "site_2", "p_adj")
  pw_df$metric <- m
  pw_df |> filter(!is.na(p_adj))
})
pairwise_tbl <- bind_rows(pairwise_list) |>
  select(metric, site_1, site_2, p_adj)
saveRDS(pairwise_tbl, file.path(tab_dir, "alpha_pairwise.rds"))

# Boxplot across all four indices
alpha_long <- alpha_div |>
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long,
                  aes(x = SampleType, y = Value, fill = SampleType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = colors_body) +
  labs(title = "Alpha diversity across five major body sites",
       subtitle = sprintf("Rarefied to %d reads/sample", min_depth),
       x = "Body Site", y = "Diversity value") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "alpha_diversity_5sites.png"),
       p_alpha, width = 10, height = 8, dpi = 300)

###############################################################################
## H2: Samples cluster by body site (beta diversity)
###############################################################################
# Note: the previous version used phyloseq::distance() / ordinate() repeatedly,
# which recomputed the Bray-Curtis distance matrix up to three times and
# sometimes hung on large objects. This version extracts the OTU matrix
# ONCE, uses vegan::vegdist directly, and reuses the same distance matrix for
# PERMANOVA, BETADISPER, PCoA (via cmdscale) and hclust.

## ---- beta_prep --------
cat("[H2] Preparing OTU matrices...\n")

otu_counts <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_counts <- t(otu_counts)   # rows = samples
storage.mode(otu_counts) <- "double"

otu_log <- log1p(otu_counts)                              # log(1+x)
otu_pa  <- (otu_counts > 0) + 0L                          # presence/absence

metadata <- data.frame(sample_data(ps_rare))
metadata <- metadata[rownames(otu_counts), , drop = FALSE]
metadata$SampleType <- factor(metadata$SampleType)

cat(sprintf("[H2] %d samples x %d taxa\n",
            nrow(otu_counts), ncol(otu_counts)))

## ---- beta_permanova_bray --------
cat("[H2] Computing Bray-Curtis distance...\n")
dmat_bray <- vegan::vegdist(otu_log, method = "bray")

cat("[H2] PERMANOVA (Bray)...\n")
set.seed(123)
perm_bray <- vegan::adonis2(dmat_bray ~ SampleType,
                            data = metadata, permutations = N_PERM)

cat("[H2] BETADISPER (Bray)...\n")
set.seed(123)
bdisp_bray <- vegan::betadisper(dmat_bray, metadata$SampleType)
bperm_bray <- vegan::permutest(bdisp_bray, permutations = N_PERM)

## ---- beta_permanova_jaccard --------
cat("[H2] Computing Jaccard distance...\n")
dmat_jaccard <- vegan::vegdist(otu_pa, method = "jaccard", binary = TRUE)

cat("[H2] PERMANOVA (Jaccard)...\n")
set.seed(123)
perm_jaccard <- vegan::adonis2(dmat_jaccard ~ SampleType,
                               data = metadata, permutations = N_PERM)

cat("[H2] BETADISPER (Jaccard)...\n")
set.seed(123)
bdisp_jaccard <- vegan::betadisper(dmat_jaccard, metadata$SampleType)
bperm_jaccard <- vegan::permutest(bdisp_jaccard, permutations = N_PERM)

permanova_tbl <- bind_rows(
  tibble(distance     = "bray",
         R2           = perm_bray$R2[1],
         F            = perm_bray$F[1],
         p_permanova  = perm_bray$`Pr(>F)`[1],
         F_betadisper = bperm_bray$tab$F[1],
         p_betadisper = bperm_bray$tab$`Pr(>F)`[1]),
  tibble(distance     = "jaccard",
         R2           = perm_jaccard$R2[1],
         F            = perm_jaccard$F[1],
         p_permanova  = perm_jaccard$`Pr(>F)`[1],
         F_betadisper = bperm_jaccard$tab$F[1],
         p_betadisper = bperm_jaccard$tab$`Pr(>F)`[1])
)
saveRDS(permanova_tbl, file.path(tab_dir, "permanova.rds"))

## ---- pcoa --------
# Use cmdscale on the already-computed Bray distance; do NOT call ordinate()
# (which would recompute the distance matrix internally).
cat("[H2] PCoA (cmdscale on Bray)...\n")
pcoa_bray <- cmdscale(dmat_bray, k = 2, eig = TRUE)

eig_vals <- pcoa_bray$eig
var1 <- round(100 * eig_vals[1] / sum(abs(eig_vals)), 1)
var2 <- round(100 * eig_vals[2] / sum(abs(eig_vals)), 1)

ord_df <- data.frame(
  Axis.1     = pcoa_bray$points[, 1],
  Axis.2     = pcoa_bray$points[, 2],
  SampleType = metadata$SampleType
)

p_pcoa <- ggplot(ord_df,
                 aes(x = Axis.1, y = Axis.2, color = SampleType)) +
  geom_point(size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = colors_body) +
  labs(title = "PCoA (Bray-Curtis) - Five body sites",
       subtitle = sprintf(
         "PERMANOVA R2 = %.3f, p = %.3g",
         permanova_tbl$R2[permanova_tbl$distance == "bray"],
         permanova_tbl$p_permanova[permanova_tbl$distance == "bray"]
       ),
       x = sprintf("PC1 (%.1f%% variance)", var1),
       y = sprintf("PC2 (%.1f%% variance)", var2),
       color = "Body Site")

ggsave(file.path(fig_dir, "pcoa_ordination_5sites.png"),
       plot = p_pcoa, width = 9, height = 7, dpi = 300)

## ---- hclust --------
cat("[H2] Hierarchical clustering (Ward D2 on Bray)...\n")
hc <- hclust(dmat_bray, method = "ward.D2")

sample_types  <- as.character(metadata$SampleType[hc$order])
sample_colors <- colors_body[sample_types]

png(file.path(fig_dir, "hierarchical_clustering_5sites.png"),
    width = 14, height = 10, units = "in", res = 300)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(8, 0.5, 0.5))

par(mar = c(0, 4, 4, 2))
plot(hc, labels = FALSE, xlab = "",
     main = sprintf("Hierarchical clustering - Five body sites (n = %d)",
                    nrow(otu_counts)),
     ylab = "Bray-Curtis distance (Ward's method)")

par(mar = c(0, 4, 0, 2))
barplot(rep(1, length(sample_colors)), col = sample_colors, border = NA,
        space = 0, axes = FALSE, ylab = "")

par(mar = c(1, 0, 0, 0))
plot.new()
legend("center", legend = names(colors_body), fill = colors_body,
       cex = 1.3, title = "Body Site", bty = "n", horiz = TRUE)
dev.off()

cat("[H2] Done.\n")


###############################################################################
## H3: Microbial co-occurrence networks
###############################################################################

## ---- network_function --------
build_cooccurrence_network <- function(ps_subset, site_name,
                                       cor_threshold = 0.5,
                                       p_threshold   = 0.01,
                                       max_taxa      = 500) {

  otu_mat <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) otu_mat <- t(otu_mat)

  # Prevalence filter: taxa must be present in at least 10% of samples
  prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
  otu_mat    <- otu_mat[, prevalence >= 0.1, drop = FALSE]

  # Cap the number of taxa by total abundance so rcorr() does not blow up
  # on very dense sites. This keeps runtime bounded regardless of input size.
  if (ncol(otu_mat) > max_taxa) {
    totals <- colSums(otu_mat)
    keep   <- names(sort(totals, decreasing = TRUE))[seq_len(max_taxa)]
    otu_mat <- otu_mat[, keep, drop = FALSE]
  }
  cat(sprintf("[H3] %s: %d samples, %d taxa (after cap)\n",
              site_name, nrow(otu_mat), ncol(otu_mat)))

  if (ncol(otu_mat) < 3) return(NULL)

  cor_result <- Hmisc::rcorr(otu_mat, type = "spearman")
  cor_matrix <- cor_result$r
  p_matrix   <- cor_result$P

  cor_matrix[abs(cor_matrix) < cor_threshold] <- 0
  cor_matrix[p_matrix >= p_threshold]         <- 0
  diag(cor_matrix) <- 0

  g <- igraph::graph_from_adjacency_matrix(
    cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  g <- igraph::delete.vertices(g, which(igraph::degree(g) == 0))

  if (igraph::vcount(g) == 0) return(g)

  # Attach taxonomy
  tax_info     <- as.data.frame(tax_table(ps_subset))
  vertex_names <- igraph::V(g)$name
  for (col in c("Phylum", "Family", "Genus")) {
    if (col %in% colnames(tax_info)) {
      igraph::V(g)$.data <- tax_info[vertex_names, col]
      names(igraph::vertex_attr(g))[length(igraph::vertex_attr_names(g))] <- col
    }
  }

  # Centrality
  igraph::V(g)$degree      <- igraph::degree(g)
  igraph::V(g)$betweenness <- igraph::betweenness(g)
  igraph::V(g)$closeness   <- igraph::closeness(g)
  igraph::V(g)$eigenvector <- igraph::eigen_centrality(g)$vector

  # Edges
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$edge_type   <- ifelse(igraph::E(g)$weight > 0,
                                     "Positive", "Negative")

  # Communities (Louvain)
  communities <- igraph::cluster_louvain(g)
  igraph::V(g)$community <- igraph::membership(communities)
  attr(g, "site") <- site_name

  g
}

## ---- build_networks --------
networks <- list()
network_stats <- data.frame()

site_vec_rare <- as.character(sample_data(ps_rare)$SampleType)
names(site_vec_rare) <- sample_names(ps_rare)

for (site in names(colors_body)) {
  site_samples <- names(site_vec_rare)[site_vec_rare == site]
  if (length(site_samples) == 0) next
  ps_site <- prune_samples(site_samples, ps_rare)
  ps_site <- prune_taxa(taxa_sums(ps_site) > 0, ps_site)

  if (nsamples(ps_site) < 20) next

  g <- build_cooccurrence_network(ps_site, site)
  if (is.null(g) || igraph::vcount(g) == 0) next

  networks[[site]] <- g
  network_stats <- rbind(network_stats, data.frame(
    BodySite    = site,
    Nodes       = igraph::vcount(g),
    Edges       = igraph::ecount(g),
    Density     = igraph::edge_density(g),
    AvgDegree   = mean(igraph::degree(g)),
    Clustering  = igraph::transitivity(g, type = "global"),
    Modularity  = igraph::modularity(igraph::cluster_louvain(g)),
    Communities = length(unique(igraph::V(g)$community))
  ))
}

saveRDS(network_stats, file.path(tab_dir, "network_stats.rds"))
saveRDS(networks,      file.path(out_dir, "networks.rds"))

## ---- network_plot --------
stats_long <- network_stats |>
  pivot_longer(cols = c(Nodes, Edges, Density, AvgDegree, Clustering,
                        Modularity, Communities),
               names_to = "Metric", values_to = "Value") |>
  mutate(Metric = factor(Metric,
                         levels = c("Nodes", "Edges", "Density", "AvgDegree",
                                    "Clustering", "Modularity", "Communities")))

p_network_stats <- ggplot(stats_long,
                          aes(x = BodySite, y = Value, fill = BodySite)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = colors_body) +
  labs(title = "Co-occurrence network properties across body sites",
       subtitle = "Spearman |r| > 0.5, p < 0.01",
       x = "Body Site", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "network_statistics_comparison.png"),
       p_network_stats, width = 12, height = 8, dpi = 300)


###############################################################################
## H4: Supervised classification of body site with tidymodels
###############################################################################
# Goal: showcase the MADA topics on model fitting, cross-validation,
# model comparison, performance metrics, regularization, ensembles, variable
# importance and uncertainty.
#
# Pipeline:
#   - features: top-300 most-abundant taxa, log(1+x), zero-variance and z-score
#   - 75/25 stratified train/test split
#   - 5-fold CV on the training data (stratified by body site)
#   - five candidate workflows compared on CV macro-AUC and accuracy:
#         (a) null baseline (predicts modal class)
#         (b) elastic-net multinomial GLM           [regularization]
#         (c) decision tree (CART)                  [interpretable single tree]
#         (d) random forest (ranger)                [bagged ensemble]
#         (e) XGBoost                               [boosted ensemble; optional]
#   - tuned hyperparameters are saved per model
#   - the winning model (highest CV macro-AUC) is refit on the full training
#     set and evaluated on the held-out test set; we additionally compute
#     bootstrap 95% CIs on test accuracy and macro-AUC, the per-class ROC
#     curves, the confusion matrix, and a permutation-based variable
#     importance ranking.

cat("[H4] Preparing feature matrix...\n")
otu_feat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_feat <- t(otu_feat)
top <- names(sort(colSums(otu_feat), decreasing = TRUE))[seq_len(min(300, ncol(otu_feat)))]
ml_df <- as_tibble(log1p(otu_feat[, top, drop = FALSE]),
                   .name_repair = ~ make.names(.x, unique = TRUE)) |>
  mutate(SampleType = factor(sample_data(ps_rare)$SampleType,
                             levels = names(colors_body)))

set.seed(123)
split      <- initial_split(ml_df, prop = 0.75, strata = SampleType)
train_data <- training(split)
folds      <- vfold_cv(train_data, v = 5, strata = SampleType)

rec <- recipe(SampleType ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

mset <- metric_set(accuracy, roc_auc)

## ---- model_specs --------
null_spec <- null_model() |>
  set_engine("parsnip") |>
  set_mode("classification")

glm_spec <- multinom_reg(penalty = tune(), mixture = tune()) |>
  set_engine("glmnet") |> set_mode("classification")

tree_spec <- decision_tree(cost_complexity = tune(),
                           tree_depth      = tune(),
                           min_n           = tune()) |>
  set_engine("rpart") |> set_mode("classification")

rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
  set_engine("ranger", importance = "permutation") |>
  set_mode("classification")

xgb_available <- requireNamespace("xgboost", quietly = TRUE)
if (xgb_available) {
  xgb_spec <- boost_tree(trees      = tune(),
                         tree_depth = tune(),
                         learn_rate = tune(),
                         min_n      = tune()) |>
    set_engine("xgboost") |> set_mode("classification")
} else {
  message("[H4] xgboost not installed; skipping the boosted-trees model.")
}

null_wf <- workflow() |> add_recipe(rec) |> add_model(null_spec)
glm_wf  <- workflow() |> add_recipe(rec) |> add_model(glm_spec)
tree_wf <- workflow() |> add_recipe(rec) |> add_model(tree_spec)
rf_wf   <- workflow() |> add_recipe(rec) |> add_model(rf_spec)
if (xgb_available) {
  xgb_wf <- workflow() |> add_recipe(rec) |> add_model(xgb_spec)
}

## ---- cv_tuning --------
ctrl <- control_grid(save_pred = TRUE, save_workflow = FALSE, verbose = FALSE)

cat("[H4] CV: null baseline...\n")
set.seed(123)
null_res <- fit_resamples(null_wf, folds, metrics = mset,
                          control = control_resamples(save_pred = TRUE))

cat("[H4] CV: elastic-net multinomial GLM (24 candidates)...\n")
set.seed(123)
glm_res <- tune_grid(
  glm_wf, folds, control = ctrl, metrics = mset,
  grid = grid_regular(penalty(c(-4, 0)), mixture(c(0, 1)),
                      levels = c(8, 3))
)

cat("[H4] CV: single decision tree (16 candidates)...\n")
set.seed(123)
tree_res <- tune_grid(
  tree_wf, folds, control = ctrl, metrics = mset,
  grid = grid_regular(cost_complexity(c(-4, -1)),
                      tree_depth(c(3L, 15L)),
                      min_n(c(2L, 20L)),
                      levels = c(4, 2, 2))
)

cat("[H4] CV: random forest (16 candidates)...\n")
set.seed(123)
rf_res <- tune_grid(
  rf_wf, folds, control = ctrl, metrics = mset,
  grid = grid_regular(mtry(c(5L, 100L)), min_n(c(2L, 20L)),
                      levels = c(4, 4))
)

if (xgb_available) {
  cat("[H4] CV: XGBoost (16 candidates)...\n")
  set.seed(123)
  xgb_res <- tune_grid(
    xgb_wf, folds, control = ctrl, metrics = mset,
    grid = grid_regular(trees(c(100L, 800L)),
                        tree_depth(c(3L, 9L)),
                        learn_rate(c(-2, -0.5)),
                        min_n(c(2L, 20L)),
                        levels = c(2, 2, 2, 2))
  )
}

## ---- cv_comparison --------
# For each model, take the hyperparameter setting with the highest CV
# macro-AUC and report mean +/- SD across the 5 folds for both metrics.
get_cv_summary <- function(res, model_name) {
  perfold <- collect_metrics(res, summarize = FALSE)
  best_id <- if (model_name == "null") {
    perfold |> distinct(.config) |> pull(.config) |> head(1)
  } else {
    select_best(res, metric = "roc_auc")$.config
  }
  perfold |>
    filter(.config == best_id) |>
    group_by(.metric) |>
    dplyr::summarize(mean = mean(.estimate, na.rm = TRUE),
                     sd   = sd(.estimate,   na.rm = TRUE),
                     n    = sum(!is.na(.estimate)),
                     .groups = "drop") |>
    mutate(model = model_name)
}

cv_summary <- bind_rows(
  get_cv_summary(null_res, "null"),
  get_cv_summary(glm_res,  "multinom_glmnet"),
  get_cv_summary(tree_res, "decision_tree"),
  get_cv_summary(rf_res,   "random_forest"),
  if (xgb_available) get_cv_summary(xgb_res, "xgboost")
)
saveRDS(cv_summary, file.path(tab_dir, "ml_cv_metrics.rds"))

cv_plot_df <- cv_summary |>
  mutate(model = factor(model,
                        levels = c("null", "multinom_glmnet", "decision_tree",
                                   "random_forest", "xgboost")))

p_cv <- ggplot(cv_plot_df,
               aes(x = model, y = mean,
                   ymin = mean - sd, ymax = mean + sd, color = model)) +
  geom_pointrange(size = 0.8, linewidth = 0.9) +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(title = "5-fold cross-validation: model comparison",
       subtitle = "Best hyperparameter setting per model (mean +/- 1 SD across folds)",
       x = NULL, y = "Metric value") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(fig_dir, "ml_cv_comparison.png"),
       p_cv, width = 9, height = 5, dpi = 300)

## ---- pick_winner --------
auc_means <- cv_summary |> filter(.metric == "roc_auc") |>
  arrange(desc(mean))
best_name <- auc_means$model[1]
cat(sprintf("[H4] Winning model on CV macro-AUC: %s (mean = %.3f)\n",
            best_name, auc_means$mean[1]))

# Save tuned hyperparameters for every (non-null) candidate
tuned_params <- bind_rows(
  tibble(model = "multinom_glmnet",
         param = names(select_best(glm_res, metric = "roc_auc")),
         value = as.character(unlist(select_best(glm_res, metric = "roc_auc")))),
  tibble(model = "decision_tree",
         param = names(select_best(tree_res, metric = "roc_auc")),
         value = as.character(unlist(select_best(tree_res, metric = "roc_auc")))),
  tibble(model = "random_forest",
         param = names(select_best(rf_res, metric = "roc_auc")),
         value = as.character(unlist(select_best(rf_res, metric = "roc_auc")))),
  if (xgb_available) tibble(model = "xgboost",
         param = names(select_best(xgb_res, metric = "roc_auc")),
         value = as.character(unlist(select_best(xgb_res, metric = "roc_auc"))))
) |> filter(param != ".config")
saveRDS(tuned_params, file.path(tab_dir, "ml_best_params.rds"))

## ---- test_metrics_all_models --------
# Refit each candidate on the full training set with its tuned hyperparams,
# then evaluate on the held-out test data.
finalize_and_test <- function(wf, par, model_name) {
  if (is.null(par)) {
    fit_obj <- wf |> last_fit(split, metrics = mset)
  } else {
    fit_obj <- finalize_workflow(wf, par) |> last_fit(split, metrics = mset)
  }
  list(fit = fit_obj,
       metrics = collect_metrics(fit_obj) |> mutate(model = model_name))
}

cat("[H4] Refit + test: all candidates...\n")
null_test <- finalize_and_test(null_wf, NULL,                                        "null")
glm_test  <- finalize_and_test(glm_wf,  select_best(glm_res,  metric = "roc_auc"),   "multinom_glmnet")
tree_test <- finalize_and_test(tree_wf, select_best(tree_res, metric = "roc_auc"),   "decision_tree")
rf_test   <- finalize_and_test(rf_wf,   select_best(rf_res,   metric = "roc_auc"),   "random_forest")
xgb_test  <- if (xgb_available) finalize_and_test(xgb_wf,
                                                 select_best(xgb_res, metric = "roc_auc"),
                                                 "xgboost") else NULL

ml_test_metrics <- bind_rows(
  null_test$metrics, glm_test$metrics, tree_test$metrics,
  rf_test$metrics,
  if (!is.null(xgb_test)) xgb_test$metrics
)
saveRDS(ml_test_metrics, file.path(tab_dir, "ml_test_metrics.rds"))

best_fit <- switch(best_name,
                   "null"            = null_test$fit,
                   "multinom_glmnet" = glm_test$fit,
                   "decision_tree"   = tree_test$fit,
                   "random_forest"   = rf_test$fit,
                   "xgboost"         = xgb_test$fit)
best_preds <- collect_predictions(best_fit)

## ---- confusion_matrix --------
cm <- conf_mat(best_preds, truth = SampleType, estimate = .pred_class)
saveRDS(list(model = best_name, cm = cm),
        file.path(tab_dir, "ml_confusion_matrix.rds"))
ggsave(file.path(fig_dir, "ml_confusion_matrix.png"),
       autoplot(cm, type = "heatmap") +
         labs(title = sprintf("Test-set confusion matrix - %s", best_name)),
       width = 7, height = 6, dpi = 300)

## ---- roc_curves --------
# One-vs-rest ROC for the winning model on the held-out test set.
prob_cols <- grep("^\\.pred_", names(best_preds), value = TRUE)
prob_cols <- setdiff(prob_cols, ".pred_class")
roc_df    <- best_preds |>
  roc_curve(truth = SampleType, all_of(prob_cols))

p_roc <- ggplot(roc_df,
                aes(x = 1 - specificity, y = sensitivity, color = .level)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  scale_color_manual(values = colors_body) +
  coord_equal() +
  labs(title = sprintf("Test-set ROC curves - %s (one-vs-rest)", best_name),
       x = "False positive rate", y = "True positive rate",
       color = "Body Site")
ggsave(file.path(fig_dir, "ml_roc_curves.png"),
       p_roc, width = 7, height = 6, dpi = 300)

## ---- bootstrap_uncertainty --------
# Bootstrap 95% CIs for accuracy and macro-AUC on the held-out test set.
B <- 1000
set.seed(123)
n  <- nrow(best_preds)
boot_metrics <- vapply(seq_len(B), function(b) {
  idx <- sample.int(n, n, replace = TRUE)
  s   <- best_preds[idx, , drop = FALSE]
  acc <- mean(s$.pred_class == s$SampleType)
  auc <- tryCatch({
    yardstick::roc_auc_vec(truth    = s$SampleType,
                           estimate = as.matrix(s[, prob_cols]))
  }, error = function(e) NA_real_)
  c(accuracy = acc, roc_auc = auc)
}, numeric(2))

boot_ci <- tibble(
  metric = c("accuracy", "roc_auc"),
  point  = c(mean(best_preds$.pred_class == best_preds$SampleType),
             yardstick::roc_auc_vec(truth = best_preds$SampleType,
                                    estimate = as.matrix(best_preds[, prob_cols]))),
  lower  = c(quantile(boot_metrics["accuracy", ], 0.025, na.rm = TRUE),
             quantile(boot_metrics["roc_auc",  ], 0.025, na.rm = TRUE)),
  upper  = c(quantile(boot_metrics["accuracy", ], 0.975, na.rm = TRUE),
             quantile(boot_metrics["roc_auc",  ], 0.975, na.rm = TRUE))
) |> mutate(model = best_name, n_boot = B)
saveRDS(boot_ci, file.path(tab_dir, "ml_test_uncertainty.rds"))
cat(sprintf(
  "[H4] Bootstrap 95%% CI on the test set (%d resamples):\n", B
))
cat(sprintf("       accuracy = %.3f  [%.3f, %.3f]\n",
            boot_ci$point[1], boot_ci$lower[1], boot_ci$upper[1]))
cat(sprintf("       roc_auc  = %.3f  [%.3f, %.3f]\n",
            boot_ci$point[2], boot_ci$lower[2], boot_ci$upper[2]))

## ---- variable_importance --------
# Top-20 features for the winning model. Uses vip if available; otherwise
# falls back to NULL and skips the figure.
imp_df <- tryCatch({
  if (requireNamespace("vip", quietly = TRUE)) {
    fit_obj <- extract_fit_parsnip(best_fit$.workflow[[1]])
    vip::vi(fit_obj) |>
      arrange(desc(Importance)) |>
      slice_head(n = 20)
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(imp_df) && nrow(imp_df) > 0) {
  saveRDS(imp_df, file.path(tab_dir, "ml_variable_importance.rds"))
  p_imp <- ggplot(imp_df,
                  aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "steelblue", alpha = 0.85) +
    coord_flip() +
    labs(title = sprintf("Top 20 features - %s", best_name),
         subtitle = "Higher = more useful for predicting body site",
         x = NULL, y = "Importance")
  ggsave(file.path(fig_dir, "ml_variable_importance.png"),
         p_imp, width = 8, height = 7, dpi = 300)
}

cat(sprintf("[H4] Test-set point estimates for %s:\n", best_name))
print(ml_test_metrics |> filter(model == best_name))

## ---- finish --------
cat("Statistical analysis complete. Outputs written to:\n")
cat("  ", fig_dir, "\n")
cat("  ", tab_dir, "\n")
cat("  ", out_dir, "\n")
