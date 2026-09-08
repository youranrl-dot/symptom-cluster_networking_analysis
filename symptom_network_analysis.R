## ============================================================
## CRC Symptom Co-occurrence Network + Louvain Community Detection
## Replicates the manuscript's Section 3.5/3.6 / Table 4 / Figure 3-4 pipeline
## (originally built in Python: networkx + python-louvain 0.16)
##
## Input: crc_symptom_predictions_with_demographics_survival.csv
##        (note-level: 2,728 rows = notes, 1,507 unique subject_id)
##
## Only base R + igraph + mclust are used (no tidyverse/qgraph dependency),
## so this runs anywhere R is installed without extra setup.
## ============================================================

## ---- 0. Packages ----------------------------------------------------------
## Force pre-compiled binaries so Windows doesn't try to build from source
## (that needs Rtools/gfortran, which most machines don't have installed).
if (.Platform$OS.type == "windows") {
  options(install.packages.compile.from.source = "never")
}
required_pkgs <- c("igraph", "mclust")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, type = "binary")

library(igraph)   # network + cluster_louvain
library(mclust)   # adjustedRandIndex

set.seed(42)  # matches manuscript's Louvain seed=42

## Working directory -- set this to the folder that holds your input CSV.
## All outputs (plot, centrality table, patient-level cluster file) are also
## written here. Example (Windows): setwd("C:/path/to/your/folder")
# setwd("path/to/your/folder")

## ---- 1. Load + aggregate to patient level ---------------------------------
## The merged file is note-level (multiple discharge notes per patient).
## Main analysis uses "any-note" logic: a symptom counts as present for a
## patient if it was flagged PRESENT in ANY of that patient's notes -- this
## reproduces the 20/46 symptoms retained at >=5% and the exact prevalence
## counts already validated against your data.

df <- read.csv("crc_symptom_predictions_with_demographics_survival.csv",
                stringsAsFactors = FALSE)

symptom_cols <- grep("^pred_", names(df), value = TRUE)

patient_level <- aggregate(df[symptom_cols],
                            by = list(subject_id = df$subject_id),
                            FUN = function(x) as.integer(max(x, na.rm = TRUE)))

n_patients <- nrow(patient_level)
cat("Patients:", n_patients, "\n")  # should be 1507

## ---- 2. Filter to symptoms with prevalence >= 5% --------------------------
prev <- colMeans(patient_level[symptom_cols]) * 100
keep_syms <- names(prev[prev >= 5])
cat(length(keep_syms), "symptoms retained at >=5% prevalence\n")
print(sort(round(prev[keep_syms], 1), decreasing = TRUE))

sym_mat <- as.matrix(patient_level[keep_syms])
colnames(sym_mat) <- sub("^pred_", "", colnames(sym_mat))

## ---- 3. Phi correlation matrix ---------------------------------------------
## phi coefficient = Pearson correlation for two binary (0/1) variables;
## this reproduces the "Symptom Phi Correlation Matrix" in your Figure 3.
phi_mat <- cor(sym_mat, method = "pearson")
diag(phi_mat) <- 0

cat("Phi range (off-diagonal):",
    round(min(phi_mat), 3), "to", round(max(phi_mat), 3), "\n")

## ---- 4. Build network (edges = phi >= threshold) --------------------------
build_network <- function(phi_mat, threshold = 0.10) {
  adj <- phi_mat
  adj[adj < threshold] <- 0
  graph_from_adjacency_matrix(adj, mode = "undirected",
                               weighted = TRUE, diag = FALSE)
}

g_main <- build_network(phi_mat, threshold = 0.10)
cat("Edges at threshold 0.10:", ecount(g_main), "\n")

isolates <- V(g_main)$name[degree(g_main) == 0]
cat("Isolated symptoms (no edge >=0.10):",
    if (length(isolates) == 0) "none" else paste(isolates, collapse = ", "), "\n")

## ---- 5. Louvain community detection ---------------------------------------
set.seed(42)
louvain_main <- cluster_louvain(g_main, weights = E(g_main)$weight)

cluster_membership <- data.frame(
  symptom = V(g_main)$name,
  cluster = as.integer(membership(louvain_main))
)
cluster_membership <- cluster_membership[order(cluster_membership$cluster,
                                                cluster_membership$symptom), ]

cat("\nNumber of clusters:", length(unique(cluster_membership$cluster)), "\n")
print(table(cluster_membership$cluster))
print(cluster_membership, row.names = FALSE)

cat("Modularity:", round(modularity(louvain_main), 3), "\n")

## ---- 6. Strength centrality (Table 4) --------------------------------------
strength_centrality <- strength(g_main, weights = E(g_main)$weight)
centrality_table <- cluster_membership
centrality_table$strength <- round(strength_centrality[centrality_table$symptom], 2)
centrality_table$prevalence_pct <- round(prev[paste0("pred_", centrality_table$symptom)], 1)
centrality_table <- centrality_table[order(centrality_table$cluster,
                                            -centrality_table$strength), ]
print(centrality_table, row.names = FALSE)

write.csv(centrality_table, "symptom_cluster_centrality_table.csv", row.names = FALSE)

## ---- 7. Network plot (base igraph; swap in qgraph if you have it) ---------
palette <- c("#F28E8E", "#8EC6F2", "#9ED9A0", "#D6B8F0", "#F2D98E")
node_colors <- palette[(as.integer(membership(louvain_main)[V(g_main)$name]) - 1) %% length(palette) + 1]

pdf("symptom_network_plot.pdf", width = 9, height = 7)
set.seed(42)
plot(g_main,
     vertex.color = node_colors,
     vertex.size = 14,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     edge.width = E(g_main)$weight * 8,
     layout = layout_with_fr(g_main),
     main = "CRC Symptom Co-occurrence Network (phi >= 0.10)")
dev.off()

## ---- 8. Threshold sensitivity: 0.10 vs 0.15 --------------------------------
g_015 <- build_network(phi_mat, threshold = 0.15)
set.seed(42)
louvain_015 <- cluster_louvain(g_015, weights = E(g_015)$weight)

ari_thresh <- adjustedRandIndex(
  as.integer(membership(louvain_main)[V(g_main)$name]),
  as.integer(membership(louvain_015)[V(g_main)$name])
)
cat("\nThreshold sensitivity (0.10 vs 0.15): ARI =", round(ari_thresh, 3),
    "| n clusters:", length(unique(membership(louvain_main))), "vs",
    length(unique(membership(louvain_015))), "\n")

## ---- 9. Seed stability: Louvain across 100 seeds ---------------------------
## igraph's Louvain has minor run-to-run variability from tie-breaking order;
## this checks how stable the partition is (your manuscript's "100 seeds" check).
n_seeds <- 100
partitions <- vector("list", n_seeds)
for (i in seq_len(n_seeds)) {
  set.seed(i)
  partitions[[i]] <- as.integer(membership(cluster_louvain(g_main, weights = E(g_main)$weight))[V(g_main)$name])
}

reference <- partitions[[1]]
ari_vals <- sapply(partitions, function(p) adjustedRandIndex(reference, p))
cat("\nLouvain 100-seed ARI (vs seed 1): mean =", round(mean(ari_vals), 3),
    "range =", round(min(ari_vals), 3), "-", round(max(ari_vals), 3), "\n")

## ---- 10. Bootstrap stability (resample patients) ---------------------------
n_boot <- 200
boot_ari <- rep(NA_real_, n_boot)
ref_membership <- as.integer(membership(louvain_main)[V(g_main)$name])
names(ref_membership) <- V(g_main)$name

for (b in seq_len(n_boot)) {
  set.seed(1000 + b)
  boot_idx <- sample(seq_len(n_patients), n_patients, replace = TRUE)
  boot_mat <- sym_mat[boot_idx, , drop = FALSE]

  boot_phi <- suppressWarnings(cor(boot_mat, method = "pearson"))
  boot_phi[is.na(boot_phi)] <- 0
  diag(boot_phi) <- 0

  g_boot <- build_network(boot_phi, threshold = 0.10)
  common_nodes <- intersect(V(g_main)$name, V(g_boot)$name)
  if (length(common_nodes) < 3) next

  g_boot_sub <- induced_subgraph(g_boot, common_nodes)
  clust_boot <- tryCatch(
    as.integer(membership(cluster_louvain(g_boot_sub, weights = E(g_boot_sub)$weight))),
    error = function(e) NULL)
  if (is.null(clust_boot)) next
  names(clust_boot) <- V(g_boot_sub)$name

  boot_ari[b] <- adjustedRandIndex(ref_membership[common_nodes], clust_boot[common_nodes])
}

cat("\nBootstrap ARI (", n_boot, "resamples): mean =", round(mean(boot_ari, na.rm = TRUE), 3),
    "SD =", round(sd(boot_ari, na.rm = TRUE), 3),
    "median =", round(median(boot_ari, na.rm = TRUE), 3),
    "range =", round(min(boot_ari, na.rm = TRUE), 3), "-", round(max(boot_ari, na.rm = TRUE), 3), "\n")
cat("Valid bootstrap replicates:", sum(!is.na(boot_ari)), "/", n_boot, "\n")

## ---- 11. Save cluster-burden scores per patient (for outcome models) -------
## e.g. for logistic regression of 1-year mortality on cluster burden, as in
## your Section 3.7 (cluster burden = count of present symptoms in that cluster)
for (cl in unique(cluster_membership$cluster)) {
  syms_in_cluster <- cluster_membership$symptom[cluster_membership$cluster == cl]
  cols_in_cluster <- paste0("pred_", syms_in_cluster)
  patient_level[[paste0("cluster", cl, "_burden")]] <-
    rowSums(patient_level[, cols_in_cluster, drop = FALSE])
}

write.csv(patient_level, "patient_level_symptom_clusters.csv", row.names = FALSE)

cat("\nDone. Outputs written:\n",
    "- symptom_cluster_centrality_table.csv\n",
    "- symptom_network_plot.pdf\n",
    "- patient_level_symptom_clusters.csv (includes cluster burden scores per patient)\n")
