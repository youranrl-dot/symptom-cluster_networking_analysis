## ============================================================
## CRC Symptom Co-occurrence Network + Louvain Community Detection
## -- CLAUDE HAIKU extraction version --
## Mirrors symptom_network_analysis.R (the Gemini 3.5 Flash version) exactly,
## but runs on predictions_claude_haiku_full.csv, so the two model's cluster
## solutions can be compared head-to-head (Section 12 below).
##
## Input: crc_symptom_predictions_claude_with_demographics_survival.csv
##        (note-level: 2,728 rows = notes, 1,507 unique subject_id --
##         identical cohort/notes to the Gemini file, same demographics)
##
## Only base R + igraph + mclust are used.
## ============================================================

## ---- 0. Packages ----------------------------------------------------------
if (.Platform$OS.type == "windows") {
  options(install.packages.compile.from.source = "never")
}
required_pkgs <- c("igraph", "mclust")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, type = "binary")

library(igraph)   # network + cluster_louvain
library(mclust)   # adjustedRandIndex

set.seed(42)  # same seed as the Gemini run, for a fair comparison

## Working directory -- same folder you used for the Gemini script. Put
## crc_symptom_predictions_claude_with_demographics_survival.csv here too.
# setwd("path/to/your/folder")

## ---- 1. Load + aggregate to patient level ---------------------------------
df <- read.csv("crc_symptom_predictions_claude_with_demographics_survival.csv",
                stringsAsFactors = FALSE)

symptom_cols <- grep("^pred_", names(df), value = TRUE)

## Claude Haiku failed to extract from 18/2728 notes (error = "missing_cc_hpi":
## no CC/HPI section found in the note text) -- Gemini's run has zero such
## failures on this same note set. safe_max() ignores NA notes for a patient
## who has at least one successfully-extracted note; a patient whose notes
## ALL failed gets NA (handled just below), instead of the -Inf/NA-by-accident
## you'd get from max(x, na.rm=TRUE) on an all-NA vector.
safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_integer_)
  as.integer(max(x))
}

patient_level <- aggregate(df[symptom_cols],
                            by = list(subject_id = df$subject_id),
                            FUN = safe_max)

n_before <- nrow(patient_level)

## Exclude patients whose notes ALL failed extraction (no usable Claude data
## at all for that patient) -- documented here rather than silently dropped,
## and rather than imputing "symptom absent" for a note that was never read.
all_na_patient <- apply(patient_level[symptom_cols], 1, function(x) all(is.na(x)))
excluded_subject_ids <- patient_level$subject_id[all_na_patient]
cat(sum(all_na_patient),
    "patient(s) excluded from the Claude network analysis (all notes failed",
    "extraction, error = 'missing_cc_hpi'):\n")
print(excluded_subject_ids)
patient_level <- patient_level[!all_na_patient, ]

n_patients <- nrow(patient_level)
cat("Patients retained for Claude network analysis:", n_patients,
    "of", n_before, "(Gemini analysis used all 1507)\n")

## ---- 2. Filter to symptoms with prevalence >= 5% --------------------------
prev <- colMeans(patient_level[symptom_cols], na.rm = TRUE) * 100
keep_syms <- names(prev[prev >= 5])
cat(length(keep_syms), "symptoms retained at >=5% prevalence (Claude)\n")
print(sort(round(prev[keep_syms], 1), decreasing = TRUE))

sym_mat <- as.matrix(patient_level[keep_syms])
colnames(sym_mat) <- sub("^pred_", "", colnames(sym_mat))

## ---- 3. Phi correlation matrix ---------------------------------------------
phi_mat <- cor(sym_mat, method = "pearson", use = "pairwise.complete.obs")
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

cat("\nNumber of clusters (Claude):", length(unique(cluster_membership$cluster)), "\n")
print(table(cluster_membership$cluster))
print(cluster_membership, row.names = FALSE)

cat("Modularity:", round(modularity(louvain_main), 3), "\n")

## NOTE: Louvain cluster IDs (1, 2, 3, ...) are arbitrary labels assigned by
## the algorithm -- they do NOT necessarily correspond to the Gemini run's
## "Systemic / Gastrointestinal / CRC Disease-Specific" numbering. Compare
## membership by SYMPTOM COMPOSITION (Section 12), not by cluster number.

## ---- 6. Strength centrality --------------------------------------
strength_centrality <- strength(g_main, weights = E(g_main)$weight)
centrality_table <- cluster_membership
centrality_table$strength <- round(strength_centrality[centrality_table$symptom], 2)
centrality_table$prevalence_pct <- round(prev[paste0("pred_", centrality_table$symptom)], 1)
centrality_table <- centrality_table[order(centrality_table$cluster,
                                            -centrality_table$strength), ]
print(centrality_table, row.names = FALSE)

write.csv(centrality_table, "symptom_cluster_centrality_table_claude.csv", row.names = FALSE)

## ---- 7. Network plot (base igraph plot window, with legend) ---------------
palette <- c("#F28E8E", "#8EC6F2", "#9ED9A0", "#D6B8F0", "#F2D98E")
mem <- as.integer(membership(louvain_main)[V(g_main)$name])
node_colors <- palette[(mem - 1) %% length(palette) + 1]

## Generic labels -- see the note in Section 5. Re-label these by hand once
## you've looked at Section 12's symptom-composition comparison to Gemini.
n_clusters <- sort(unique(mem))
legend_labels <- paste("Cluster", n_clusters)
legend_colors <- palette[(n_clusters - 1) %% length(palette) + 1]

set.seed(42)
lay <- layout_with_fr(g_main, niter = 2000)

draw_network <- function() {
  par(mar = c(1, 1, 3, 1))
  plot(g_main,
       vertex.color = node_colors,
       vertex.size = 10,
       vertex.label.cex = 0.75,
       vertex.label.color = "black",
       vertex.label.dist = 1.1,
       vertex.label.degree = -pi / 2,
       edge.width = E(g_main)$weight * 8,
       edge.color = adjustcolor("gray50", alpha.f = 0.5),
       layout = lay,
       main = "CRC Symptom Co-occurrence Network -- Claude Haiku (phi >= 0.10)")
  legend("bottomleft", legend = legend_labels, pt.bg = legend_colors,
         pch = 21, pt.cex = 2, cex = 0.85, bty = "n")
}
draw_network()

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
for (cl in unique(cluster_membership$cluster)) {
  syms_in_cluster <- cluster_membership$symptom[cluster_membership$cluster == cl]
  cols_in_cluster <- paste0("pred_", syms_in_cluster)
  patient_level[[paste0("cluster", cl, "_burden")]] <-
    rowSums(patient_level[, cols_in_cluster, drop = FALSE])
}

write.csv(patient_level, "patient_level_symptom_clusters_claude.csv", row.names = FALSE)

## ---- 12. Cross-model comparison: Claude vs. Gemini cluster structure ------
## Requires symptom_cluster_centrality_table.csv (the GEMINI run's output,
## from symptom_network_analysis.R Section 6) to already exist in this same
## working directory. This is the formal "which model's clusters are more
## valid / are the clusters robust across extraction models" comparison.

gemini_path <- "symptom_cluster_centrality_table.csv"
if (!file.exists(gemini_path)) {
  cat("\n[Section 12 skipped] '", gemini_path, "' not found in the working ",
      "directory -- run symptom_network_analysis.R (Gemini version) first ",
      "so its centrality table is written here.\n", sep = "")
} else {
  gemini_membership <- read.csv(gemini_path, stringsAsFactors = FALSE)[, c("symptom", "cluster")]
  names(gemini_membership) <- c("symptom", "cluster_gemini")

  claude_membership <- cluster_membership
  names(claude_membership) <- c("symptom", "cluster_claude")

  common <- merge(gemini_membership, claude_membership, by = "symptom")
  cat("\nSymptoms retained (>=5%) by both models:", nrow(common), "\n")
  cat("Retained by Gemini only:",
      paste(setdiff(gemini_membership$symptom, claude_membership$symptom), collapse = ", "), "\n")
  cat("Retained by Claude only:",
      paste(setdiff(claude_membership$symptom, gemini_membership$symptom), collapse = ", "), "\n")

  cat("\nCross-tab of cluster assignment (rows = Gemini cluster, cols = Claude cluster):\n")
  print(table(Gemini = common$cluster_gemini, Claude = common$cluster_claude))

  ari_cross <- adjustedRandIndex(common$cluster_gemini, common$cluster_claude)
  cat("\nGemini vs. Claude cluster concordance (ARI, on", nrow(common),
      "shared symptoms):", round(ari_cross, 3), "\n")
  cat("(ARI = 1 means identical partitions; ARI ~ 0 means no more agreement\n",
      " than chance. This is the number to report as evidence that the\n",
      " symptom cluster structure is -- or isn't -- robust across LLM\n",
      " extraction methods.)\n", sep = "")

  write.csv(common, "symptom_cluster_gemini_vs_claude_comparison.csv", row.names = FALSE)
}

cat("\nDone. Outputs written (Claude):\n",
    "- symptom_cluster_centrality_table_claude.csv\n",
    "- patient_level_symptom_clusters_claude.csv\n",
    "- symptom_cluster_gemini_vs_claude_comparison.csv (if Gemini table was found)\n")
