# symptom-cluster_networking_analysis

R pipeline for building a symptom co-occurrence network from binary
symptom-presence data and identifying symptom clusters via Louvain
community detection. Built for a colorectal cancer (CRC) symptom-extraction
study comparing two LLM extraction methods (Gemini 3.5 Flash and Claude
Haiku) on MIMIC-IV discharge notes, but each script only assumes a
patient-level (or note-level) table of `pred_*` binary symptom columns
plus a `subject_id` column, so it generalizes to any binary symptom
dataset.

## Scripts

- `symptom_network_analysis.R` — main pipeline (Gemini 3.5 Flash
  extraction), described below.
- `symptom_network_analysis_claude.R` — identical pipeline run on the
  Claude Haiku extraction, plus a Section 12 that computes the
  cross-model Adjusted Rand Index between the two models' Louvain
  partitions (a robustness check: are the clusters an artifact of one
  extraction model, or do they reproduce across models?).
- `symptom_cluster_outcome_regression.R` — logistic regression of three
  clinical outcomes (in-hospital mortality, 30-day readmission, 1-year
  mortality) on the per-cluster symptom burden scores written by the
  network script, adjusted for age and sex.

## What the network script does

1. Aggregates note-level symptom predictions to patient level using
   "any-note" logic (a symptom counts as present if it was flagged in any
   of that patient's notes).
2. Excludes patients whose notes all failed extraction (no detectable
   CC/HPI section) rather than letting a missing-input note masquerade as
   a genuine "symptom absent" — see the note at the top of each script.
3. Filters to symptoms with prevalence >= 5% (configurable) to avoid
   unstable correlation estimates from very rare symptoms.
4. Builds a phi-correlation matrix across the retained symptoms (phi = the
   Pearson correlation of two binary 0/1 variables).
5. Constructs a symptom co-occurrence network (nodes = symptoms, edges =
   phi >= threshold, default 0.10) and plots it with cluster-colored
   nodes and a legend.
6. Runs Louvain community detection (`igraph::cluster_louvain`) to identify
   symptom clusters, and computes strength centrality per symptom.
7. Runs sensitivity/stability checks: edge-threshold sensitivity (0.10 vs
   0.15), Louvain seed stability (100 seeds), and bootstrap stability
   (200 patient resamples), each summarized with the Adjusted Rand Index
   (`mclust::adjustedRandIndex`) against the main clustering.
8. Writes per-patient cluster "burden" scores (count of present symptoms
   per cluster) for downstream outcome modeling.

## Requirements

- R >= 4.0
- Packages: `igraph`, `mclust` (only these two; no tidyverse/qgraph
  dependency, so it runs anywhere R is installed without extra setup)

The scripts install missing packages automatically. On Windows they force
binary installs (`type = "binary"`) so they don't try to compile from
source, which requires Rtools/gfortran that most machines don't have.

## Usage

```r
# 1. Edit the setwd() call near the top of each script to point at the
#    folder containing your input CSV(s).
# 2. Make sure your CSV has:
#      - one column named subject_id (or a note_id like "<subject_id>-DS-#"
#        that the script parses)
#      - symptom columns prefixed pred_ (0/1 binary presence/absence)
# 3. Run the network script first (Source in RStudio, or Rscript from a shell):
Rscript symptom_network_analysis.R
# 4. Optionally run the Claude-extraction counterpart for a cross-model
#    comparison (needs symptom_cluster_centrality_table.csv from step 3
#    in the same folder):
Rscript symptom_network_analysis_claude.R
# 5. Run the outcome regression (needs patient_level_symptom_clusters.csv
#    from step 3, plus an outcomes CSV with subject_id/age/gender/outcome
#    columns):
Rscript symptom_cluster_outcome_regression.R
```

Outputs written to the working directory (network script):
- `symptom_cluster_centrality_table.csv` — symptom, cluster assignment,
  strength centrality, prevalence
- `patient_level_symptom_clusters.csv` — one row per patient with symptom
  columns plus a `clusterN_burden` score per identified cluster

The Claude-extraction script additionally writes
`symptom_cluster_centrality_table_claude.csv`,
`patient_level_symptom_clusters_claude.csv`, and (when the Gemini table is
present) `symptom_cluster_gemini_vs_claude_comparison.csv`.

The regression script writes `section37_cluster_burden_results_N*.csv`
with OR / 95% CI / p-value per cluster per outcome.

## Data

No patient-level data is included in this repository. MIMIC-IV is a
credentialed, restricted-access dataset (PhysioNet), and its data use
agreement prohibits redistributing note content or patient-level derived
data outside PhysioNet-approved environments. Bring your own binary
symptom-presence table to use these scripts.
