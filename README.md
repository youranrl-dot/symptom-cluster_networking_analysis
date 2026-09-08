# symptom-cluster_networking_analysis

R pipeline for building a symptom co-occurrence network from binary
symptom-presence data and identifying symptom clusters via Louvain
community detection. Built for a colorectal cancer (CRC) symptom-extraction
study (LLM-derived symptom labels from MIMIC-IV discharge notes), but the
script only assumes a patient-level (or note-level) table of `pred_*`
binary symptom columns plus a `subject_id` column, so it generalizes to any
binary symptom dataset.

## What it does

1. Aggregates note-level symptom predictions to patient level using
   "any-note" logic (a symptom counts as present if it was flagged in any
   of that patient's notes).
2. Filters to symptoms with prevalence >= 5% (configurable) to avoid
   unstable correlation estimates from very rare symptoms.
3. Builds a phi-correlation matrix across the retained symptoms (phi = the
   Pearson correlation of two binary 0/1 variables).
4. Constructs a symptom co-occurrence network (nodes = symptoms, edges =
   phi >= threshold, default 0.10).
5. Runs Louvain community detection (`igraph::cluster_louvain`) to identify
   symptom clusters, and computes strength centrality per symptom.
6. Runs sensitivity/stability checks: edge-threshold sensitivity (0.10 vs
   0.15), Louvain seed stability (100 seeds), and bootstrap stability
   (200 patient resamples), each summarized with the Adjusted Rand Index
   (`mclust::adjustedRandIndex`) against the main clustering.
7. Writes per-patient cluster "burden" scores (count of present symptoms
   per cluster) for downstream outcome modeling (e.g., logistic regression
   of 1-year mortality / readmission on cluster burden).

## Requirements

- R >= 4.0
- Packages: `igraph`, `mclust` (only these two; no tidyverse/qgraph
  dependency, so it runs anywhere R is installed without extra setup)

The script installs missing packages automatically. On Windows it forces
binary installs (`type = "binary"`) so it doesn't try to compile from
source, which requires Rtools/gfortran that most machines don't have.

## Usage

```r
# 1. Edit the setwd() call near the top of the script to point at the
#    folder containing your input CSV.
# 2. Make sure your CSV has:
#      - one column named subject_id (or a note_id like "<subject_id>-DS-#"
#        that the script parses)
#      - symptom columns prefixed pred_ (0/1 binary presence/absence)
# 3. Run the whole script (e.g. Source in RStudio, or Rscript from a shell):
Rscript symptom_network_analysis.R
```

Outputs written to the working directory:
- `symptom_cluster_centrality_table.csv` — symptom, cluster assignment,
  strength centrality, prevalence
- `symptom_network_plot.pdf` — network plot colored by cluster
- `patient_level_symptom_clusters.csv` — one row per patient with symptom
  columns plus a `clusterN_burden` score per identified cluster

## Data

No patient-level data is included in this repository. MIMIC-IV is a
credentialed, restricted-access dataset (PhysioNet), and its data use
agreement prohibits redistributing note content or patient-level derived
data outside PhysioNet-approved environments. Bring your own binary
symptom-presence table to use this script.
