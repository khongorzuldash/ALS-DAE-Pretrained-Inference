# ALS DAE Pretrained Inference Demo

Inference-only demonstration using a pretrained denoising autoencoder.
No model training is performed.

## Contents
- model/FULL_HVG3000_DAE_best.keras: pretrained DAE
- model/gene_min.npy: original normalization minimum
- model/gene_range.npy: original normalization range
- data/sample_logNorm_100cells.npz: 100 example cells
- data/HVG3000_genes.csv: ordered 3000 input genes
- test.py: reconstruction error and lowest-20-percent selection

## Run
Install: pip install -r requirements.txt
Execute: python test.py

## Output
- reconstruction_errors.csv
- selected_cell_indices.csv

## Notes
Input genes must match the original HVG3000 order.
The input is normalized using saved training parameters.
This 100-cell example demonstrates inference, not independent validation.
Original dataset: GEO GSE244263.
Only load trusted model files.

## DEPA reproducibility (ALS progression)

DEPA analyzes 1,272 genes across 40 samples.

### Original pipeline

1. Gene-wise z-score normalization.
2. DTW with absolute-distance cost and window=5.
3. Two-dimensional metric MDS, random_state=42.
4. K-means: K=6, random_state=42, n_init=200.

### Verification

Run: `python test_depa.py`

Expected results:
- Genes: 1,272
- DTW pairs tested: 45
- DTW maximum error: approximately 1.93e-6
- K-means ARI: 1.0
- Return code: 0

The test recalculates DTW for 45 gene pairs and reproduces
K-means using the archived MDS coordinates.

It does not recompute the complete DTW matrix or MDS embedding.
Independent MDS recomputation produced a pairwise-distance
correlation of 0.999818 and clustering ARI of 0.989757.

### Data files

The `data/depa/` directory contains the gene-wise z-score input,
DTW distance matrix, MDS coordinates, and K=6 cluster assignments.

The 100-cell pretrained DAE inference example is separate
and does not regenerate the 1,272-gene DEPA input.

## Running the complete DEPA pipeline

The repository provides two DEPA execution modes:

### 1. Verify archived results (recommended first)

```bash
python test_depa.py
```

This verifies 45 DTW gene pairs and reproduces the archived
K=6 cluster memberships (expected ARI = 1.0).

### 2. Recompute DTW, MDS and K-means

```bash
python run_depa_full.py
```

Input: `data/depa/DTW_input_P005_gene_by_sample_zscore.csv`

Output directory: `results/depa_full/`

Generated files:

- `DTW_distance_matrix.csv`
- `MDS_coordinates.csv`
- `K6_gene_clusters.csv`

The complete calculation may be computationally expensive.
The full pipeline was successfully executed on all 1,272 genes.
See the end-to-end verification results below.

MDS results may differ slightly across scikit-learn versions.
Consequently, the regenerated K=6 labels may not exactly
match the archived cluster memberships.

### Pretrained DAE inference

```bash
python test.py
```

This uses the supplied pretrained DAE model and sample data.
Model training is not required.

## Full DEPA end-to-end verification (1,272 genes)

The complete DTW-MDS-K-means (K=6) workflow was executed
successfully on all 1,272 candidate genes.

- Execution time: 8.89 minutes
- DTW maximum absolute error vs. archived matrix: 0.0
- DTW mean absolute error vs. archived matrix: 0.0
- MDS pairwise-distance correlation: 0.9998177859
- K-means adjusted Rand index (ARI): 0.98975674
- K-means normalized mutual information (NMI): 0.98545060
- Matched cluster assignments after label alignment: 1,266/1,272 (99.53%)

The complete pipeline executed successfully, but the
recomputed MDS and K-means results are not numerically
identical to the archived outputs.

The original DTW distance matrix was reproduced exactly.
