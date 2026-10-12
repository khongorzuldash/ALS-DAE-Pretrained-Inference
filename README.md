# ALS DAE–DEPA Analysis Pipeline

## Research code and example data

This repository provides the computational code used in our amyotrophic lateral sclerosis (ALS) study, together with a pretrained denoising autoencoder (DAE) model, supporting files, and example data for demonstrating the analysis workflow.

The original study was conducted using the full discovery dataset (GEO: **GSE244263**), comprising **458,226 cells**. Cells were ranked according to their DAE reconstruction errors, and **91,646 cells** with the lowest reconstruction errors (20%) were retained for the primary downstream analysis.

The example data included in this repository are provided to demonstrate code execution and do not replace the full dataset used in the original study.

## Analysis overview

The study workflow comprises two main components:

1. **DAE-based cell prioritization:** A DAE was trained on the full discovery dataset to calculate cell-level reconstruction errors. Cells with the lowest errors were prioritized for downstream discovery analysis. Alternative cell-selection thresholds were evaluated in sensitivity analyses.
2. **Dynamic Expression Pattern Analysis (DEPA):** Candidate gene-expression profiles were analyzed using dynamic time warping (DTW), multidimensional scaling (MDS), and K-means clustering to identify progression-associated transcriptional programs.

## Repository contents

### Pretrained DAE inference example

| File | Description |
| --- | --- |
| `model/FULL_HVG3000_DAE_best.keras` | Pretrained DAE model supplied for the inference example |
| `model/gene_min.npy` | Saved gene-wise normalization parameter |
| `model/gene_range.npy` | Saved gene-wise normalization parameter |
| `data/HVG3000_genes.csv` | Ordered input-gene list for the pretrained model |
| `data/sample_logNorm_100cells.npz` | Example expression data for 100 cells |
| `test.py` | Loads the pretrained model, calculates reconstruction errors, and demonstrates low-error cell selection |

The pretrained `.keras` model is used for the supplied inference example; running `test.py` does **not** retrain the DAE. The 100-cell example is for demonstrating execution and is not the full discovery dataset.

### DEPA analysis

| File or directory | Description |
| --- | --- |
| `run_depa_full.py` | DTW–MDS–K-means analysis pipeline |
| `test_depa.py` | Utility for checking the supplied DEPA files and workflow |
| `data/depa/` | Processed expression matrix and reference files for DEPA |

The DEPA input contains **1,272 candidate genes across 40 samples**. It is provided as a processed input for the gene-level analysis and is separate from the 100-cell DAE inference example.

## Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

## Run the pretrained DAE example

```bash
python test.py
```

The script loads the supplied pretrained model and example data, calculates per-cell reconstruction errors, and demonstrates selection of the lowest-error 20% of example cells.

Expected output files:

- `reconstruction_errors.csv`
- `selected_cell_indices.csv`

**Important:** Input genes must follow the order specified in `data/HVG3000_genes.csv`. Use the normalization files provided with the pretrained model. The selected cells from this example are not the 91,646 cells reported in the original full-data analysis.

## Run DEPA

The supplied DEPA pipeline analyzes the processed **1,272-gene × 40-sample** expression matrix using DTW, MDS, and K-means clustering (**K = 6**).

Run the analysis:

```bash
python run_depa_full.py
```

Input file:

```text
data/depa/DTW_input_P005_gene_by_sample_zscore.csv
```

Output directory:

```text
results/depa_full/
```

Expected output files:

- `DTW_distance_matrix.csv`
- `MDS_coordinates.csv`
- `K6_gene_clusters.csv`

The full DTW–MDS–K-means calculation may require substantial computing time. The repository also includes the corresponding reference files.

An optional utility can be run to check the supplied DEPA files:

```bash
python test_depa.py
```

## Full-data DAE methodology

In the original full-data analysis, the DAE was trained on a cell-by-gene expression matrix constructed from the discovery dataset using **3,000 highly variable genes**. The full-data DAE training implementation used **PyTorch**, with an encoder of **1,024 and 256 hidden units** and a **64-dimensional latent representation**. Training used **20% dropout**, **15% input masking noise**, the **Adam** optimizer, and **mean squared error** reconstruction loss. The training settings were **25 epochs**, **batch size 512**, **learning rate 0.001**, and **weight decay 0.00001**.

After training, reconstruction errors were calculated using unmasked input profiles. Cells were ranked by reconstruction error, and the lowest-error 20% were selected for primary discovery analysis. Sensitivity analyses considered other reconstruction-error subsets.

These full-data training settings describe the original study methodology. The separately supplied pretrained `.keras` file is provided for the public inference demonstration; its file format and example execution should not be confused with the full-data PyTorch training procedure.

## Data availability

The discovery dataset is identified by GEO accession **GSE244263**. The example data and processed DEPA inputs supplied here support execution of the included analysis code. They do not replace the complete discovery data or constitute an independent validation cohort.

## Code and archive

- **GitHub:** https://github.com/khongorzuldash/ALS-DAE-Pretrained-Inference
- **Zenodo:** https://doi.org/10.5281/zenodo.23233282

For the study design, full analytical methods, and biological interpretation, please refer to the associated manuscript and supplementary information.
