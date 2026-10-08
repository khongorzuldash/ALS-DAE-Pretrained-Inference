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
