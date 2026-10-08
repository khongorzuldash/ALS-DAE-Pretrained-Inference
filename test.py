from pathlib import Path
import numpy as np
import os
os.environ['KERAS_BACKEND'] = 'tensorflow'
import tensorflow as tf

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / 'model/FULL_HVG3000_DAE_best.keras'
DATA = ROOT / 'data/sample_logNorm_100cells.npz'

X = np.load(DATA)['X'].astype(np.float32)
gene_min = np.load(ROOT / 'model/gene_min.npy')
gene_range = np.load(ROOT / 'model/gene_range.npy')

assert X.shape[1] == 3000
assert gene_min.shape == gene_range.shape == (1, 3000)

X_scaled = ((X - gene_min) / gene_range).astype(np.float32)
model = tf.keras.models.load_model(MODEL, compile=False)
reconstructed = model.predict(X_scaled, batch_size=32, verbose=0)

errors = np.mean((X_scaled - reconstructed) ** 2, axis=1)
n_select = int(np.ceil(len(errors) * 0.20))
selected = np.argsort(errors)[:n_select]

np.savetxt(ROOT / 'reconstruction_errors.csv', errors, delimiter=',', header='reconstruction_mse', comments='')
np.savetxt(ROOT / 'selected_cell_indices.csv', selected, fmt='%d', delimiter=',', header='cell_index_zero_based', comments='')

print('PASS: pretrained DAE loaded')
print('INPUT:', X.shape)
print('RECONSTRUCTION:', reconstructed.shape)
print('SELECTED CELLS:', len(selected), '/', len(errors))
print('MSE RANGE:', float(errors.min()), 'to', float(errors.max()))
print('SELECTED INDICES:', selected.tolist())
