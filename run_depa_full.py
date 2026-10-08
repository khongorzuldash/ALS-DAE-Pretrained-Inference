import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.manifold import MDS
from sklearn.cluster import KMeans

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "depa"
OUT = ROOT / "results" / "depa_full"

def dtw_distance(a, b, window=5):
    n, m = len(a), len(b)
    w = max(window, abs(n-m))
    dp = np.full((n+1, m+1), np.inf)
    dp[0, 0] = 0.0
    for i in range(1, n+1):
        for j in range(max(1, i-w), min(m, i+w)+1):
            cost = abs(a[i-1] - b[j-1])
            dp[i, j] = cost + min(dp[i-1, j], dp[i, j-1], dp[i-1, j-1])
    return dp[n, m]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    X = pd.read_csv(DATA / "DTW_input_P005_gene_by_sample_zscore.csv", index_col=0)
    genes = X.index.astype(str).tolist()
    values = X.to_numpy(dtype=float)
    n = len(genes)
    D = np.zeros((n, n), dtype=np.float32)
    print("Calculating DTW:", n, "genes", flush=True)
    for i in range(n):
        if i % 100 == 0: print("DTW progress:", i, "/", n, flush=True)
        for j in range(i+1, n):
            D[i, j] = dtw_distance(values[i], values[j], window=5)
            D[j, i] = D[i, j]
    pd.DataFrame(D, index=genes, columns=genes).to_csv(OUT / "DTW_distance_matrix.csv")
    print("Running MDS...", flush=True)
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, normalized_stress="auto")
    coords = mds.fit_transform(D)
    result = pd.DataFrame({"gene_name": genes, "MDS1": coords[:, 0], "MDS2": coords[:, 1]})
    result.to_csv(OUT / "MDS_coordinates.csv", index=False)
    print("Running K-means...", flush=True)
    km = KMeans(n_clusters=6, random_state=42, n_init=200)
    result["cluster"] = km.fit_predict(coords) + 1
    result.to_csv(OUT / "K6_gene_clusters.csv", index=False)
    print("DONE:", OUT, flush=True)

if __name__ == "__main__":
    main()
