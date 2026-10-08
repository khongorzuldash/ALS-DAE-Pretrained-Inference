import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "depa"

def dtw_distance(s1, s2, window=5):
    n, m = len(s1), len(s2)
    window = max(window, abs(n-m))
    dp = np.full((n+1, m+1), np.inf, dtype=float)
    dp[0, 0] = 0.0
    for i in range(1, n+1):
        for j in range(max(1, i-window), min(m, i+window)+1):
            cost = abs(s1[i-1] - s2[j-1])
            dp[i, j] = cost + min(dp[i-1, j], dp[i, j-1], dp[i-1, j-1])
    return dp[n, m]

def main():
    X = pd.read_csv(DATA / "DTW_input_P005_gene_by_sample_zscore.csv", index_col=0)
    D = pd.read_csv(DATA / "DTW_distance_matrix_P005.csv", index_col=0)
    mds = pd.read_csv(DATA / "DTW_MDS_coordinates_P005.csv")
    original = pd.read_csv(DATA / "DTW_MDS_kmeans_k6_gene_clusters_P005.csv")

    assert X.shape == (1272, 40), f"Unexpected input shape: {X.shape}"
    assert D.shape == (1272, 1272), f"Unexpected DTW shape: {D.shape}"

    genes = X.index[:10].tolist()
    errors = [
        abs(dtw_distance(X.loc[a].to_numpy(), X.loc[b].to_numpy()) - D.loc[a, b])
        for i, a in enumerate(genes) for b in genes[i+1:]
    ]
    dtw_pass = np.allclose(errors, 0, atol=1e-4)

    df = mds.merge(
        original[["gene_name", "DTW_kmeans_cluster"]],
        on="gene_name", validate="one_to_one"
    )
    km = KMeans(n_clusters=6, random_state=42, n_init=200)
    predicted = km.fit_predict(df[["MDS1", "MDS2"]].to_numpy())
    ari = adjusted_rand_score(df["DTW_kmeans_cluster"], predicted)

    print("DEPA GENES:", len(df))
    print("DTW PAIRS TESTED:", len(errors))
    print("DTW MAX ABS ERROR:", max(errors))
    print("DTW TEST:", "PASS" if dtw_pass else "FAIL")
    print("KMEANS ARI:", round(ari, 8))
    print("KMEANS TEST:", "PASS" if np.isclose(ari, 1.0) else "FAIL")

    if not dtw_pass or not np.isclose(ari, 1.0):
        raise SystemExit("DEPA reproducibility test failed")

if __name__ == "__main__":
    main()
