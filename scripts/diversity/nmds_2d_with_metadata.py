import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

def nmds_plot(df, index_col=None, metric="braycurtis", random_state=42, annotate=True):
    """
    Minimal, general NMDS plotter.
    - df: features x samples table (rows=features, cols=samples)
    - index_col: set as index if your sample/feature ID column exists (or None)
    - metric: dissimilarity metric (default: braycurtis)
    - annotate: write sample names next to points
    Returns (coords, fig, ax)
    """
    if index_col:
        df = df.set_index(index_col)
    X = df.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")

    dist = pdist(X.T.values, metric=metric)
    dist_mat = squareform(dist)

    mds = MDS(n_components=2, dissimilarity="precomputed", metric=False, random_state=random_state)
    coords = mds.fit_transform(dist_mat)

    fig, ax = plt.subplots(figsize=(7,6))
    for i, name in enumerate(X.columns):
        ax.scatter(coords[i,0], coords[i,1], s=80, edgecolor="black")
        if annotate:
            ax.text(coords[i,0], coords[i,1], str(name), fontsize=8, ha="left", va="bottom")
    ax.set_xlabel("NMDS1"); ax.set_ylabel("NMDS2"); ax.set_title(f"NMDS ({metric})")
    plt.tight_layout()
    return coords, fig, ax

