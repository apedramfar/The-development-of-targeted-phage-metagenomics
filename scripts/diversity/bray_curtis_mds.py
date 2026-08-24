#!/usr/bin/env python3
import os
import argparse
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from itertools import cycle

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Input CSV")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--index-col", default="vOTUs", help="Column with sample IDs (default: vOTUs)")
    ap.add_argument("--group-regex", default=r"^([A-Za-z]+)", help="Regex to extract group from sample ID")
    ap.add_argument("--metric", default="braycurtis", help="Distance metric (default: braycurtis)")
    ap.add_argument("--n-components", type=int, default=2, help="MDS dimensions (default: 2)")
    ap.add_argument("--random-state", type=int, default=42, help="MDS random state")
    ap.add_argument("--palette", default="", help='JSON dict mapping group->color (e.g. \'{"TVB":"#E41A1C"}\')')
    ap.add_argument("--shapes", default="", help='JSON dict mapping group->marker (e.g. \'{"TVB":"v"}\')')
    ap.add_argument("--png-name", default="mds_plot.png")
    return ap.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.csv)
    if args.index_col not in df.columns:
        raise KeyError(f"Missing index column: {args.index_col}")
    df = df.copy()
    df.index = df[args.index_col].astype(str)
    X = df.drop(columns=[args.index_col])

    # keep numeric only
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.fillna(0.0)

    # distance + MDS
    dist_vec = pdist(X.values, metric=args.metric)
    dist_mat = squareform(dist_vec)
    dist_df = pd.DataFrame(dist_mat, index=X.index, columns=X.index)
    dist_df.to_csv(os.path.join(args.out, f"{args.metric}_matrix.csv"))

    mds = MDS(n_components=args.n_components, dissimilarity="precomputed", random_state=args.random_state)
    coords = mds.fit_transform(dist_mat)
    mds_df = pd.DataFrame(coords, index=X.index, columns=[f"MDS{i+1}" for i in range(args.n_components)])

    # groups
    groups = df[args.index_col].astype(str).str.extract(args.group_regex, expand=True)[0]
    groups.index = X.index
    mds_df["Group"] = groups

    # palettes/shapes (fallbacks)
    default_colors = cycle(["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"])
    default_markers = cycle(["o","s","^","v","D","P","X","<",">","h"])
    user_palette = json.loads(args.palette) if args.palette else {}
    user_shapes  = json.loads(args.shapes) if args.shapes else {}

    color_map = {}
    shape_map = {}
    for g in mds_df["Group"].dropna().unique():
        color_map[g] = user_palette.get(g, next(default_colors))
        shape_map[g] = user_shapes.get(g, next(default_markers))

    # plot
    plt.figure(figsize=(10, 8))
    for g, sub in mds_df.groupby("Group"):
        plt.scatter(
            sub["MDS1"], sub["MDS2"],
            label=str(g),
            c=color_map.get(g, "#333333"),
            marker=shape_map.get(g, "o"),
            s=90,
            edgecolor="black",
            linewidth=0.6,
            alpha=0.9
        )
    plt.xlabel("MDS1"); plt.ylabel("MDS2")
    plt.title(f"MDS ({args.metric})")
    plt.legend(title="Group", loc="best", frameon=True)
    plt.tight_layout()
    out_png = os.path.join(args.out, args.png_name)
    plt.savefig(out_png, dpi=300)
    plt.close()

    # save coords
    mds_df.to_csv(os.path.join(args.out, "mds_coordinates.csv"))
    print("Done.")
    print("Matrix:", os.path.join(args.out, f"{args.metric}_matrix.csv"))
    print("Coords:", os.path.join(args.out, "mds_coordinates.csv"))
    print("Plot:", out_png)

if __name__ == "__main__":
    main()
