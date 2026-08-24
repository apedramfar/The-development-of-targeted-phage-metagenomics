#!/usr/bin/env python3

import os
import argparse
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px


# ------------------------------------------------------------
# Defaults
# ------------------------------------------------------------

DEFAULT_GROUP_A = "TVB"
DEFAULT_GROUP_B = "MEPB"
DEFAULT_N_REPS = 3
DEFAULT_PSEUDO = 0.01
DEFAULT_BINS = 60

POINT_SIZE = 40
POINT_ALPHA = 0.9

COLOR_ENRICHED = "red"
COLOR_DEPLETED = "blue"
COLOR_SIMILAR = "gray"
HIGHLIGHT_COLOR = "green"


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

def find_replicate_cols(
    df: pd.DataFrame,
    base: str,
    n_reps: int = 3
) -> List[str]:
    """
    Find replicate columns beginning with, for example:
    TVB1_RPK, TVB2_RPK, TVB3_RPK.

    Additional text after the prefix is permitted.
    """
    cols = []

    for i in range(1, n_reps + 1):
        prefix = f"{base}{i}_RPK"
        matches = [c for c in df.columns if c.startswith(prefix)]

        if len(matches) == 0:
            raise KeyError(
                f"No column found beginning with '{prefix}'"
            )

        if len(matches) > 1:
            raise KeyError(
                f"Multiple columns found for '{prefix}': {matches}"
            )

        cols.append(matches[0])

    return cols


def compute_clr_per_sample(
    df: pd.DataFrame,
    sample_cols: List[str],
    pseudo: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate CLR values independently for each sample.

    For each sample:
        x -> x + pseudocount
        log2(x)
        CLR = log2(x) - mean(log2(x))

    The mean is calculated across all vOTUs within each sample.
    """

    dfo = df.copy()

    dfo[sample_cols] = dfo[sample_cols].apply(
        pd.to_numeric,
        errors="coerce"
    )

    invalid_mask = (
        dfo[sample_cols].isna().any(axis=1)
        | (dfo[sample_cols] < 0).any(axis=1)
    )

    dropped = dfo.loc[
        invalid_mask,
        ["vOTUs"] + sample_cols
    ].copy()

    dfo = dfo.loc[~invalid_mask].copy()

    x = dfo[sample_cols].to_numpy(dtype=float) + pseudo

    if np.any(x <= 0):
        raise ValueError(
            "Non-positive values remained after pseudocount addition."
        )

    log_x = np.log2(x)

    clr = log_x - log_x.mean(axis=0, keepdims=True)

    clr_cols = [f"{c}_CLR" for c in sample_cols]
    dfo[clr_cols] = clr

    return dfo, dropped


def compute_group_means(
    df: pd.DataFrame,
    group_a_cols: List[str],
    group_b_cols: List[str],
    group_a: str,
    group_b: str
) -> pd.DataFrame:
    """
    Calculate mean CLR for each group and the CLR enrichment difference.

    M = mean(group B CLR) - mean(group A CLR)
    A = average of the two group means
    """

    dfo = df.copy()

    group_a_clr = [f"{c}_CLR" for c in group_a_cols]
    group_b_clr = [f"{c}_CLR" for c in group_b_cols]

    mean_a_col = f"Mean_{group_a}_CLR"
    mean_b_col = f"Mean_{group_b}_CLR"

    dfo[mean_a_col] = dfo[group_a_clr].mean(axis=1)
    dfo[mean_b_col] = dfo[group_b_clr].mean(axis=1)

    dfo["M"] = dfo[mean_b_col] - dfo[mean_a_col]
    dfo["enrichment"] = dfo["M"]

    dfo["A"] = (
        dfo[mean_a_col] + dfo[mean_b_col]
    ) / 2

    return dfo


def assign_categories(
    df: pd.DataFrame,
    lower_threshold: float,
    upper_threshold: float
) -> pd.DataFrame:
    """
    Classify vOTUs according to CLR enrichment.

    Enriched:
        M > mean + 2 SD

    Depleted:
        M < mean - 2 SD

    Similar abundance:
        between the two thresholds
    """

    dfo = df.copy()

    dfo["category"] = np.select(
        [
            dfo["M"] > upper_threshold,
            dfo["M"] < lower_threshold
        ],
        [
            "Enriched",
            "Depleted"
        ],
        default="Similar abundance"
    )

    return dfo


def write_category_summary(
    df: pd.DataFrame,
    out_path: str
) -> None:

    order = [
        "Enriched",
        "Depleted",
        "Similar abundance"
    ]

    summary = (
        df["category"]
        .value_counts()
        .reindex(order, fill_value=0)
        .rename_axis("category")
        .reset_index(name="count")
    )

    summary["percent"] = (
        summary["count"] / len(df) * 100
    ).round(2)

    summary.to_csv(out_path, index=False)


# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def scatter_group_means(
    df: pd.DataFrame,
    group_a: str,
    group_b: str,
    out_path: str,
    highlight_votu: Optional[str] = None
) -> None:

    xcol = f"Mean_{group_a}_CLR"
    ycol = f"Mean_{group_b}_CLR"

    colors = {
        "Enriched": COLOR_ENRICHED,
        "Depleted": COLOR_DEPLETED,
        "Similar abundance": COLOR_SIMILAR
    }

    plt.figure(figsize=(7, 7))

    for category, color in colors.items():

        sub = df[df["category"] == category]

        if sub.empty:
            continue

        plt.scatter(
            sub[xcol],
            sub[ycol],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            c=color,
            edgecolor="black",
            linewidth=0.4,
            label=f"{category} (n={len(sub)})"
        )

    mn = min(df[xcol].min(), df[ycol].min())
    mx = max(df[xcol].max(), df[ycol].max())

    plt.plot(
        [mn, mx],
        [mn, mx],
        "--",
        color="black",
        linewidth=1
    )

    if highlight_votu:

        hi = df[
            df["vOTUs"].astype(str)
            == str(highlight_votu)
        ]

        if not hi.empty:

            plt.scatter(
                hi[xcol],
                hi[ycol],
                s=80,
                c=HIGHLIGHT_COLOR,
                edgecolor="black",
                linewidth=1,
                zorder=10
            )

            plt.annotate(
                highlight_votu,
                (
                    hi[xcol].iloc[0],
                    hi[ycol].iloc[0]
                ),
                xytext=(10, 10),
                textcoords="offset points"
            )

    plt.xlabel(xcol)
    plt.ylabel(ycol)

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        out_path,
        dpi=300
    )

    plt.close()


def plot_ma(
    df: pd.DataFrame,
    lower_threshold: float,
    upper_threshold: float,
    global_mean: float,
    out_path: str,
    highlight_votu: Optional[str] = None
) -> None:

    colors = {
        "Enriched": COLOR_ENRICHED,
        "Depleted": COLOR_DEPLETED,
        "Similar abundance": COLOR_SIMILAR
    }

    plt.figure(figsize=(8, 6))

    for category, color in colors.items():

        sub = df[df["category"] == category]

        if sub.empty:
            continue

        plt.scatter(
            sub["A"],
            sub["M"],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            c=color,
            edgecolor="black",
            linewidth=0.4,
            label=f"{category} (n={len(sub)})"
        )

    plt.axhline(
        global_mean,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axhline(
        upper_threshold,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axhline(
        lower_threshold,
        linestyle="--",
        color="black",
        linewidth=1
    )

    if highlight_votu:

        hi = df[
            df["vOTUs"].astype(str)
            == str(highlight_votu)
        ]

        if not hi.empty:

            plt.scatter(
                hi["A"],
                hi["M"],
                s=80,
                c=HIGHLIGHT_COLOR,
                edgecolor="black",
                linewidth=1,
                zorder=10
            )

            plt.annotate(
                highlight_votu,
                (
                    hi["A"].iloc[0],
                    hi["M"].iloc[0]
                ),
                xytext=(10, 10),
                textcoords="offset points"
            )

    plt.xlabel("A (mean CLR)")
    plt.ylabel("M (CLR enrichment)")
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        out_path,
        dpi=300
    )

    plt.close()


def write_histogram(
    df: pd.DataFrame,
    out_dir: str,
    bins: int,
    lower_threshold: float,
    upper_threshold: float
) -> None:

    color_map = {
        "Enriched": COLOR_ENRICHED,
        "Depleted": COLOR_DEPLETED,
        "Similar abundance": COLOR_SIMILAR
    }

    fig = px.histogram(
        df,
        x="M",
        color="category",
        nbins=bins,
        color_discrete_map=color_map
    )

    fig.update_traces(
        marker_line_color="black",
        marker_line_width=1
    )

    fig.add_vline(
        x=upper_threshold,
        line_dash="dash",
        annotation_text="Mean + 2SD"
    )

    fig.add_vline(
        x=lower_threshold,
        line_dash="dash",
        annotation_text="Mean - 2SD"
    )

    fig.update_layout(
        xaxis_title="CLR enrichment (M)",
        yaxis_title="Count",
        title="Distribution of CLR enrichment"
    )

    fig.write_html(
        os.path.join(
            out_dir,
            "clr_enrichment_histogram.html"
        )
    )

    try:
        fig.write_image(
            os.path.join(
                out_dir,
                "clr_enrichment_histogram.png"
            )
        )
    except Exception:
        pass


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "CLR-based enrichment analysis of "
            "minicell-enriched vOTUs."
        )
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Input CSV containing vOTU RPK values."
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--group-a",
        default=DEFAULT_GROUP_A,
        help="Reference group prefix (default: TVB)."
    )

    parser.add_argument(
        "--group-b",
        default=DEFAULT_GROUP_B,
        help="Enriched group prefix (default: MEPB)."
    )

    parser.add_argument(
        "--n-reps",
        type=int,
        default=DEFAULT_N_REPS,
        help="Number of replicates per group (default: 3)."
    )

    parser.add_argument(
        "--pseudo",
        type=float,
        default=DEFAULT_PSEUDO,
        help="Pseudocount added before CLR transformation (default: 0.01)."
    )

    parser.add_argument(
        "--bins",
        type=int,
        default=DEFAULT_BINS,
        help="Histogram bins (default: 60)."
    )

    parser.add_argument(
        "--highlight-votu",
        default=None,
        help="Optional vOTU ID to highlight in plots."
    )

    args = parser.parse_args()

    os.makedirs(
        args.out,
        exist_ok=True
    )

    df = pd.read_csv(args.csv)
    df.columns = df.columns.str.strip()

    if "vOTUs" not in df.columns:
        raise KeyError(
            "Input CSV must contain a 'vOTUs' column."
        )

    group_a_cols = find_replicate_cols(
        df,
        args.group_a,
        args.n_reps
    )

    group_b_cols = find_replicate_cols(
        df,
        args.group_b,
        args.n_reps
    )

    sample_cols = (
        group_a_cols
        + group_b_cols
    )

    print("Reference group columns:")
    print(group_a_cols)

    print("Enriched group columns:")
    print(group_b_cols)

    print(
        f"Pseudocount: {args.pseudo}"
    )

    clr_df, dropped = compute_clr_per_sample(
        df[["vOTUs"] + sample_cols],
        sample_cols,
        args.pseudo
    )

    if not dropped.empty:

        dropped.to_csv(
            os.path.join(
                args.out,
                "dropped_invalid_rows.csv"
            ),
            index=False
        )

    clr_df.to_csv(
        os.path.join(
            args.out,
            "clr_table.csv"
        ),
        index=False
    )

    result = compute_group_means(
        clr_df,
        group_a_cols,
        group_b_cols,
        args.group_a,
        args.group_b
    )

    global_mean = result["M"].mean()

    global_sd = result["M"].std(
        ddof=1
    )

    upper_threshold = (
        global_mean
        + 2 * global_sd
    )

    lower_threshold = (
        global_mean
        - 2 * global_sd
    )

    result = assign_categories(
        result,
        lower_threshold,
        upper_threshold
    )

    result.to_csv(
        os.path.join(
            args.out,
            "clr_enrichment_results.csv"
        ),
        index=False
    )

    write_category_summary(
        result,
        os.path.join(
            args.out,
            "category_summary.csv"
        )
    )

    result[
        result["category"]
        == "Enriched"
    ].sort_values(
        "M",
        ascending=False
    ).to_csv(
        os.path.join(
            args.out,
            "enriched_vOTUs.csv"
        ),
        index=False
    )

    result[
        result["category"]
        == "Depleted"
    ].sort_values(
        "M"
    ).to_csv(
        os.path.join(
            args.out,
            "depleted_vOTUs.csv"
        ),
        index=False
    )

    result[
        result["category"]
        == "Similar abundance"
    ].to_csv(
        os.path.join(
            args.out,
            "similar_abundance_vOTUs.csv"
        ),
        index=False
    )

    summary_stats = pd.DataFrame(
        [{
            "pseudo": args.pseudo,
            "log_base": "log2",
            "global_mean_M": global_mean,
            "global_sd_M": global_sd,
            "lower_mean_minus_2sd": lower_threshold,
            "upper_mean_plus_2sd": upper_threshold,
            "n_total": len(result),
            "n_enriched": (
                result["category"]
                == "Enriched"
            ).sum(),
            "n_depleted": (
                result["category"]
                == "Depleted"
            ).sum(),
            "n_similar": (
                result["category"]
                == "Similar abundance"
            ).sum()
        }]
    )

    summary_stats.to_csv(
        os.path.join(
            args.out,
            "summary_2SD.csv"
        ),
        index=False
    )

    scatter_group_means(
        result,
        args.group_a,
        args.group_b,
        os.path.join(
            args.out,
            "clr_group_mean_scatter.png"
        ),
        args.highlight_votu
    )

    plot_ma(
        result,
        lower_threshold,
        upper_threshold,
        global_mean,
        os.path.join(
            args.out,
            "clr_MA_plot.png"
        ),
        args.highlight_votu
    )

    write_histogram(
        result,
        args.out,
        args.bins,
        lower_threshold,
        upper_threshold
    )

    print()
    print("=== CLR enrichment summary ===")
    print(f"Total vOTUs: {len(result)}")

    print(
        "Enriched:",
        (
            result["category"]
            == "Enriched"
        ).sum()
    )

    print(
        "Depleted:",
        (
            result["category"]
            == "Depleted"
        ).sum()
    )

    print(
        "Similar abundance:",
        (
            result["category"]
            == "Similar abundance"
        ).sum()
    )

    print(
        f"Mean M: {global_mean:.4f}"
    )

    print(
        f"SD M: {global_sd:.4f}"
    )

    print(
        f"Lower threshold: {lower_threshold:.4f}"
    )

    print(
        f"Upper threshold: {upper_threshold:.4f}"
    )

    print()
    print(
        f"Outputs written to: {args.out}"
    )


if __name__ == "__main__":
    main()
