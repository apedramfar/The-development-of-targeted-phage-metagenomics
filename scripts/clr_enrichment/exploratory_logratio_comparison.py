import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
def categorize_from_m1(m1, thr):
    return np.where(m1 > thr, 'Enriched',
                    np.where(m1 < -thr, 'Depleted', 'Similar abundance'))
 
def export_summary(df, out_dir, label):
    summary = df['category'].value_counts().rename_axis('category').reset_index(name='count')
    total = len(df)
    summary['percent'] = (summary['count'] / total * 100).round(1)
    summary = pd.concat([summary, pd.DataFrame([{'category': 'All', 'count': total, 'percent': 100.0}])], ignore_index=True)
    summary.to_csv(os.path.join(out_dir, f'category_summary_{label}.csv'), index=False)
 
    enriched = df[df['category'] == 'Enriched'].sort_values('M1', ascending=False)
    enriched.to_csv(os.path.join(out_dir, f'all_enriched_{label}.csv'), index=False)
 
    df.to_csv(os.path.join(out_dir, f'all_labeled_{label}.csv'), index=False)
 
def scatter(df, xcol, ycol, xlabel, ylabel, title, out_path):
    plt.figure(figsize=(7,7))
    for cat, color in [('Enriched','red'),('Depleted','green'),('Similar abundance','blue')]:
        sub = df[df['category'] == cat]
        plt.scatter(sub[xcol], sub[ycol], c=color, label=f"{cat} (n={len(sub)})",
                    edgecolor='black', linewidth=0.3, alpha=0.8)
    mn = float(np.nanmin([df[xcol].min(), df[ycol].min()]))
    mx = float(np.nanmax([df[xcol].max(), df[ycol].max()]))
    plt.plot([mn,mx],[mn,mx],'--',color='gray')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Input CSV")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--pseudo", type=float, default=1e-3)
    ap.add_argument("--threshold", type=float, default=5.0)
    ap.add_argument("--col-a", default="TVB_Mean", help="Denominator column (default: TVB_Mean)")
    ap.add_argument("--col-b", default="MEPB_Mean", help="Numerator column (default: MEPB_Mean)")
    args = ap.parse_args()
 
    os.makedirs(args.out, exist_ok=True)
 
    df = pd.read_csv(args.csv)
    df.columns = df.columns.str.strip()
    if args.col_a not in df.columns or args.col_b not in df.columns:
        raise KeyError(f"Required columns not found. Available: {df.columns.tolist()}")
 
    df[args.col_a] = pd.to_numeric(df[args.col_a], errors='coerce')
    df[args.col_b] = pd.to_numeric(df[args.col_b], errors='coerce')
    df = df.dropna(subset=[args.col_a, args.col_b]).reset_index(drop=True)
 
    # ── Approach A: drop zeros
    df_A = df[(df[args.col_a] != 0) & (df[args.col_b] != 0)].copy()
    df_dropped = df[(df[args.col_a] == 0) | (df[args.col_b] == 0)].copy()
    df_dropped.to_csv(os.path.join(args.out, 'dropped_zero_rows.csv'), index=False)
 
    df_A['M1'] = np.log2(df_A[args.col_b] / df_A[args.col_a])
    df_A[f'log2_{args.col_a}'] = np.log2(df_A[args.col_a])
    df_A[f'log2_{args.col_b}'] = np.log2(df_A[args.col_b])
    df_A['category'] = categorize_from_m1(df_A['M1'], args.threshold)
    df_A.to_csv(os.path.join(args.out, 'drop_zeros_keep_negatives.csv'), index=False)
    scatter(
        df_A,
        xcol=f'log2_{args.col_a}', ycol=f'log2_{args.col_b}',
        xlabel=f'log₂({args.col_a})', ylabel=f'log₂({args.col_b})',
        title='Approach A: Zeros dropped',
        out_path=os.path.join(args.out, 'scatter_drop_zeros.png')
    )
 
    # ── Approach B: zeros → pseudo
    df_B = df.copy()
    df_B[f'{args.col_a}_adj']  = np.where(df_B[args.col_a]==0,  args.pseudo, df_B[args.col_a])
    df_B[f'{args.col_b}_adj']  = np.where(df_B[args.col_b]==0,  args.pseudo, df_B[args.col_b])
    df_B['M1'] = np.log2(df_B[f'{args.col_b}_adj'] / df_B[f'{args.col_a}_adj'])
    df_B['log2_A'] = np.log2(df_B[f'{args.col_a}_adj'])
    df_B['log2_B'] = np.log2(df_B[f'{args.col_b}_adj'])
    df_B['category'] = categorize_from_m1(df_B['M1'], args.threshold)
    df_B.to_csv(os.path.join(args.out, 'zero2pseudo_all.csv'), index=False)
    scatter(
        df_B,
        xcol='log2_A', ycol='log2_B',
        xlabel=f'log₂({args.col_a}) [zeros→pseudo]', ylabel=f'log₂({args.col_b}) [zeros→pseudo]',
        title='Approach B: Zeros → pseudo',
        out_path=os.path.join(args.out, 'scatter_zero2pseudo.png')
    )
 
    # ── Approach C: shift all by pseudo
    df_C = df.copy()
    df_C[f'{args.col_a}_shift']  = df_C[args.col_a] + args.pseudo
    df_C[f'{args.col_b}_shift']  = df_C[args.col_b] + args.pseudo
    df_C['M1'] = np.log2(df_C[f'{args.col_b}_shift'] / df_C[f'{args.col_a}_shift'])
    df_C['log2_A'] = np.log2(df_C[f'{args.col_a}_shift'])
    df_C['log2_B'] = np.log2(df_C[f'{args.col_b}_shift'])
    df_C['category'] = categorize_from_m1(df_C['M1'], args.threshold)
    df_C.to_csv(os.path.join(args.out, 'shift_all_by_pseudo.csv'), index=False)
    scatter(
        df_C,
        xcol='log2_A', ycol='log2_B',
        xlabel=f'log₂({args.col_a}_shift)', ylabel=f'log₂({args.col_b}_shift)',
        title='Approach C: All + pseudo',
        out_path=os.path.join(args.out, 'scatter_shift_all_by_pseudo.png')
    )
 
    # ── SD of M1
    sd_A = float(df_A['M1'].std(ddof=1))
    sd_B = float(df_B['M1'].std(ddof=1))
    sd_C = float(df_C['M1'].std(ddof=1))
    print(f"SD M1 (A): {sd_A:.3f}, (B): {sd_B:.3f}, (C): {sd_C:.3f}")
 
    # ── Summaries
    export_summary(df_A, args.out, 'drop')
    export_summary(df_B, args.out, 'zero2pseudo')
    export_summary(df_C, args.out, 'shift_all')
    print("✔ Done.")
 
if __name__ == "__main__":
    main()
