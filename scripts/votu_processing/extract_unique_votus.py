#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_fasta(path):
    seqs = {}

    with open(path) as f:
        ident, seq_lines = None, []

        for line in f:
            if line.startswith(">"):
                if ident:
                    seqs[ident] = "".join(seq_lines)

                ident = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())

        if ident:
            seqs[ident] = "".join(seq_lines)

    return seqs


def get_representatives(
    tsv_file,
    fasta_file,
    ids_out,
    fasta_out,
    method="longest"
):
    seqs = parse_fasta(fasta_file)
    tsv = pd.read_csv(tsv_file, sep="\t", header=None)

    reps = []

    for _, row in tsv.iterrows():
        contigs = str(row[0]).split(",")

        if method == "longest":
            rep = max(
                contigs,
                key=lambda c: len(seqs.get(c, ""))
            )
        else:
            rep = contigs[0]

        reps.append(rep)

    with open(ids_out, "w") as f:
        for rep in reps:
            f.write(rep + "\n")

    with open(fasta_out, "w") as f:
        for rep in reps:
            if rep in seqs:
                f.write(f">{rep}\n{seqs[rep]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract representative vOTU sequences from clustered contigs."
    )

    parser.add_argument(
        "--clusters",
        required=True,
        help="Input TSV containing clustered contig IDs"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA containing contig sequences"
    )

    parser.add_argument(
        "--ids-out",
        required=True,
        help="Output file containing representative IDs"
    )

    parser.add_argument(
        "--fasta-out",
        required=True,
        help="Output FASTA containing representative sequences"
    )

    parser.add_argument(
        "--method",
        choices=["longest", "first"],
        default="longest",
        help="Representative selection method (default: longest)"
    )

    args = parser.parse_args()

    get_representatives(
        args.clusters,
        args.fasta,
        args.ids_out,
        args.fasta_out,
        args.method
    )


if __name__ == "__main__":
    main()
