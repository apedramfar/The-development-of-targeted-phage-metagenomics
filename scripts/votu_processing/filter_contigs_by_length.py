#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def filter_sequences(input_fasta, output_fasta, min_length=10000):
    records = (
        rec for rec in SeqIO.parse(input_fasta, "fasta")
        if len(rec.seq) >= min_length
    )

    count = SeqIO.write(records, output_fasta, "fasta")

    print(f"{count} sequences >= {min_length} bp written to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by minimum length."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )

    parser.add_argument(
        "-m", "--min-length",
        type=int,
        default=10000,
        help="Minimum sequence length in bp (default: 10000)"
    )

    args = parser.parse_args()

    filter_sequences(
        args.input,
        args.output,
        args.min_length
    )


if __name__ == "__main__":
    main()
