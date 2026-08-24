#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO


def change_fasta_headers(
    input_file,
    output_file,
    pattern,
    replacement,
    remove_pattern
):
    """
    Modify FASTA headers by replacing an initial pattern
    and removing an unwanted substring.
    """

    records = []

    for record in SeqIO.parse(input_file, "fasta"):
        header = record.id

        if re.match(pattern, header):
            header = re.sub(pattern, replacement, header)

        header = re.sub(remove_pattern, "", header)

        record.id = header
        record.description = header
        records.append(record)

    SeqIO.write(records, output_file, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Modify FASTA sequence headers using regular expressions."
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
        "--pattern",
        required=True,
        help="Regular expression to match in the FASTA header"
    )

    parser.add_argument(
        "--replacement",
        default="",
        help="Replacement string for the matched pattern"
    )

    parser.add_argument(
        "--remove-pattern",
        default="",
        help="Regular expression describing an additional substring to remove"
    )

    args = parser.parse_args()

    change_fasta_headers(
        args.input,
        args.output,
        args.pattern,
        args.replacement,
        args.remove_pattern
    )


if __name__ == "__main__":
    main()
