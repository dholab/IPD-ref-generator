#!/usr/bin/env python3

"""
This script removes some of the peculiar IPD formatting conventions from an
input Genbank file so that it can be processed more fearlessly downstream.
"""

import argparse
from contextlib import ExitStack
from datetime import date
from Bio import SeqIO

def parse_command_line_args() -> tuple[str, str, str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--animal", "-a",
                        required=True,
                        type=str,
                        help="The species in question, e.g. mamu or mafa.")
    parser.add_argument("--gene", "-g",
                        required=True,
                        type=str,
                        help="The immunological gene of interest, e.g. MHC or KIR.")
    parser.add_argument("--file", "-f",
                        required=True,
                        type=str,
                        help="The input Genbank file to be cleaned")
    args = parser.parse_args()
    return args.animal, args.gene, args.file

def clean_ipd_genbank_file(animal: str, gene: str, file_path: str):
    """
    Remove X characters from IPD and output only sequences that
    are at least 100 nucleotides long.
    """

    # use an exit stack to concisely open input and output handles without too many indents
    with ExitStack() as stack:
        input_handle = stack.enter_context(open(file_path, "r", encoding="utf-8"))
        output_handle = stack.enter_context(open(f"ipd-{gene}-{animal}-{date.today()}_cleaned.gbk",
                                                 "a", encoding="utf-8"))

        # iterate through each record in the genbank file
        for seq_record in SeqIO.parse(input_handle, "genbank"):

            # strip X characters from IPD
            seq_record.seq = seq_record.seq.rstrip("X")

            # use an assertion to make sure the script isn't accidentally being used on
            # amino acid sequences
            allowed_characters = set("ATGCN-")
            assert all(char in allowed_characters for char in seq_record.seq), \
                "The sequence contains invalid characters. Only A, T, G, C, N, and '-' are allowed."

            # write to the output handle if the sequence is at least 100 bases long
            if len(seq_record.seq) >= 100:
                SeqIO.write(seq_record, output_handle, "genbank")

def main():
    """Daisy-chain the above functions when run together as a script"""

    # retrieve command line arguments
    animal, gene, file = parse_command_line_args()

    # supply arguments to genbank file cleaner function
    clean_ipd_genbank_file(animal, gene, file)

if __name__ == "__main__":
    main()
