#!/usr/bin/env python3
"""
This script assumes that there are many IPD allele EMBL files in the current working
directory, where each file contains a single record.
"""

import os
import argparse
from contextlib import ExitStack
from datetime import date

def parse_command_line_args() -> str:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", "-g",
                        required=True,
                        type=str,
                        help="The immunological gene of interest, e.g. MHC or KIR.")
    args = parser.parse_args()
    return args.gene

def collate_proteins(gene: str, current_date: str):
    """
    Work through the provided NHP protein records in EMBL format,
    sort them by animal, and then write them out to an all-NHP file 
    plus one file for each of the three major macaque species used 
    in research colonies.
    """

    # define gene string for file name
    gene_str = gene.lower()

    # define stack of file objects for each animal
    with ExitStack() as stack:
        all_nhp = stack.enter_context(open(f"ipd-{gene_str}-prot-nhp-{current_date}.fasta",
                                           "a", encoding="utf-8"))
        mamu = stack.enter_context(open(f"ipd-{gene_str}-prot-mamu-{current_date}.fasta",
                                        "a", encoding="utf-8"))
        mafa = stack.enter_context(open(f"ipd-{gene_str}-prot-mafa-{current_date}.fasta",
                                        "a", encoding="utf-8"))
        mane = stack.enter_context(open(f"ipd-{gene_str}-prot-mane-{current_date}.fasta",
                                        "a", encoding="utf-8"))

        # build a list of files to collate
        files = [f for f in os.listdir(".") if f.endswith('.fasta')]

        # iterate over all NHP files
        for file in files:
            with open(file, "r", encoding="utf-8") as input_handle:

                # read the file
                ipd_fasta = input_handle.read()

                # write to all NHP fasta
                all_nhp.write(ipd_fasta)

                # if rhesus sequence
                if 'Mamu' in ipd_fasta:
                    mamu.write(ipd_fasta)

                # if cyno sequence
                if 'Mafa' in ipd_fasta:
                    mafa.write(ipd_fasta)

                # if mane sequence
                if 'Mane' in ipd_fasta:
                    mane.write(ipd_fasta)

def main():
    """Gather information from the command line and run
    the provided functions"""

    # define the expected number of IPD alleles
    gene = parse_command_line_args()

    # store a string of the current date for file naming
    todays_date = str(date.today())

    # cross reference with previous download to bring in any
    # new alleles from the latest release
    collate_proteins(gene, todays_date)

if __name__ == "__main__":
    main()
