#!/usr/bin/env python3
"""
This script assumes that there are many IPD allele EMBL files in the current working
directory, where each file contains a single record.
"""

import os
import argparse
from contextlib import ExitStack
from typing import IO, Optional
from datetime import date
from Bio import SeqIO

def parse_command_line_args() -> tuple[Optional[str], str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--previous_database", "-d",
                        required=False,
                        default=None,
                        type=str,
                        help="A previous database for the gene in question, if available.")
    parser.add_argument("--gene", "-g",
                        required=True,
                        type=str,
                        help="The immunological gene of interest, e.g. MHC or KIR.")
    args = parser.parse_args()
    return args.previous_database, args.gene

def parse_ipd_file(file_path: str):
    """
    IPD uses slightly unconventional specs for their EMBL files, so we need
    to oversee the parsing more carefully than using SeqIO.parse alone.
    """

    with open(file_path, "r", encoding="utf-8") as input_handle:

        # read the file
        ipd_embl = input_handle.read()

        # IPD MHC uses non-standard ID line
        # need to remove first two semicolons in ID line
        ipd_line = ipd_embl.splitlines() # split by line
        id_line = ipd_line[0] # get ID line
        id_line_split = id_line.split(';') # split elements by semicolon

        # reconstruct ID line in EMBL format that can be parsed by biopython
        ipd_line[0] = f"{id_line_split[0]} {id_line_split[1]} {id_line_split[2]}; {id_line_split[3]}; {id_line_split[4]}; {id_line_split[5]}"

        # join lines to create embl file
        corrected_embl = ('\n').join(ipd_line)

    return corrected_embl

def parse_and_write_alleles(embl_file: str, variables: dict[str, IO[str]]) -> list[str]:
    """
    Re-shuffle some of the information in the embl file and write out
    each record to a species-specific Genbank file. These Genbank files
    will constitute the new reference databases.
    """

    # make a dictionary the specify which files to write to based on the species
    prefix_to_file = {
        'Mamu': variables['mamu'],
        'Mafa': variables['mafa'],
        'Mane': variables['mane']
    }

    # make an empty list that will keep track of which alleles are
    # new or updated
    updated_accessions = []

    # read corrected EMBL file with SeqIO.parse and export as Genbank
    for record in SeqIO.parse(embl_file, "embl"):

        # shuffle around some of the information and strip
        # IPD's X characters
        record.description = record.name
        record.name = record.annotations['keywords'][0]
        record.seq = record.seq.rstrip("X")

        # ignore if the allele is truncated and therefore useless
        if len(record.seq) < 100:
            continue

        # Write to the 'all_nhp' file regardless
        SeqIO.write(record, variables['all_nhp'], "genbank")

        # determine the species and attempt to write out
        species = record.name[:4]
        file_to_write = prefix_to_file.get(species, None)
        if file_to_write:
            SeqIO.write(record, file_to_write, "genbank")

        # keep track of which alleles were added or updated in the latest release
        updated_accessions.append(record.name)

    return updated_accessions

def collate_alleles(old_database: Optional[str], gene: str, current_date: str) -> Optional[int]:
    """
    Work through the provided NHP allele records in EMBL format,
    sort them by animal, reformat them to Genbank format, and then
    write them out to an all-NHP file plus one file for each of the
    three major macaque species used in research colonies.
    """

    # define gene string for file name
    gene_str = gene.lower()

    # name a temporary embl file
    temp_name = "response.embl"

    # define stack of file objects for each animal
    with ExitStack() as stack:
        temp_embl = stack.enter_context(open(temp_name, "a", encoding="utf-8"))
        all_nhp = stack.enter_context(open(f"ipd-{gene_str}-nhp-{current_date}.gbk","a", encoding="utf-8"))
        mamu = stack.enter_context(open(f"ipd-{gene_str}-mamu-{current_date}.gbk", "a", encoding="utf-8"))
        mafa = stack.enter_context(open(f"ipd-{gene_str}-mafa-{current_date}.gbk", "a", encoding="utf-8"))
        mane = stack.enter_context(open(f"ipd-{gene_str}-mane-{current_date}.gbk", "a", encoding="utf-8"))

        # make an empty list that will keep track of which alleles are
        # new or updated
        updated_accessions = []

        # build a list of files to collate
        files = [f for f in os.listdir(".") if f.endswith('.embl')]

        # iterate over all NHP files
        for file in files:

            # prepare IPD formatting for conversion
            corrected_embl = parse_ipd_file(file)

            # create temporary file in correct EMBL format
            temp_embl.write(corrected_embl)

        # keep track of which alleles were added or updated in the latest release
        # as the alleles are parsed
        updated_accessions = parse_and_write_alleles(temp_name, locals())

        # reconcile with a previous database, if provided
        if not old_database:
            return len(updated_accessions)

        # make a dictionary the specify which files to write to
        # based on the species
        write_lookup = {
            'Mamu': mamu,
            'Mafa': mafa,
            'Mane': mane
        }

        # iterate through each record and write to the correct file
        for record in SeqIO.parse(old_database, "genbank"):

            # skip the record if it was recently updated
            if record.name in updated_accessions:
                continue

            # Write to the 'all_nhp' file regardless
            SeqIO.write(record, all_nhp, "genbank")

            # Determine which file to write to based on the species
            species = record.name[:4]
            file_to_write = write_lookup.get(species, None)
            if file_to_write:
                SeqIO.write(record, file_to_write, "genbank")

        return len(updated_accessions)

def main():
    """Gather information from the command line and run
    the provided functions"""

    # define the expected number of IPD alleles
    old_db, gene = parse_command_line_args()

    # store a string of the current date for file naming
    todays_date = str(date.today())

    # cross reference with previous download to bring in any
    # new alleles from the latest release
    updated_count = collate_alleles(old_db, gene, todays_date)

    # tell the user how many new alleles were downloaded
    print(f"{updated_count} updated alleles were downloaded and added to the reference database.")

if __name__ == "__main__":
    main()
