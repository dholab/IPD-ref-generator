#!/usr/bin/env python3

import sys
from contextlib import ExitStack
import csv
from datetime import date
import requests
from Bio import SeqIO

def main():
    """
    Barebones switch to better python practices, but more to come soon.
    """

    samplesheet_path = str(sys.argv[1])
    current_date = date.today()

    # define stack of file objects for each animal
    with ExitStack() as stack:
        csvfile = stack.enter_context(open(samplesheet_path, 
                                           "r", encoding="utf-8"))
        all_nhp = stack.enter_context(open(f"ipd-mhc-nhp-{current_date}_added.gbk",
                                           "a", encoding="utf-8"))
        mamu = stack.enter_context(open(f"ipd-mhc-mamu-{current_date}_added.gbk",
                                        "a", encoding="utf-8"))
        mafa = stack.enter_context(open(f"ipd-mhc-mafa-{current_date}_added.gbk",
                                        "a", encoding="utf-8"))
        mane = stack.enter_context(open(f"ipd-mhc-mane-{current_date}_added.gbk",
                                        "a", encoding="utf-8"))

        samplesheet = csv.reader(csvfile)

        # loop through samplesheet CSV
        for accession,formal_name,informal_name,ipd_id,animal_id in samplesheet:

            # get record
            db_fetch = requests.get(f"https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id={accession}&style=raw",
                                    timeout=30)

            # IPD MHC uses non-standard ID line
            # need to remove first two semicolons in ID line
            ipd_embl = db_fetch.text # read content

            # getting rid of unexpected characters
            ipd_embl = ipd_embl.replace('<', '')

            # handle missing records
            # these don't have identifiers and can be skipped
            if not ipd_embl.startswith('ID'):
                continue

            # need to remove first two semicolons in ID line
            ipd_line = ipd_embl.splitlines() # split by line
            id_line = ipd_line[0] # get ID line
            id_line_split = id_line.split(';') # split elements by semicolon

            # reconstruct ID line in EMBL format that can be parsed by biopython
            ipd_line[0] = f"{id_line_split[0]} {id_line_split[1]} {id_line_split[4]}; {id_line_split[3]}; {id_line_split[5]}; {id_line_split[6]}"

            # join lines to create embl file
            embl_file = ('\n').join(ipd_line)

            # create temporary file in correct EMBL format
            # i tried to use tempfile but couldn't get it to work
            with open("response.embl", "w", encoding="utf-8") as temp_handle:
                temp_handle.write(embl_file)

            # read EMBL file and export as Genbank
            for record in SeqIO.parse("response.embl", "embl"):

                record.description = record.name
                record.name = formal_name
                print(record.name + ' - ' + record.description)

                record.seq = record.seq.rstrip("X")

                if len(record.seq) > 100:
                    SeqIO.write(record, all_nhp, "genbank")

                # if rhesus sequence
                if formal_name.startswith('Mamu'):
                    SeqIO.write(record, mamu, "genbank")

                # if cyno sequence
                if formal_name.startswith('Mafa'):
                    SeqIO.write(record, mafa, "genbank")

                # if mane sequence
                if formal_name.startswith('Mane'):
                    SeqIO.write(record, mane, "genbank")

if __name__ == "__main__":
    main()
