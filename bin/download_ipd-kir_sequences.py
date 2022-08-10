#!/usr/bin/env python3

import sys
from Bio import SeqIO
import requests
from datetime import datetime

ipd_number = int(sys.argv[1])

# create output genbank file for entire database
with open("ipd-kir-nhp-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".gbk", "a") as all_nhp:

  # create output genbank file for rhesus only
  with open("ipd-kir-mamu-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".gbk", "a") as mamu:

    # create output genbank file for cyno only
    with open("ipd-kir-mafa-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".gbk", "a") as mafa:

      # create output genbank file for mane only
      with open("ipd-kir-mane-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".gbk", "a") as mane:

        # add leading zeros and nhp prefix
        # this is the EBI dbfetch format
        nhp_id = 'NHP' + str(ipd_number).zfill(5)
        
        # get record
        u = requests.get("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdnhkir;id=" + nhp_id + ";style=raw")

        # IPD MHC uses non-standard ID line
        # need to remove first two semicolons in ID line
        ipd_embl = u.text # read content

        # handle missing records
        # these don't have identifiers and can be skipped
        if not ipd_embl.startswith('ID'):
          all_nhp.close()
          mamu.close()
          mafa.close()
          mane.close()
        else:
          ipd_line = ipd_embl.splitlines() # split by line
          id_line = ipd_line[0] # get ID line
          id_line_split = id_line.split(';') # split elements by semicolon
        
          # reconstruct ID line in EMBL format that can be parsed by biopython
          ipd_line[0] = id_line_split[0] + ' ' + id_line_split[1] + ' ' + id_line_split[2] + '; ' + id_line_split[3] + '; ' + id_line_split[4] + '; ' + id_line_split[5]
        
          # join lines to create embl file
          embl_file = ('\n').join(ipd_line)
        
          # create temporary file in correct EMBL format
          # i tried to use tempfile but couldn't get it to work
          with open("response.embl", "w") as f:
              f.write(embl_file)
        
          # read EMBL file and export as Genbank
          for record in SeqIO.parse("response.embl", "embl"):
              record.description = record.name
              
              record.name = record.annotations['keywords'][0]
              print(record.name + ' - ' + record.description)
              SeqIO.write(record, all_nhp, "genbank")
        
              # if rhesus sequence
              if record.name.startswith('Mamu'):
                SeqIO.write(record, mamu, "genbank")
        
              # if cyno sequence
              if record.name.startswith('Mafa'):
                SeqIO.write(record, mafa, "genbank")
        
              # if mane sequence
              if record.name.startswith('Mane'):
                SeqIO.write(record, mane, "genbank")
