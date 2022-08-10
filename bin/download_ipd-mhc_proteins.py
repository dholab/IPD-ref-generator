#!/usr/bin/env python3

import sys
from Bio import SeqIO
import requests
from datetime import datetime

ipd_number = int(sys.argv[1])

# create output genbank file for entire database
with open("ipd-mhc-nhp-prot-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".fasta", "a") as all_nhp:

  # create output genbank file for rhesus only
  with open("ipd-mhc-mamu-prot-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".fasta", "a") as mamu:

    # create output genbank file for cyno only
    with open("ipd-mhc-mafa-prot-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".fasta", "a") as mafa:

      # create output genbank file for mane only
      with open("ipd-mhc-mane-prot-" + datetime.today().strftime('%Y-%m-%d') + "_" + str(ipd_number) + ".fasta", "a") as mane:
        
        # add leading zeros and nhp prefix
        # this is the EBI dbfetch format
        nhp_id = 'NHP' + str(ipd_number).zfill(5)
        
       # get record
        u = requests.get("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdmhcpro;id=" + nhp_id + ";style=raw")
        
        ipd_fasta = u.text # read content
        
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

