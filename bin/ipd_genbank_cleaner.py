#!/usr/bin/env python3

import os
import sys
from datetime import datetime
from Bio import SeqIO

file_dir = os.getcwd()
animal = sys.argv[1]
locus = sys.argv[2]
input_file = sys.argv[3]

with open(input_file, "r") as input_handle, open(os.path.join(file_dir, "ipd-" +locus + "-" + animal + "-" + datetime.today().strftime('%Y-%m-%d') + "_cleaned.gbk"), "w") as output_handle:

	for seq_record in SeqIO.parse(input_handle, "genbank"):
		seq_record.seq = seq_record.seq.rstrip("X")
		if len(seq_record.seq) > 100:
			SeqIO.write(seq_record, output_handle, "genbank")
