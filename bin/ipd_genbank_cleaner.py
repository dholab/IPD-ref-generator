import os
# import glob
from Bio import SeqIO

# for file in glob.glob(os.getcwd() + "*.gbk"):

with open("ipd-mhc-mafa-2022-07-11.gbk", "rU") as input_handle, open("ipd-mhc-mafa-2022-07-11_cleaned.gbk", "w") as output_handle:

	for seq_record in SeqIO.parse(input_handle, "genbank"):
		seq_record.seq = seq_record.seq.rstrip("X")
		if len(seq_record.seq) > 100:
			SeqIO.write(seq_record, output_handle, "genbank")

with open("ipd-mhc-mane-2022-07-11.gbk", "rU") as input_handle, open("ipd-mhc-mane-2022-07-11_cleaned.gbk", "w") as output_handle:

	for seq_record in SeqIO.parse(input_handle, "genbank"):
		seq_record.seq = seq_record.seq.rstrip("X")
		if len(seq_record.seq) > 100:
			SeqIO.write(seq_record, output_handle, "genbank")

with open("ipd-mhc-mamu-2022-07-11.gbk", "rU") as input_handle, open("ipd-mhc-mamu-2022-07-11_cleaned.gbk", "w") as output_handle:

	for seq_record in SeqIO.parse(input_handle, "genbank"):
		seq_record.seq = seq_record.seq.rstrip("X")
		if len(seq_record.seq) > 100:
			SeqIO.write(seq_record, output_handle, "genbank")

with open("ipd-mhc-nhp-2022-07-11.gbk", "rU") as input_handle, open("ipd-mhc-nhp-2022-07-11_cleaned.gbk", "w") as output_handle:

	for seq_record in SeqIO.parse(input_handle, "genbank"):
		seq_record.seq = seq_record.seq.rstrip("X")
		if len(seq_record.seq) > 100:
			SeqIO.write(seq_record, output_handle, "genbank")
