#!/usr/bin/env python3

import sys
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_genbank = sys.argv[1]

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

def separateExon2(in_gbk, exon2_fasta):

  # determine the coordinates for exon 2 in 
  # each allele, trim to those coordinates,
  # and save the sequence to a new file.

  with open(in_gbk) as input_handle, open(exon2_fasta, "w") as output_handle:
    
    # iterate over sequences
    sequences = SeqIO.parse(input_handle, "genbank")

    for record in sequences:
      
        # Find exon 2 feature
        for feature in record.features:
          
          if feature.type == "exon" and 'number' in feature.qualifiers and feature.qualifiers['number'] == ['2']:
          
            # Extract exon 2 sequence
            start = feature.location.start
            end = feature.location.end
            exon2_seq = record.seq[start:end]
          
            # Add exon 2 sequence to new record
            exon2_record = record[:]
            exon2_record.seq = exon2_seq
            exon2_record.id = exon2_record.name + "_exon2"
            # exon2_record.id = removeSpecialCharacters(exon2_record.id)
          
            # Write exon 2 record to output file
            SeqIO.write(exon2_record, output_handle, "fasta")
          
            # Stop searching for exon 2 features in this record
            break
          
    return output_handle
  
  
def deduplicate_fasta(in_fasta, out_fasta):
  '''manually deduplicate FASTA records with identical sequences
  create new FASTA file named with all matching sequences'''

  # create dictionary to store names and sequences

  deduplicated = dict()
  unique_seqs = set()

  ipd = SeqIO.parse(in_fasta, "fasta")

  # create a set of unique sequences
  for record in ipd:
    unique_seqs.add(str(record.seq))

  # create another biopython object
  ipd = SeqIO.parse(in_fasta, "fasta")

  # iterate over object.
  # add list elements when sequences match
  for record in ipd:
      for i in unique_seqs:
        if record.seq == i:
          # if key does not exist, add it
          if i not in deduplicated:
            deduplicated[i] = [record.id]
          else:
            deduplicated[i].append(record.id)

  # create deduplicated FASTA file from dictionary
  with open(out_fasta, "w") as output_handle:
    for item in deduplicated.items():
      # concatenate identical sequence names
      seq_name = ','.join(item[1])
      # get sequence
      deduplicated_sequence = item[0]
      # create biopython object
      record = SeqRecord(
      Seq(deduplicated_sequence),
      id=seq_name,
      description=''
      )
      SeqIO.write(record, output_handle, "fasta")

    return out_fasta


def splitByAlleleClass(in_fasta, classI_fasta, classII_fasta):
    
    '''
    Here we split the FASTA of exon-2-only sequences into separate
    databases for class I and class II MHC alleles. We do this 
    because these alleles sometimes have to be treated differently.
    '''
    
    # list of MHC class II loci to look for in the sequence deflines
    classII_alleles = ["DRA", "DRB", "DPA", "DPB", "DQA", "DQB"]
    
    with open(classII_fasta, "w") as classII, open(classI_fasta, "w") as classI:
        
        # for each record
        for record in SeqIO.parse(in_fasta, "fasta"):
            
            found = False
            
            # determine if the defline contains one of the class II alleles
            for locus in classII_alleles:
                
                if locus in record.id:
                    
                    found = True
                    break
        
            # If a sequence is class II, save it in the class II fasta.
            # Otherwise, save it in a class I fasta.
            if found:
                SeqIO.write(record, classII, "fasta")
            else:
                SeqIO.write(record, classI, "fasta")
    
        return classI, classII


def createExon2Fasta(in_gbk):
  '''given an IPD genbank file chain together above functions 
  to create FASTA file trimmed to exon 2 only, both for class
  I and class II alleles.
  '''

  # create output file names
  IPD_BASENAME = os.path.basename(in_gbk[:-4])
  IPD_EXON2_FASTA = IPD_BASENAME + '_exon2.fasta'
  IPD_DEDUP_FASTA = IPD_BASENAME + '_exon2_deduplicated.fasta'
  IPD_CLASS_I_FASTA = IPD_BASENAME + '_exon2_deduplicated_classI.fasta'
  IPD_CLASS_II_FASTA = IPD_BASENAME + '_exon2_deduplicated_classII.fasta'
  
  # make sure file names are available to subprocess environment
  os.environ['IPD_BASENAME'] = IPD_BASENAME
  os.environ['IPD_EXON2_FASTA'] = IPD_EXON2_FASTA
  os.environ['IPD_DEDUP_FASTA'] = IPD_DEDUP_FASTA
  os.environ['IPD_CLASS_I_FASTA'] = IPD_CLASS_I_FASTA
  os.environ['IPD_CLASS_II_FASTA'] = IPD_CLASS_II_FASTA

  separateExon2(in_gbk, IPD_EXON2_FASTA)

  # deduplicate cDNA
  deduplicate_fasta(IPD_EXON2_FASTA, IPD_DEDUP_FASTA)
  
  # separate out Class I and class II
  splitByAlleleClass(IPD_DEDUP_FASTA, IPD_CLASS_I_FASTA, IPD_CLASS_II_FASTA)

# Run all the steps
createExon2Fasta(input_genbank)
