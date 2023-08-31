#!/usr/bin/env python3

import sys
import os
import subprocess
from datetime import datetime
import logging
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def printStatus(status):
    '''print timestamped status update'''
    print('--[' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '] ' + status + '--')
    log.info(status)
    
def createOutputFolder(cwd):
    '''create timestamped output folder at specified location'''
    
    # fetch current time
    CURRENT_TIME = datetime.now().strftime("%Y%m%d%H%M%S")

    # path to output folder
    OUTPUT_FOLDER = cwd + '/' + CURRENT_TIME

    # create folder if it doesn't already exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # print output folder name
    printStatus('Output folder: ' + OUTPUT_FOLDER)
    
    return OUTPUT_FOLDER

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

def createFastaFromGenbank(in_gbk, out_fasta):
  '''create FASTA file from Genbank file'''

  printStatus('Create FASTA file from Genbank file')

  with open(in_gbk) as input_handle, open(out_fasta, "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "genbank")
    SeqIO.write(sequences, output_handle, "fasta") 

  return out_fasta

def vsearchMapToReference(query, ref, out_sam):
    '''map FASTA sequences to exemplar FASTA for exon 2 trimming'''
    # map IPD FASTA reads to exemplar sequences
    # SAM output CIGAR strings define where matches begin and end
    # which allows trimming IPD sequences to exon 2 length
    # write samheader so pysam parsing works correctly
    
    printStatus('Map to reference with vsearch usearch_global algorithm')
    
    os.environ['query'] = query
    os.environ['ref'] = ref
    os.environ['out_sam'] = out_sam
    
    subprocess.run('vsearch --usearch_global $query --db $ref --id 0.8 --samheader --samout $out_sam', shell = True)

    return out_sam

def trimQueryByMapping(in_sam):
    '''use CIGAR coordinates to trim SAM file to reference length
    return python list with trimmed sequences and sequence names
    '''

    printStatus('Deterine location of mapped regions of query')

    # import SAM file as pysam object
    samfile = pysam.AlignmentFile(in_sam, 'r', check_sq=False)
    
    os.environ['in_sam'] = in_sam

    # remove temporary sam file
    subprocess.run('rm $in_sam', shell = True)

    # create empty list of lists to store trimmed sequences
    trimmed_sequence = []

    # iterate over each read to determine trim locations
    # save trimmed sequence to list

    for read in samfile:
        if not read.is_unmapped:
            cigar = read.cigar

            # get number of bases to trim from left
            left_trim = 0 # default when there is no left trim
            
            # if the first entry in cigar string is insertion
            # that number of bases needs to be trimmed from sequence
            if cigar[0][0] == 1:
                left_trim = (cigar[0][1])
                # remove bases from left
                left_trimmed = read.seq[left_trim:]
            else:
                left_trimmed = read.seq

            # get number of bases to trim from right
            right_trim = 0 # default when there is no right trim

            if cigar[-1][0] == 1:
                right_trim = (cigar[-1][1])
                # remove bases from right
                right_trimmed = left_trimmed[:-right_trim]
            else:
                right_trimmed = left_trimmed

            # add to list
            trimmed_sequence.append([read.query_name, right_trimmed])

    return trimmed_sequence

def replaceGenbankSequence(in_gbk, trimmed_sequence_list, out_fasta):
    '''substitute trimmed sequences for original sequences in Genbank file
    expects trimmed_sequence_list object.
    outputs FASTA, since original genbank coordinates make no sense'''

    printStatus('Create trimmed FASTA')

    # load IPD Genbank file as Biopython object
    ipd = SeqIO.parse(in_gbk, "genbank")

    # replace original IPD sequence with trimmed sequence for miSeq genotyping
    # output as FASTA file
    # (Genbank annotation positions make no sense after trimming)
    # edit FASTA headers to be more informative
    with open(out_fasta, "w") as output_handle:
        for record in ipd:
            for i in trimmed_sequence_list:
                if i[0] == record.description:
                    record.seq = Seq(i[1])
                    record.id = removeSpecialCharacters(record.name)
                    SeqIO.write(record, output_handle, "fasta")

    return out_fasta

def deduplicate_fasta(in_fasta, out_fasta):
    '''manually deduplicate FASTA records with identical sequences
    create new FASTA file named with all matching sequences'''

    printStatus('Remove duplicates in FASTA sequence and rename to retain duplicate names')

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
            seq_name = '|'.join(item[1])
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

def separateGdnaFromCdna(in_gbk, gdna_fasta, cdna_fasta):
    '''
    parse genbank file into FASTA files with and without introns
    '''

    printStatus('Create gDNA and cDNA files from Genbank file')

    # test if a record has any intron annotations
    # if so, write to output gDNA FASTA file
    # if not, write to output cDNA FASTA file

    with open(in_gbk) as input_handle, open(gdna_fasta, "w") as gdna_handle, open(cdna_fasta, "w") as cdna_handle:
      
      # iterate over sequences
      sequences = SeqIO.parse(input_handle, "genbank")

      for record in sequences:
            # flag for presence of intron
            has_intron = 0
            
            for feature in record.features:
                if feature.type == 'intron':
                    has_intron = 1

            # if intron flag == 1, write to gDNA file
            # otherwise write to cDNA file
            if has_intron == 1:
                record.description = record.id
                record.id = removeSpecialCharacters(record.name)
                SeqIO.write(record, gdna_handle, "fasta") 
            else:
                SeqIO.write(record, cdna_handle, "fasta")

      return gdna_fasta, cdna_fasta

def main():
    '''given an IPD genbank file and an exemplar fasta
    chain together above functions to create FASTA file trimmed to exemplar sequences
    here exemplar should contain complete exon 2 sequences.
    then combine gDNA and cDNA sequences to single database
    '''

    animal = sys.argv[1]
    in_gbk = sys.argv[2]
    exemplar_fasta = sys.argv[3]

    # create output file names
    IPD_BASENAME = os.path.basename(in_gbk[:-4])
    IPD_TRIMMED_FASTA = IPD_BASENAME + '.exon2.trimmed.fasta'
    IPD_EXON2_FASTA = IPD_BASENAME + '.exon2.fasta'
    IPD_IMMUNOWES_FASTA = IPD_BASENAME + '.immunowes.fasta'
    IPD_CDNA_FASTA = IPD_BASENAME + '.cdna.fasta'
    IPD_GDNA_FASTA = IPD_BASENAME + '.gdna.fasta'
    
    # make sure file names are available to subprocess environment
    os.environ['IPD_BASENAME'] = IPD_BASENAME
    os.environ['IPD_TRIMMED_FASTA'] = IPD_TRIMMED_FASTA
    os.environ['IPD_EXON2_FASTA'] = IPD_EXON2_FASTA
    os.environ['IPD_IMMUNOWES_FASTA'] = IPD_IMMUNOWES_FASTA
    os.environ['IPD_CDNA_FASTA'] = IPD_CDNA_FASTA
    os.environ['IPD_GDNA_FASTA'] = IPD_GDNA_FASTA

    # make gdna and cdna files
    separateGdnaFromCdna(in_gbk, IPD_GDNA_FASTA, IPD_CDNA_FASTA)

    # map to reference
    vsearch = vsearchMapToReference(IPD_CDNA_FASTA, exemplar_fasta, 'tmp.sam')

    # calculate trim
    trim_list = trimQueryByMapping('tmp.sam')

    # apply trim to FASTA
    replaceGenbankSequence(in_gbk, trim_list, IPD_TRIMMED_FASTA)

    # deduplicate cDNA
    deduplicate_fasta(IPD_TRIMMED_FASTA, IPD_EXON2_FASTA)

    # merge gDNA and cDNA files
    subprocess.run('cat $IPD_GDNA_FASTA $IPD_EXON2_FASTA > $IPD_IMMUNOWES_FASTA', shell = True)
    
    # remove unnecessary intermediate files
    subprocess.run('rm $IPD_TRIMMED_FASTA', shell = True)

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()
