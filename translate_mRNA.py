#!/usr/bin/env python3

import gzip
import sys
import biotools as bt
import argparse

# Use argparse
# Write a program that translates an mRNA
# Assume the protein encoded is the longest ORF

parser = argparse.ArgumentParser(
        description='Translates an mRNA.')
# required arguments
parser.add_argument('--file', required=True, type=str,
        metavar='<str>', help='path to protein file')

# finalization
arg = parser.parse_args()

def longest_orf(seq):
    # find all ATGs
    atgs = []
    for i in range(len(seq)-2):
        if seq[i:i+3] == 'ATG': atgs.append(i)
    
    # for each atg, find nearest in frame stop; find longest...
    max_len = 0
    max_seq = None
    for atg in atgs:
        stop = None
        for i in range(atg, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                stop = i
                break
        if stop != None:
            cds_len = stop - atg + 3
            if cds_len > max_len:
                max_len = cds_len
                max_seq = seq[atg:atg+cds_len]

    if max_seq == None: return None

    # translate longest ORF into protein
    return bt.translate(max_seq)

for name, seq in bt.read_fasta(arg.file):
    protein = longest_orf(seq)
    if protein != None:
        print(f'>{name}')
        print(protein)
"""
python3 translate_mRNA.py --file ../Lesson05/transcripts.fasta.gz
>CBG00001.1
MTFCENKNLPKPPSDRCQVVVISILSMILDFYLKYNPDKHWAHLFYGASPILEILVIFGMLANSVYGNKLAMFACVLDLVSGVFCLLTLPVISVAENATGVRLHLPYISTFHSQFSFQVSTPVDLFYVATFLGFVSTILILLFLILDALKFMKLRKLRNEDLEKEKKMNPIEKV*
>CBG00006.1
MNGVEKVNKYFDIKDKRDFLYHFGFGVDTLDIKAVFGDTKFVCTGGSPGRFKLYAEWFAKETSIPCSENLSRSDRFVIYKTGPVCWINHGMGTPSLSIMLVESFKLMHHAGVKNPTFIRLGTSGGVGVPPGTVVVSTGAMNAELGDTYVQVIAGKRIERPTQLDATLREALCAVGKEKNIPVETGKTMCADDFYEGQMRLDGYFCDYEEEDKYAFLRKLNSLGVRNIEMESTCFASFTCRAGFPSAIVCVTLLNRMDGDQVQIDKEKYIEYEERPFRLVTAYIRQQTGV*
etc.
"""
