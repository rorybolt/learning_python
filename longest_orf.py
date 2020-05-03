#!/usr/bin/env python3

import gzip
import sys

# Write a program that finds all of the open reading frames in a transcript
# ORFs start with ATG and end with a stop codon (TAA, TAG, TGA)
# Command line: python3 orfs.py
# Output:
#	Gene_name seq_length cds_length

def read_fasta(filename):
    name = None
    seqs = []

    fp = None
    if filename == '-':
        fp = sys.stdin
    elif filename.endswith('.gz'):
        fp = gzip.open(filename, 'rt')
    else:
        fp = open(filename)

    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith('>'):
            if len(seqs) > 0:
                seq = ''.join(seqs)
                yield(name, seq)
                name = line[1:]
                seqs = []
            else:
                name = line[1:]
        else:
            seqs.append(line)
    yield(name, ''.join(seqs))
    fp.close()

def getOrfLen(seq):
    orfs = []
    for f in range(3):
        atgs = []
        for i in range(f, len(seq) - f - 2, 3):
            codon = seq[i:i+3]
            if codon == 'ATG':
                atgs.append(i)
        for start in atgs:
            for i in range(start, len(seq) - f - 2, 3):
                codon = seq[i:i+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orfs.append(i - start + 3)
#                    print(i-start+3,seq[start:i+3])
                    break
    orfs.sort()
    return orfs[-1]

assert(len(sys.argv)==2)

for name, seq in read_fasta(sys.argv[1]):
    print(name, len(seq), getOrfLen(seq))

"""
python3 longest_orf.py transcripts.fasta.gz  | head
CBG00001.1 720 525
CBG00006.1 1310 870
CBG00032.1 1985 1737
CBG00035.1 1867 1716
CBG00041.1 1075 891
CBG00049.1 693 189
CBG00050.1 1797 1245
CBG00056.1 1545 1392
CBG00059.1 1567 1488
CBG00060.1 1488 1296

python3 fasta_stats.py transcripts.fasta.gz fasta_stats.py -
Count: 232
Total: 278793
Min: 603
Max: 1991
Mean: 1201.7
NTs: 0.291 0.218 0.210 0.281

python3 rand_fasta.py 232 603 1991 0.291 0.218 0.210 0.281 | python3 fasta_stats.py -
Count: 232
Total: 311538
Min: 605
Max: 1990
Mean: 1342.8
NTs: 0.291 0.219 0.210 0.280

python3 rand_fasta.py 232 603 1991 0.291 0.218 0.210 0.281 | python3 longest_orf.py - | head
seq-0 1833 165
seq-1 1759 387
seq-2 1677 81
seq-3 1765 144
seq-4 1347 135
seq-5 1141 180
seq-6 1127 240
seq-7 1532 162
seq-8 1376 78
seq-9 793 213
"""
