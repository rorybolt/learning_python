#!/usr/bin/env python3

import gzip
import sys

# Write a program that reports the longest coding sequence
# Translate the sequence to amino acids
# Use a generator function somewhere in your program
# Check both strands
# See below for command line and expected output

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

AAtable = dict({
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'F', 'CTC':'F', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
})

def printProtein(cds, seq):
    protein = "" 
    for i in range(cds[1],cds[1]+cds[0],3):
        codon = seq[i:i+3]
        protein += AAtable[codon]
    print(protein)

def getStart(seq, frame):
    for i in range(frame, len(seq) -frame -2, 3):
        codon = seq[i:i+3]
        if codon == 'ATG':
            yield i

def getLongestCDS(seq):
    orfs = []
    
    for f in range(3):
        for start in getStart(seq, f):
            for i in range(start, len(seq) - f - 2, 3):
                codon = seq[i:i+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                    orfs.append((i - start + 3, start))
#                    print(i-start+3,seq[start:i+3])
                    break
    orfs.sort()
    return orfs[-1]

assert(len(sys.argv)==2)

inverse = dict({ord('A'):'T',
                ord('C'):'G',
                ord('G'):'C',
                ord('T'):'A'})

for name, seq in read_fasta(sys.argv[1]):
    print(f'>{name}')
    cds = getLongestCDS(seq)
    partner = seq.translate(inverse)
    partner_cds = getLongestCDS(partner)
    if partner_cds[0] > cds[0]:
        printProtein(partner_cds, partner)
    else:
        printProtein(cds,seq)

"""
python3 longest_orf.py transcripts.fasta.gz
>CBG00001.1
MTFCENKNLPKPPSDRCQVVVISILSMILDFYLKYNPDKHWAHLFYGASPILEILVIFGMLANSVYGNKLAMFACVLDLVSGVFCLLTLPVISVAENATGVRLHLPYISTFHSQFSFQVSTPVDLFYVATFLGFVSTILILLFLILDALKFMKLRKLRNEDLEKEKKMNPIEKV*
>CBG00006.1
MNGVEKVNKYFDIKDKRDFLYHFGFGVDTLDIKAVFGDTKFVCTGGSPGRFKLYAEWFAKETSIPCSENLSRSDRFVIYKTGPVCWINHGMGTPSLSIMLVESFKLMHHAGVKNPTFIRLGTSGGVGVPPGTVVVSTGAMNAELGDTYVQVIAGKRIERPTQLDATLREALCAVGKEKNIPVETGKTMCADDFYEGQMRLDGYFCDYEEEDKYAFLRKLNSLGVRNIEMESTCFASFTCRAGFPSAIVCVTLLNRMDGDQVQIDKEKYIEYEERPFRLVTAYIRQQTGV*
etc.
"""
