#!/usr/bin/env python3

import gzip
import sys

# Write a program that predicts if a protein is trans-membrane
# Trans-membrane proteins have the following properties
#	Signal peptide: https://en.wikipedia.org/wiki/Signal_peptide
#	Hydrophobic regions(s): https://en.wikipedia.org/wiki/Transmembrane_protein
#	No prolines (alpha helix)
# Hydrophobicity is measued via Kyte-Dolittle
#	https://en.wikipedia.org/wiki/Hydrophilicity_plot
# For our purposes:
#	Signal peptide is 8 aa long, KD > 2.5, first 30 aa
#	Hydrophobic region is 11 aa long, KD > 2.0, after 30 aa


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

KDtable = dict({'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5})

# Compute the average KD value
def computeKD(seq, start, len):
    kd = 0
    for i in range(len):
        val = KDtable.get(seq[start+i])
        if val is None:
#            print("Bad data: ", seq[start:start+len])
            return 0 # Bad data...
        kd += val
    return kd/len

for name, seq in read_fasta("proteins.fasta.gz"):
    cond1 = False
    cond2 = False
    # test for signal condition KD > 2.5
    for i in range(0,(30 - 8) + 1):
        assert(i < 23)
        if computeKD(seq, i, 8) > 2.5:
            cond1 = True
            break
    if not cond1: continue


    # test for hydrophobic condition KD > 2.0 and no peptide
    for i in range(30,(len(seq) - 11)): # leave out trailing *
        if computeKD(seq, i, 11) > 2.0:
            if 'P' in seq[i:i+11]:
                continue
            cond2 = True
            break
    if cond2:
        print(name) 

"""
18w
Dtg
Krn
Lac
Mcr
PRY
Pxt
Pzl
QC
Ror
S1P
S2P
Spt
apn
bai
bdl
bou
bug
cue
drd
ft
grk
knk
ksh
m
nac
ort
rk
smo
thw
tsg
waw
zye
"""
