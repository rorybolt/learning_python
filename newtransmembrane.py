#!/usr/bin/env python3

import argparse
import biotools as bt

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


# setup
parser = argparse.ArgumentParser(
	description='Predicts transmembrane proteins.')
# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<str>', help='file name of protein file')

# optional arguments with default parameters
parser.add_argument('--win1', required=False, type=int, default=8,
	metavar='<int>', help='length of signal peptide [%(default)i]')
parser.add_argument('--win2', required=False, type=int, default=11,
	metavar='<int>', help='length of trans-membrane region [%(default)i]')
parser.add_argument('--kd1', required=False, type=float, default=2.5,
	metavar='<float>', help='kd value for signal peptide [%(default)f]')
parser.add_argument('--kd2', required=False, type=float, default=2.0,
	metavar='<float>', help='kd value for hydrophobic region [%(default)f]')

# finalization
arg = parser.parse_args()

for name, seq in bt.read_fasta(arg.file):
    cond1 = False
    cond2 = False
    # test for signal condition KD > 2.5
    for i in range(0,(30 - arg.win1) + 1):
        assert(i < 23)
        if bt.computeKD(seq, i, arg.win1) > arg.kd1:
            cond1 = True
            break
    if not cond1: continue


    # test for hydrophobic condition KD > 2.0 and no peptide
    for i in range(30,(len(seq) - arg.win2)): # leave out trailing *
        if bt.computeKD(seq, i, arg.win2) > arg.kd2:
            if 'P' in seq[i:i+arg.win2]:
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
