#!/usr/bin/env python3

import argparse
import biotools as bt

# Write a program that masks areas of low complexity sequence
# Use argparse for command line arguments (see example below)
# Use read_fasta() from biotools.py

# setup
parser = argparse.ArgumentParser(
        description='Low complexity sequence masker.')

# required arguments
parser.add_argument('--input', required=True, type=str,
        metavar='<str>', help='fasta file')

# optional arguments with default parameters
parser.add_argument('--window', required=False, type=int, default=15,
        metavar='<int>', help='window size [%(default)i]')
parser.add_argument('--threshold', required=False, type=float, default=1.1,
        metavar='<float>', help='entropy threshold [%(default)f]')

# switches
parser.add_argument('--lowercase', action='store_true',
        help='report lower case instead of N')

# finalization
arg = parser.parse_args()

for name, seq in bt.read_fasta(arg.input):
    print(f'>{name}')
    output = []
    sequence = []
    in_low = False
    low_count = 0
    for i in range(len(seq)-arg.window+1):
        sequence = seq[i:i+arg.window]
        entropy = bt.shannon(sequence)
        if entropy <= arg.threshold:
            if in_low: low_count += 1
            else:
                low_count = arg.window
                in_low = True
        else:
            in_low = False
        if low_count > 0:
            if arg.lowercase:
                output.append(sequence[0].lower())
            else:
                output.append('N')
            low_count -= 1
        else:
            output.append(sequence[0])

    # build the remainder of the last window
    for i in range(1,arg.window):
        if low_count > 0:
            if arg.lowercase:
                output.append(sequence[i].lower())
            else:
                output.append('N')
            low_count -= 1
        else:
            output.append(sequence[i])

    # format 60 characters per line
    for i in range(0,len(output)-59,60):
        print(''.join(output[i:i+60]))
    remainder = len(output) % 60
    if remainder:
        print(''.join(output[-remainder:]))

"""
python3 entropy_filter.py --help
usage: entropy_filter.py [-h] --input <path> [--window <int>]
                         [--threshold <float>] [--lowercase]

Low complexity sequence masker.

optional arguments:
  -h, --help           show this help message and exit
  --input <path>       fasta file
  --window <int>       optional integer argument [15]
  --threshold <float>  entropy threshold [1.100000]
  --lowercase          report lower case instead of N


python3 entropy_filter.py --input genome.fa.gz | head -20
>I
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAAAAATTGAGATAAGAAAACATTTTACTTTTTCAAAATTGTTTTCATGC
TAAATTCAAAACNNNNNNNNNNNNNNNAAGCTTCTAGATATTTGGCGGGTACCTCTAATT
TTGCCTGCCTGCCAACCTATATGCTCCTGTGTTTAGGCCTAATACTAAGCCTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAAAAGAATATGGTAGCTACAGAAACGGTAGTACACTCTTCTGNNNNNNNNNNNNNN
NTGCAATTTTTATAGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAAT
TGCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAAACANNNNNNNNNNNNNNNGAAAT
GAATATCGTAGCTACAGAAACGGTTGTGCACTCATCTGAAANNNNNNNNNNNNNNNNNNN
NNGCACTTTGTGCAGAATTCTTGATTCTTGATTCTTGCAGAAATTTGCAAGAAAATTCGC
"""
