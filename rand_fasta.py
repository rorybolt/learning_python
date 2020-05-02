#!/usr/bin/env python3

import gzip
import sys
import math
import random

# Write a program that creates random fasta files
# Create a function that makes random DNA sequences
# Parameters include length and frequencies for A, C, G, T
# Command line: python3 rand_fasta.py <count> <min> <max> <a> <c> <g> <t>

def rand_dna(length, pa, pc, pg, pt):
    dna = []
    for i in range(length):
        r = random.random()
        if r < pa:		dna.append('A')
        elif r < pa+pc:		dna.append('C')
        elif r < pa+pc+pg:	dna.append('G')
        else:			dna.append('T')
    return ''.join(dna)

assert(len(sys.argv) == 8)
count = int(sys.argv[1])
assert(count > 0)
min = int(sys.argv[2])
assert(min > 0)
max = int(sys.argv[3])
assert(max >= min)
pa = float(sys.argv[4])
assert(pa >=0 and pa <= 1)
pc = float(sys.argv[5])
assert(pc >=0 and pc <= 1)
pg = float(sys.argv[6])
assert(pg >=0 and pg <= 1)
pt = float(sys.argv[7])
assert(pt >=0 and pa <= 1)
assert(math.isclose(1,pa+pc+pg+pt))

for i in range(count):
    r = random.randint(min,max)
    dna = rand_dna(r, pa, pc, pg, pt)
    print(f'>{i}')
    print(dna)

"""

"""
