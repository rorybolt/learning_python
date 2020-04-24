#!/usr/bin/env python3

import random
import sys

# Create a program that generates random sequences in FASTA format
# Each name should be unique
# Length should have a minimum and maximum
# GC% should be a parameter
# Use assert() to check bounds of command line values
# When creating sequences, append and join
# Command line:
#	python3 rand_seq.py <# of seqs> <min> <max> <gc>

assert(len(sys.argv) == 5)
num = int(sys.argv[1])
min = int(sys.argv[2])
max = int(sys.argv[3])
gc = float(sys.argv[4])

assert(num > 0)
assert(min > 0)
assert(max > min)
assert(gc >= 0 and gc <= 1)

for i in range(num):
  print(f'>seq-{i}')
  len = random.randint(min, max)
  seq = []
  for j in range(len):
    r = random.random()
    if r < gc:
      r = random.random()
      if r < 0.5: nt = 'G'
      else      : nt = 'C'
    else:
      r = random.random()
      if r < 0.5: nt = 'A'
      else      : nt = 'T'
    seq.append(nt)
  print(''.join(seq))

"""
python3 rand_seq.py 3 10 20 0.5
>seq-0
GCGCGACCTTAT
>seq-1
ATCCTAGAAGT
>seq-2
CTTCGCTCGTG
"""

