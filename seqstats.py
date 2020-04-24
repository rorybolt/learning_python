#!/usr/bin/env python3

import fileinput

# Write a program that computes typical sequence stats
# No, you cannot import any other modules!
# Use rand_seq to generate the sequences
# Expected output is shown below

count=0;
lengths = []
sum = 0
na = 0 	# number of A's
nc = 0 	# number of C's
ng = 0 	# number of G's
nt = 0 	# number of T's

for line in fileinput.input():
  if line.startswith('>'): continue   # skip names
  i = len(line)
  lengths.append(i)
  sum += i
  count += 1
  for j in range(i):
    if line[j] == 'A'  : na += 1
    elif line[j] == 'T': nt += 1
    elif line[j] == 'G': ng += 1
    else               : nc += 1

lengths.sort()

# compute N50...
N50sum = sum
for i in range(count):
  N50sum -= lengths[count-i-1]
  if N50sum <= sum/2 : break
N50 = lengths[count-i-1]
  
print(f'Number of sequences: {count}')
print(f'Number of letters: {sum}')
print(f'Minimum length: {lengths[0]}')
print(f'Maximum length: {lengths[count-1]}')
print(f'N50: {N50}')
print(f'Composition: A={na/sum:.3f} C={nc/sum:.3f} G={ng/sum:.3f} T={nt/sum:.3f}')

"""
python3 rand_seq.py 100 100 100000 0.35 | python3 seqstats.py
Number of sequences: 100
Number of letters: 4957689
Minimum length: 219
Maximum length: 99853
N50: 67081
Composition: A=0.325 C=0.175 G=0.175 T=0.325
"""
