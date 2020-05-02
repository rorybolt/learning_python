#!/usr/bin/env python3

import math
import fileinput

# Write a Shannon entropy calculator: H = -sum(pi * log(pi))
# Use fileinput to get the data from nucleotides.txt
# Make sure that the values are probabilities
# Make sure that the distribution sums to 1
# Report with 3 decimal figures

sum = 0
count = 0
data = []
for line in fileinput.input():
  if line.startswith('#'): continue   # skip comments
  f = float(line[1:])
  assert( f >= 0 and f <= 1)
  sum += f
  data.append(f)              # add to data array
  count += 1

assert(math.isclose(sum, 1))

h = 0
for i in range(count):
    h -= data[i] * math.log(data[i],2)
 
print(f'{h:.3f}')


"""
python3 entropy.py nucleotides.txt
1.846
"""
