#!/usr/bin/env python3

import random
import sys

# Write a program that simulates random BAC coverage over a genome
# Command line arguments include
# 	Genome size (e.g. 1000)
# 	X coverage (e.g. 5)
# Use assert() to check parameter bounds
# Report min, max, and histogram of coverage
# Note that your output may vary due to random function

assert(len(sys.argv) == 3)
size = int(sys.argv[1])
coverage = int(sys.argv[2])

assert(size > 0)
assert(coverage > 0)

genome = [0] * size

for  i in range (size * coverage):
    r = random.randint(0,size-1)
    genome[r] += 1

print(f'Size: {size}')
print(f'X: {coverage:.1f}')
print(f'BACs: {size*coverage}')

genome.sort()
min = genome[0]
max = genome[size-1]
print(f'Min: {min}')
print(f'Max: {max}')
print("Counts:")

histogram = [0] * (max+1)
for i in range(size):
    histogram[genome[i]] += 1

for i in range(max+1):
    print(i, histogram[i])

"""
Size: 1000
X: 5.0
BACs: 5000
Min: 0
Max: 13
Counts:
0 5
1 39
2 88
3 144
4 175
5 150
6 151
7 116
8 59
9 40
10 20
11 5
12 6
13 2
"""
