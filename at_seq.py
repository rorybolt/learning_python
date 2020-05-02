#!/usr/bin/env python3

import random
#random.seed(1) # comment-out this line to change sequence each time

# Write a program that stores random DNA sequence in a string
# The sequence should be 30 nt long
# On average, the sequence should be 60% AT
# Calculate the actual AT fraction while generating the sequence
# Report the length, AT fraction, and sequence

dna = []
at = 0
for i in range(30):
    r = random.random()
    if r < 0.6:
        r = random.random()
        if r < 0.5: dna.append('A')
        else      : dna.append('T')
        at += 1
    else:
        r = random.random()
        if r < 0.5: dna.append('C')
        else      : dna.append('G')
i += 1
print(i, at/i, ''.join(dna))
"""
30 0.6666666666666666 ATTACCGTAATCTACTATTAAGTCACAACC
"""
