#!/usr/bin/env python3

# Write a program that prints out the position, frame, and letter of the DNA
# Try coding this with a single loop
# Try coding this with nested loops

dna = 'ATGGCCTTT'


"""
0 0 A
1 1 T
2 2 G
3 0 G
4 1 C
5 2 C
6 0 T
7 1 T
8 2 T
"""

#Frame is index modulo 3...

print("single loop")
for index in range(0,len(dna)):
  print (index, index % 3, dna[index])

# two loops... 
# NOTE: this does no range checking... only works if the
# length of dna is a multiple of 3.
#
# Outer loop iterates over index
# Inner loop iterates over frame

print("\nnested loops")
framesize = 3
for index in range(0, len(dna), framesize):
  for frameOffset in range(0, framesize):
    print (index+frameOffset, frameOffset, dna[index+frameOffset])
