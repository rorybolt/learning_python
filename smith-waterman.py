#!/usr/bin/env python3


# Smith-Waterman algorithm for local alignment of nucleotide sequences.


import argparse
import biotools

# These default scores are taken from Wikipedia.
# en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
# setup
parser = argparse.ArgumentParser(
        description='Computes optimal pairwise alignment between two sequences.')

# required arguments
#parser.add_argument('--input', required=True, type=str,
#        metavar='<str>', help='fasta file')

# optional arguments with default parameters
parser.add_argument('--match', required=False, type=int, default=2,
        metavar='<int>', help='match value` [%(default)i]')
parser.add_argument('--mismatch', required=False, type=int, default=-1,
        metavar='<int>', help='mismatch value` [%(default)i]')
parser.add_argument('--gap', required=False, type=int, default=-1,
        metavar='<int>', help='gap value` [%(default)i]')

# finalization
arg = parser.parse_args()


def calc_score(matrix, x, y):
    if seq1[x - 1] == seq2[y - 1]:
        sim_score = arg.match
    else:
        sim_score = arg.mismatch

    diag_score = matrix[x - 1][y - 1] + sim_score
    up_score   = matrix[x - 1][y] + arg.gap
    left_score = matrix[x][y - 1] + arg.gap

    return max(0, diag_score, up_score, left_score)


def create_score_matrix(rows, cols):
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]

    # create the scoring matrix.
    max_score = 0
    max_pos   = None    # The coordinates of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)

            score_matrix[i][j] = score

    assert(max_pos is not None)

    return score_matrix, max_pos

# directions for matrix traversal
stop = 0
diag = 1
up = 2
left = 3

def next_move(score_matrix, x, y):
    diag_score = score_matrix[x - 1][y - 1]
    up_score   = score_matrix[x - 1][y]
    left_score = score_matrix[x][y - 1]
    # diag move is prefered for ties....
    if diag_score >= up_score and diag_score >= left_score:
        if diag_score != 0: return diag
        else:               return stop	
    # up move is preferred for ties...
    elif up_score > diag_score and up_score >= left_score:
        if up_score != 0: return up
        else:             return stop
    # only option left is left!
    else:
        if left_score != 0: return left
        else:               return stop


def traceback(score_matrix, start_pos):
    aligned1 = []
    aligned2 = []
    x, y         = start_pos
    move         = next_move(score_matrix, x, y)
    while move != stop:
        if move == diag:
            aligned1.append(seq1[x - 1])
            aligned2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == up:
            aligned1.append(seq1[x - 1])
            aligned2.append('-')
            x -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[y - 1])
            y -= 1

        move = next_move(score_matrix, x, y)

    aligned1.append(seq1[x - 1])
    aligned2.append(seq1[y - 1])

    return x-1,''.join(reversed(aligned1)), y-1,''.join(reversed(aligned2))



def alignment_string(aligned1, aligned2):
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for i in range(len(aligned1)):
        nt1=aligned1[i]
        nt2=aligned2[i]
        if nt1 == nt2:
            alignment_string.append('|')
            idents += 1
        elif nt1 == '-' or nt2 == '-':
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches




seq1     = "CTATCACCTGACCTCCAGGCCGATGCCCCTTCCGGC"
seq2     = "GCGAGTTCATCTATCACGACCGCGGTCG"

print("seq1: ",seq1)
print("seq2: ",seq2)

# Allow space for gap in the scoring matrix (e.g. dimension+1)
rows = len(seq1) + 1
cols = len(seq2) + 1

# Initialize the scoring matrix.
score_matrix, start_pos = create_score_matrix(rows, cols)

# Find the optimal path through the scoring matrix.
# This gives the optimal local alighnment
pos1, aligned1, pos2, aligned2 = traceback(score_matrix, start_pos)

# Print in BLAST style format...
alignment, idents, gaps, mismatches = alignment_string(aligned1, aligned2)
slen = len(aligned1)
print(f' Score = {(idents * arg.match)+(mismatches * arg.mismatch)+(gaps * arg.gap)}')  
print(f' Identities = {idents}/{slen} ({idents/slen:.1%}), Gaps = {gaps}/{slen} ({gaps/slen:.1%})')
print()
for i in range(0, slen, 60):
    slice1 = aligned1[i:i+60]
    slice2 = aligned2[i:i+60]
    print(f'Query  {pos1+i+1:<4}  {slice1}  {pos1+i+len(slice1):<4}')
    print(f'             {alignment[i:i+60]}')
    print(f'Sbjct  {pos2+i+1:<4}  {slice2}  {pos2+i+len(slice2):<4}')


