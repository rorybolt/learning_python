#!/usr/bin/env python3

import sys
import gzip
import random
import math

def read_fasta(filename):
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)

def randseq(l,gc):
    dna=[]
    for i in range(l):
        r=random.random()
        if r < gc:
            r = random.random()
            if r < 0.5 : dna.append('G')
            else       : dna.append('C')
        else:
            r = random.random()
            if r < 0.5 : dna.append('A')
            else       : dna.append('T')
    return(''.join(dna))

KDtable = {
    'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,
    'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,
    'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,
    'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5
}

# Compute the average KD value
def computeKD(seq, start, len):
    kd = 0
    for i in range(len):
        val = KDtable.get(seq[start+i])
        if val is None:
#            print("Bad data: ", seq[start:start+len])
            return 0 # Bad data...
        kd += val
    return kd/len

AAtable = dict({
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'F', 'CTC':'F', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
})

def translate(seq):
    protein = []
    for i in range(0,len(seq)-2,3):
        codon = seq[i:i+3]
        protein.append(AAtable[codon])
    return ''.join(protein)

def skew(seq):
    #(g-c/(g+c)
    g = 0
    c = 0
    for nt in seq:
        if nt == 'G'  : g+=1
        elif nt == 'C': c+=1
    return (g-c)/(g+c)


def shannon(seq):
    # cout frequencies in seq
    a = 0
    t = 0
    g = 0
    c = 0
    for nt in seq:
        if nt == 'G'  : g+=1
        elif nt == 'C': c+=1
        elif nt == 'A': a+=1
        elif nt == 'T': t+=1
    h = 0
    cnt = len(seq)
    if a > 0:
        h -= (a/cnt) * math.log(a/cnt, 2)
    if t > 0:
        h -= (t/cnt) * math.log(t/cnt, 2)
    if c > 0:
        h -= (c/cnt) * math.log(c/cnt, 2)
    if g > 0:
        h -= (g/cnt) * math.log(g/cnt, 2)
    return h
