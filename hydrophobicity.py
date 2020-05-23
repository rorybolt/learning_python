#!/usr/bin/env python3
import argparse
import biotools as bt

# Write a program that computes hydrophobicity in a window
# Let the user choose the method (see below)
# https://en.wikipedia.org/wiki/Hydrophilicity_plot
# https://en.wikipedia.org/wiki/Hydrophobicity_scales

# setup
parser = argparse.ArgumentParser(
        description='Compute hydrophobicity in a window.')

# required arguments
parser.add_argument('--input', required=True, type=str,
        metavar='<str>', help='fasta file')

# optional arguments with default parameters
parser.add_argument('--window', required=False, type=int, default=15,
        metavar='<int>', help='window size [%(default)i]')
parser.add_argument('--method', choices=['kd','is','os','ois','cc'],
        required=False, type=str, default='kd',
        metavar='<str>', help='calculation method (kd, is, os, ois, cc) [%(default)s]')

# finalization
arg = parser.parse_args()

KDtable = {
    'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,
    'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,
    'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,
    'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5
}
IStable = {
    'I':-0.31,'V':0.07,'L':-0.56,'F':-1.13,'C':-0.24,
    'M':-0.23,'A':0.17,'G':0.01,'T':0.14,'S':0.13,
    'W':-1.85,'Y':-0.94,'P':0.45,'H':0.96,'E':2.02,
    'Q':0.58,'D':1.23,'N':0.42,'K':0.99,'R':0.81
}
OStable = {
    'I':-1.12,'V':-0.46,'L':-1.25,'F':-1.71,'C':-0.02,
    'M':-0.67,'A':0.50,'G':1.15,'T':0.25,'S':0.46,
    'W':-2.09,'Y':-0.71,'P':0.14,'H':2.33,'E':3.63,
    'Q':0.77,'D':3.64,'N':0.85,'K':2.80,'R':1.81
}
OIStable = {
    'I':-0.81,'V':-0.53,'L':-0.69,'F':-0.58,'C':0.22,
    'M':-0.44,'A':0.33,'G':1.14,'T':0.11,'S':0.33,
    'W':-0.24,'Y':0.23,'P':-0.31,'H':1.37,'E':1.61,
    'Q':0.19,'D':2.41,'N':0.43,'K':1.81,'R':1.00
}
CCtable = {
    'I':-0.528,'V':-0.308,'L':-0.342,'F':-0.370,'C':0.081,
    'M':-0.324,'A':-0.495,'G':0.386,'T':0.853,'S':0.963,
    'W':-0.270,'Y':1.667,'P':-0.322,'H':2.029,'E':3.173,
    'Q':2.176,'D':9.573,'N':2.354,'K':2.101,'R':4.383
}

def computeHD(seq,start,len,scale):
    hd = 0
    if scale == 'kd':
        dict = KDtable
    elif scale == 'is':
        dict = IStable
    elif scale == 'os':
        dict = OStable
    elif scale == 'ois':
        dict = OIStable
    elif scale == 'cc':
        dict = CCtable
    else:
        print(f'unsupported scale {scale}')
        return None

    for i in range(len):
        val = dict.get(seq[start+i])
        if val is None:
            print("Bad data: ", seq[start:start+len])
            return None # Bad data...
        hd += val
    return hd/len

for name, seq in bt.read_fasta(arg.input):
    print(f'>{name}')
    num_below = 0
    for i in range(0,len(seq)-arg.window):
        hd = computeHD(seq, i, arg.window, arg.method)
        if hd == None:
            continue
        if hd < 1:
            num_below += 1
    print(num_below)


"""
python3 hydrophobicity.py --input proteins.fasta.gz --window 11 --method kd
"""
