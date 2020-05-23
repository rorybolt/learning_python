#!/usr/bin/env python3

# Modify entropy_fast() however you like to make it faster
# Ideally, your method is faster at all ranges of window size

import math
import time
import random

def entropy_slow(seq, w, th):
	t0 = time.perf_counter()
	low_H_count = 0
	
	for i in range(len(seq) - w + 1):
		win = seq[i:i+w]
		a, c, g, t = 0, 0, 0, 0
		for nt in win:
			if   nt == 'A': a += 1
			elif nt == 'C': c += 1
			elif nt == 'G': g += 1
			elif nt == 'T': t += 1
		total = a + c + g + t
		h = 0
		pa, pc, pg, pt = a/total, c/total, g/total, t/total
	
		if a != 0: h -= pa * math.log2(pa)
		if c != 0: h -= pc * math.log2(pc)
		if g != 0: h -= pg * math.log2(pg)
		if t != 0: h -= pt * math.log2(pt)
	
		if h < th: low_H_count += 1
	
	t1 = time.perf_counter()
	return low_H_count, t1-t0

#	"#EveryNanosecondCounts!"
# Overall strategies:
#	Avoid redoing any work already performed
#		Dictionary of pre computed values...
#		log2() is expensive!
#	Avoid unnecessary statements inside loops
#		The parameters to "for" are evaluated every loop!
#		Two for loops are used to avoid if i < w : continue
#		It is faster to unwind the ha computation loop
def entropy_fast(seq, w, th):
	t0 = time.perf_counter()
	low_H_count = 0
	num = {"A":0,"C":0,"G":0,"T":0} # dictionary of observed counts
	has = {"A":0,"C":0,"G":0,"T":0} # dictionary of computed HA values

	# build the first window
	for i in range(w): num[seq[i]] += 1

	# now compute the HA value for the 1st window
	h = 0
	# Note: I timed a clean "for nt in num:" loop and it is slower
	# so this is unwound...
	# check to see if we have already computed a HA value for a
	# probability: if so, reuse, else compute and add to dictionary
	a, c, g, t = num['A'], num['C'], num['G'], num['T']
	if a != 0: 
		if a in has: h -= has[a]
		else:
			p = a/w
			hp = p * math.log2(p)
			has[a] = hp
			h -= hp
	if c != 0: 
		if c in has: h -= has[c]
		else:
			p = c/w
			hp = p * math.log2(p)
			has[c] = hp
			h -= hp
	if g != 0: 
		if g in has: h -= has[g]
		else:
			p = g/w
			hp = p * math.log2(p)
			has[g] = hp
			h -= hp
	if t != 0: 
		if t in has: h -= has[t]
		else:
			p = t/w
			hp = p * math.log2(p)
			has[t] = hp
			h -= hp
	if h < th: low_H_count += 1

	# this surprised me... I measured a difference with the computation
	# in the for statement. Since seq is not altered, len(seq) only
	# needs to be computed once
	# gcc would have done this for you...
	l = len(seq) - w + 1
	for i in range(w,l):		# advance the window
		num[seq[i]] += 1	# entering window
		num[seq[i - w]] -= 1	# leaving window

		# Does all the following look familiar? it is only needed
		# to avoid a single "if i < w: continue" statement...
		# but this loop executes a lot so it matters.

		# now compute the HA value for the current window
		h = 0
		# Note: I timed a clean "for nt in num:" loop and it is slower
		# so this is unwound...
		# check to see if we have already computed a HA value for a
		# probability: if so, reuse, else compute and add to dictionary
		
		a, c, g, t = num['A'], num['C'], num['G'], num['T']
		if a != 0: 
			if a in has: h -= has[a]
			else:
				p = a/w
				hp = p * math.log2(p)
				has[a] = hp
				h -= hp
		if c != 0: 
			if c in has: h -= has[c]
			else:
				p = c/w
				hp = p * math.log2(p)
				has[c] = hp
				h -= hp
		if g != 0: 
			if g in has: h -= has[g]
			else:
				p = g/w
				hp = p * math.log2(p)
				has[g] = hp
				h -= hp
		if t != 0: 
			if t in has: h -= has[t]
			else:
				p = t/w
				hp = p * math.log2(p)
				has[t] = hp
				h -= hp

		if h < th: low_H_count += 1
	t1 = time.perf_counter()
	return low_H_count, t1-t0

# random.seed(0)
# create a random chromosome
seq = []
alph = ['A', 'C', 'G', 'T']
for i in range(int(1e5)):
	seq.append(alph[random.randint(0,3)])
seq = ''.join(seq)

# test speed at various word sizes
W = [7, 15, 100]
for w in W:
	cs, ts = entropy_slow(seq, w, 1)
	cf, tf = entropy_fast(seq, w, 1)
	assert(cs == cf)
	print(tf / ts)
"""
	atf = tf
	for i in range(49):
		cf, tf = entropy_fast(seq, w, 1)
		atf += tf
	print((atf/50) / ts)
"""
