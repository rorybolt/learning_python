#!/usr/bin/env python3

from math import sqrt
import fileinput

# Write a program that computes typical stats
# Count, Min, Max, Mean, Std. Dev, Median
# No, you cannot import any other modules!

data = []
sum = 0
for line in fileinput.input():
  if line.startswith('#'): continue   # skip comments
  f = float(line)
  sum += f
  data.append(f)	      # add to data array

data.sort()

count = len(data)
print(f'Count: {count}')
print(f'Minimum: {data[0]}')
print(f'Maximum: {data[count-1]}')
mean = sum/count
print(f'Mean: {mean}')

sum = 0
# compute standard deviation
for f in data:
  diff = f - mean
  diff *= diff
  sum += diff

sum /= count
sd = sqrt(sum)
print(f'Std. dev: {sd:.3f}')

# now compute median
median = 0
midpoint = int(count / 2) 
if (count % 2 == 1):
  median = data[midpoint]
else:
  median = (data[midpoint] + data[midpoint - 1]) / 2
print(f'Median: {median}')

"""
python3 stats.py numbers.txt
Count: 10
Minimum: -1.0
Maximum: 256.0
Mean: 29.147789999999997
Std. dev: 75.777
Median 2.35914
"""
