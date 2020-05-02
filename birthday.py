#!/usr/bin/env python3

import sys
import random

# You are probably well aware of the 'birthday paradox'
# https://en.wikipedia.org/wiki/Birthday_problem
# Let's try simulating it
# We will have a variable number of bins (can be months or days)
# And some number of trials for the simulation
# And some number of people whose have random birthdays
# Use assert() to check parameters
# On the command line:
#	python3 birthday.py <bins> <trials> <people>

assert(len(sys.argv) == 4)
bins = int(sys.argv[1])
trials = int(sys.argv[2])
people = int(sys.argv[3])
assert(bins > 0)
assert(trials > 0)
assert(people > 0)

collisions = 0
for i in range(0, trials):
    calendar = [0] * bins
    for j in range(people):
        birthday = random.randint(0,bins-1)
        if calendar[birthday] > 0:
            collisions += 1
            break
        else:
            calendar[birthday] += 1
print(collisions/trials)

"""
python3 birthday.py 365 1000 23
0.520
"""

