from module import *
# module = reload(module)
import sys
import csv
import numpy as np

'''
From arXiv:1511.07428
R: the number of species
r_x: number of occurrences of specie x
n: number of samples
counts: represents the sampled population distribution, assuming an hypergeometric model
Phi: array, whose i-th location counts how many species have been seen i times
t: fold-change from the initial sampling of size n (at most log n)
sd_samples: numer of times to sample the smoothing distro
'''

n = int(sys.argv[1])
t = int(sys.argv[2])
sd_samples = int(sys.argv[3])
step = float(sys.argv[4])

theFile = open("prevalences.tsv", "r")
Phi = []
for val in theFile.read().split():
    Phi.append(int(val))
theFile.close()

x = []
y = []

for tvar in np.arange(0.0, t, step):
    x.append(tvar)
    if tvar == 0: value = sum(Phi)
    if tvar > 0: value = y[0] + smoothedGT(Phi,n,tvar,sd_samples)
    y.append(value)

outfile = open('accumulation_curve.tsv', 'wa') 
for item in zip(x,y):
  print>>outfile, item

print 'Accumulation curve generated!'

