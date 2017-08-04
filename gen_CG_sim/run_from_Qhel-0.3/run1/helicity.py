#!/usr/bin/python

from math import sqrt
from math import exp

Q = []
var = 25.0
r0 = 0.5
Qi = 0.00
Ni = 4.0

filename = open("distance_bonds_1-4.xvg","r")
if filename:
   for line in filename:
       line = line.strip()
       Qi = exp( -var*(float(line.split()[1])-r0)**2 )
       Q.append(Qi)

i = 0
filename = open("distance_bonds_2-5.xvg","r")
if filename:
  for line in filename.readlines ():
    line = line.strip()
    Qi = exp( -var*(float(line.split()[1])-r0)**2 )
    Q[i] += Qi     
    i += 1

i = 0
filename = open("distance_bonds_3-6.xvg","r")
if filename:
  for line in filename.readlines ():
    line = line.strip()
    Qi = exp( -var*(float(line.split()[1])-r0)**2 )
    Q[i] += Qi
    i += 1

i = 0
filename = open("distance_bonds_4-7.xvg","r")
if filename:
  for line in filename.readlines ():
    line = line.strip()
    Qi = exp( -var*(float(line.split()[1])-r0)**2 )
    Q[i] += Qi
    i += 1

N = len(Q)
for i in range(0,N):
  
  Q[i] /= Ni
  print str(Q[i])

