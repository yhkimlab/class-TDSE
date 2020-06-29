import numpy as np
from matplotlib import pyplot as plt
from modl import  Input


st=0
n = Input.ngx
nstep = Input.nstep

L=[]
file=open('wave.txt', 'r')

while (1):
    line=file.readline()

    try:escape=line.index('\n')
    except:escape=len(line)

    if line:
        L.append(line[0:escape])
    else:
        break

file.close()

L2=[]
file=open('Potential.txt', 'r')

while (1):
    line=file.readline()

    try:escape=line.index('\n')
    except:escape=len(line)

    if line:
        L2.append(line[0:escape])
    else:
        break

file.close()

tmp=L2[0].split()
a=np.size(tmp)
pot=np.zeros(a)
for i in range(0, a):
    pot[i]=tmp[i]

x = np.linspace(st, Input.Lx, n-1)
nL = np.size(L)
store = np.zeros((nL-1,n))
a = np.max(pot[5:np.size(pot)-5])

for i in range(1, nL):
    tmp=L[i].split()
    for j in range(0, n):
        store[i-1,j]=np.float64(tmp[j])

for i in range(0,nstep):
     plt.plot(x, store[i*1,1:n], label= 'Time = %.6f fs' %(store[i*1,0]))
     if a != 0 :
     	plt.plot(x, pot[0:n-1]/a*0.15, label='Potential = %.3f eV ' %(a*27.211))
     plt.xlim(0,Input.Lx)
     plt.ylim(0,0.0052)
     plt.legend()
     plt.yticks([], [])
     plt.xlabel('Box [Angstrom]')
     plt.savefig('case_%03d.png' %i)
     plt.clf()

import os, sys

### Generate a movie from PNG files ###

os.system('convert -delay 10 -loop 1000 -quality 100 -resize 500x500 *.png movie.gif')
os.system('mkdir pngfile')
os.system('mv case_* pngfile')
