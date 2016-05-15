'''
This program calculates the magnetization of the 2D Ising model for various
values of the parameter Jkbt = J / (k_B*T). The code runs through all values
in fList for one value in N, records the magnetization at equilibrium for 
each fList value, then moves onto the next value in N. The plot is of len(N)
lines each one composed of len(fList) points modelling the change in
magnetization for values of Jkbt and N.
'''
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import product

N = [5, 10, 20, 30, 50]             ## Lattice Sizes
Tc = 0.5 * np.log(1 + np.sqrt(2))   ## Critical Temperature
SW = 10000                          ## Number of Sweeps

avgMagTot = []
avgSusTot = []
fList = []

## Creates the list used to plot the data.
def fListDefine():
    f = 0.6
    fList = []
    while f < 1.41:
        fList.append(f)
        f +=0.05
    return fList

## Runs one sweep through the lattice and records and returns the magnetization
def Ising2D(lattice, Jkbt, N):
    magnetization = 0

    for p in product(range(N), range(N)):
        i, j = p
        energy = (lattice[i][(j - 1) % N] + lattice[i][(j + 1) % N] + \
                          lattice[(i - 1) % N][j] + lattice[(i + 1) % N][j]) * \
                          lattice[i][j]
        magnetization += lattice[i][j]

        ## Metropolis update
        if energy <= 0 or np.random.random() < np.exp(-2 * Jkbt * energy):
            lattice[i][j] *= -1

    return lattice, magnetization

## For a lattice size L, runs SW sweeps of the lattice, recording the 
## magnetization, then increases the T/Tc value from 0.6-1.4
def MagVsT(N):
    factor = 0.6
    lattice = np.ones(N**2).reshape((N,N))
    avgMagTemps = []
    avgSusTemps = []

    while factor < 1.41:
        magTotSweeps = []
        magTotSqSweeps = []

        ## Creates LxL array of random distribution of spin up and spin down.
        for i, j in product(range(N), range(N)):
            lattice[i][j] = random.randint(0, 1)
            if lattice[i][j] == 0:
                lattice[i][j] = -1

        ## Runs SW sweeps and stores the absolute magnetization from each
        ## sweep in magList
        for i in range(SW):
            lattice, magnetization = Ising2D(lattice, Tc / factor, N)
            if i > 3 * SW / 4:
                magTotSweeps.append(abs(magnetization))
                magTotSqSweeps.append(magnetization**2)

        avgMag = np.average(magTotSweeps)
        avgSus = (np.average(magTotSqSweeps) - np.average(magTotSweeps)**2) / (Tc * factor)

        avgMagTemps.append(avgMag / N**2)
        avgSusTemps.append(avgSus / N**2)

        ## Increase the T/Tc values by 0.05
        factor += 0.05

    return avgMagTemps, avgSusTemps

## Runs the functions for various lattice size in the list N
for i in N:
    avgMagTemps, avgSusTemps = MagVsT(i)
    print('Done with N = %d' %i)

    ## All of the data is stored in these two lists
    avgSusTot.append(avgSusTemps)
    avgMagTot.append(avgMagTemps)

## Saves the data collected as .txt files.
with open('avgMagTot.txt', 'w') as f:
    for i in range(len(avgMagTot)):
        f.write('N = %d:\n' % N[i])
        f.write(', '.join(str(x) for x in avgMagTot[i]))
        f.write('\n')
with open('avgSusTot.txt', 'w') as f:
    for i in range(len(avgSusTot)):
        f.write('N = %d:\n' % N[i])
        f.write(", ".join(str(x) for x in avgSusTot[i]))
        f.write('\n')

fList = fListDefine()

## Plots the absolute magnetization versus temperature (T/Tc) for each of the
## lattice sizes
plt.figure()
for i in range(len(avgMagTot)):
    plt.plot(fList, avgMagTot[i], label = r'$N=%d$' % N[i], lw = '2')
plt.ylim(0, 1)
plt.xlim(0.6, 1.4)
plt.title('Magnetization for $N$x$N$ Lattices at Different Temperatures' % N)
plt.legend()
plt.xlabel(r'$T/T_c$')
plt.ylabel('Absolute Magnetization per Site')
plt.show()

## Plots the susceptibility versus temperature (T/Tc) for each of the
## lattice sizes
plt.figure()
for i in range(len(avgSusTot)):
    plt.plot(fList, avgSusTot[i], label = r'$N=%d$' % N[i], lw = '2')
plt.xlim(0.6, 1.4)
plt.title('Susceptibility for $N$x$N$ Lattices at Different Temperatures' % N)
plt.legend()
plt.xlabel(r'$T/T_c$')
plt.ylabel('Susceptilibity per Site')
plt.show()