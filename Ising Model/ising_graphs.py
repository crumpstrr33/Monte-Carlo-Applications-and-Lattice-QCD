'''
This program calculates and plots the evolution of the 2 dimensional Ising
model system for various values of Jkbt with respect to the measurement
of the magnetization.
'''
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import product

L = 50
## Crit Temp is at Jkbt = 0.5 * np.log(1 + np.sqrt(2))
mcsmax = 1000
magList = []
magTot = []
JkbtVals = [0.25, 0.375, 0.5, 0.75, 1]
lattice = np.ones(L**2).reshape((L,L))
fig = plt.figure()
ax1 = fig.add_subplot(111)

def Ising2D(energies, lattice, Jkbt):
    ## Finds the magnitization of the lattice after each 10*L potential flips
    magnetization = 0
    for i in range(L):
        for j in range(L):
           magnetization += lattice[i][j]

    for p in product(range(L), range(L)):
        i, j = p
        energies[i][j] = (lattice[i][(j - 1) % L] + lattice[i][(j + 1) % L] + \
                          lattice[(i - 1) % L][j] + lattice[(i + 1) % L][j]) * \
                          lattice[i][j]
        if energies[i][j] <= 0 or np.random.random() < np.exp(-2 * Jkbt * energies[i][j]):
            lattice[i][j] *= -1
                    
    return lattice, magnetization

for i in range(len(JkbtVals)):
    magList = []
    energies = np.zeros(L**2).reshape((L,L))

    ## Creates LxL array of random distribution of spin up and spin down.
    for row in range(L):
        for col in range(L):
            lattice[row][col] = random.randint(0, 1)
            if lattice[row][col] == 0:
                lattice[row][col] = -1

    for _ in range(mcsmax):
        lattice, magnetization = Ising2D(energies, lattice, JkbtVals[i])
        magList.append(magnetization/L**2)

    magTot.append(magList)

for i in range(len(magTot)):
    ax1.plot(magTot[i], lw = '2', label = JkbtVals[i])

ax1.legend()
plt.show()