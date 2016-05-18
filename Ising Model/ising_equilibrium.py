'''
This program calculates and plots the evolution of the 2 dimensional Ising
model system for various values of Jkbt with respect to the measurement
of the magnetization.
'''
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import product

N = 50 ## Size of lattice
SWEEPS = 100 ## Total number of sweeps per lattice
JKBT_VALS = [0.25, 0.375, 0.5, 0.75, 1] ## Temperature values 


'''
Creates an NxN lattice where each site is randomnly give a state of spin up (+1)
or a state of spin down (-1).
'''
def build_random_lat():
    lattice = np.zeros(N**2).reshape((N,N))

    for i, j in product(range(N), range(N)):
        lattice[i][j] = random.randint(0, 1)
        if lattice[i][j] == 0:
            lattice[i][j] = -1
    return lattice


'''
Completes one sweep of the lattice, checking every lattice site whether or not
it should flip and acting accordingly. It also records the total magnetization
of the permutation that the sweep creates.
'''
def ising_2d(lattice, jkbt):
    magnetization = 0

    for p in product(range(N), range(N)):
        i, j = p
        ## Calculates energy change at site (i,j)
        energies = ((lattice[i][(j - 1) % N] + lattice[i][(j + 1) % N] +
                     lattice[(i - 1) % N][j] + lattice[(i + 1) % N][j]) *
                     lattice[i][j])
        magnetization += lattice[i][j]

        ## Metropolis update, determines if site should flip
        if energies <= 0 or np.random.random() < np.exp(-2 * jkbt * energies):
            lattice[i][j] *= -1

    return lattice, magnetization


'''
With a value for temperature of Jkbt, SWEEPS sweeps are completed on the
lattice and the magnetization per site for each sweep is recorded as an element
of magList.
'''
def sweep_run(jkbt):
    mag_list = []
    lat = build_random_lat(N)

    for _ in range(SWEEPS):
        lat, magnetization = ising_2d(lat, jkbt)
        mag_list.append(magnetization/N**2)
    
    return mag_list


'''
Main function
'''
def main():
    ## Calculates magnetization per sweep for each temperature value
    magTot = []
    for i in JKBT_VALS:
        magTot.append(sweep_run(i))

    ## Plots each element of magTot
    fig = plt.figure()
    ax1 = fig.add_subplot(111)   
    for i in range(len(magTot)):
        ax1.plot(magTot[i], lw = '2', label = JKBT_VALS[i])
    ax1.legend()
    plt.show()


if __name__ == "__main__":
    main()