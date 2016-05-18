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

N = [5, 10, 20, 30, 50]  ## Lattice Sizes
TC = 0.5 * np.log(1 + np.sqrt(2))  ## Critical Temperature
SWEEPS = 10000  ## Number of Sweeps


'''
Creates the list used to plot the data.
'''
def fList_define():
    f = 0.6
    f_list = []

    while f < 1.41:
        f_list.append(f)
        f += 0.05

    return f_list


'''
Creates an NxN lattice where each site is randomnly give a state of spin up (+1)
or a state of spin down (-1).
'''
def build_random_lat(N):
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
def ising_2d(lattice, jkbt, N):
    magnetization = 0

    for p in product(range(N), range(N)):
        i, j = p
        ## Calculates energy change at site (i,j)
        energy = (lattice[i][(j - 1) % N] + lattice[i][(j + 1) % N] + \
                          lattice[(i - 1) % N][j] + lattice[(i + 1) % N][j]) * \
                          lattice[i][j]
        magnetization += lattice[i][j]

        ## Metropolis update, determines if site should flip
        if energy <= 0 or np.random.random() < np.exp(-2 * jkbt * energy):
            lattice[i][j] *= -1

    return lattice, magnetization


'''
For a lattice size L, runs SWEEPS sweeps of the lattice, recording the 
magnetization, then increases the T/TC value from 0.6-1.4 by increments of
0.05.
'''
def mag_vs_t(N):
    factor = 0.6
    lattice = np.ones(N**2).reshape((N,N))
    avg_mag_temp = []
    avg_sus_temp = []

    while factor < 1.41:
        mag_tot_sweeps = []
        mag_tot_sq_sweeps = []

        lattice = build_random_lat(N)

        ## Stores the absolute magnetization from each sweep
        for i in range(SWEEPS):
            lattice, magnetization = ising_2d(lattice, TC / factor, N)
            if i > 3 * SWEEPS / 4:
                mag_tot_sweeps.append(abs(magnetization))
                mag_tot_sq_sweeps.append(magnetization**2)

        avgMag = np.average(mag_tot_sweeps)
        ## Susceptibility calculated by variance of the magnetization
        avgSus = (np.average(mag_tot_sq_sweeps) - 
                  np.average(mag_tot_sweeps)**2) / (TC * factor)

        ## Normalizing the data
        avg_mag_temp.append(avgMag / N**2)
        avg_sus_temp.append(avgSus / N**2)

        ## Increase the T/Tc values by 0.05
        factor += 0.05

    return avg_mag_temp, avg_sus_temp


'''
Main function
'''
def main():
    avg_mag_tot = []
    avg_sus_tot = []

    ## Runs the functions for various lattice size in the list N
    for i in N:
        avgMagTemps, avgSusTemps = mag_vs_t(i)
        print('Done with N = %d' %i)
    
        ## All of the data is stored in these two lists
        avg_sus_tot.append(avgSusTemps)
        avg_mag_tot.append(avgMagTemps)
    
    ## Saves the data collected as .txt files.
    with open('avgMagTot.txt', 'w') as f:
        for i in range(len(avg_mag_tot)):
            f.write('N = %d:\n' % N[i])
            f.write(', '.join(str(x) for x in avg_mag_tot[i]))
            f.write('\n')
    with open('avgSusTot.txt', 'w') as f:
        for i in range(len(avg_sus_tot)):
            f.write('N = %d:\n' % N[i])
            f.write(", ".join(str(x) for x in avg_sus_tot[i]))
            f.write('\n')

    '''
    Plots the absolute magnetization and susceptibility versus temperature
    (T/Tc), respectively. Each plot contains one line for each lattice size
    and the plots are of temperature vs. the magnetization and the
    susceptibility at equilibrium for the respective size/temperatures.
    '''   
    fList = fList_define()    
    
    plt.figure()
    for i in range(len(avg_mag_tot)):
        plt.plot(fList, avg_mag_tot[i], label = r'$N=%d$' % N[i], lw = '2')
    plt.ylim(0, 1)
    plt.xlim(0.6, 1.4)
    plt.title('Magnetization for $N$x$N$ Lattices at Different Temperatures' % N)
    plt.legend()
    plt.xlabel(r'$T/T_c$')
    plt.ylabel('Absolute Magnetization per Site')
    plt.show()

    plt.figure()
    for i in range(len(avg_sus_tot)):
        plt.plot(fList, avg_sus_tot[i], label = r'$N=%d$' % N[i], lw = '2')
    plt.xlim(0.6, 1.4)
    plt.title('Susceptibility for $N$x$N$ Lattices at Different Temperatures' % N)
    plt.legend()
    plt.xlabel(r'$T/T_c$')
    plt.ylabel('Susceptilibity per Site')
    plt.show()


if __name__ == "__main__":
    main()