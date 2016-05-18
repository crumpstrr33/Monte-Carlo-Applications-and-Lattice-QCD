'''
This program models the Ising model. It runs through mcsmax frames of the
model for a NxN lattice with a parameter of Jkbt. Once done, it uses animate
from matplotlib to run through every frame as an animation visually
showing the evolution of the Ising model. The ith frame is saved in the ith
element of i_plot.
'''
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from itertools import product

N = 20 ## Lattice size
JKBT = 0.7 ## Temperature
SWEEPS = 100 ## Total number of sweeps


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
def ising_2d(lattice):
    plot = []

    for p in product(range(N), range(N)):
        i, j = p
        ## Calculates energy change at site (i,j)
        energy = (lattice[i][(j - 1) % N] + lattice[i][(j + 1) % N] + \
                  lattice[(i - 1) % N][j] + lattice[(i + 1) % N][j]) * \
                  lattice[i][j]

        ## Metropolis update, determines if site should flip
        if energy <= 0 or np.random.random() < np.exp(-2 * JKBT * energy):
            lattice[i][j] *= -1

        plot.append(lattice[i][j])

    return plot


'''
Animates each sweep of the lattice with each sweep being one frame on
the plot.
'''  
def ising_plots(intr, i_plot):
    plt.title(intr + 1)
    plt.imshow(i_plot[intr], cmap = cm.Greys, extent = (0, N, 0, N),
               interpolation = 'nearest')


'''
Main function
'''
def main():
    lattice = build_random_lat()

    ## Each frame of the Ising model is stored in an element of i_plot
    i_plot = [0] * SWEEPS
    for i in range(SWEEPS):
        i_plot[i] = ising_2d(lattice)
    i_plot = np.array(i_plot).reshape(SWEEPS, N, N)

    ## Animates the data from i_plot
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.clear()
    anim = animation.FuncAnimation(fig, ising_plots, frames = SWEEPS,
                                   fargs = [i_plot], repeat=False, interval=5)
    
    ## Saves the animation as an .mp4 in the below path.
    plt.rcParams['animation.ffmpeg_path'] = \
                'C:\\Program Files (x86)\\Anaconda3\\Library\\bin\\ffmpeg.exe'
    mywriter = animation.FFMpegWriter(fps=20)
    anim.save('Ising Model (%dx%d), %d frames.mp4' % (N, N, SWEEPS), 
              writer = mywriter)


if __name__ == "__main__":
    main()