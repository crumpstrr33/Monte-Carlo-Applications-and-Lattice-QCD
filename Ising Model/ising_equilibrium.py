'''
This program models the Ising model. It runs through mcsmax frames of the
model for a NxN lattice with a parameter of Jkbt. Once done, it uses animate
from matplotlib to run through every frame as an animation visually
showing the evolution of the Ising model. The ith frame is save in the ith
element of iPlot.
'''
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from itertools import product

N = 1000
Jkbt = 0.5 * np.log(1 + np.sqrt(2))
mcsmax = 1000
magList = []
lattice = np.zeros(N**2).reshape((N,N))
fig = plt.figure()
ax1 = fig.add_subplot(111)

## Creates LxL array of random distribution of spin up and spin down.
for i, j in product(range(N), range(N)):
    lattice[i][j] = random.randint(0, 1)
    if lattice[i][j] == 0:
        lattice[i][j] = -1

## Completes one sweep of the lattice, flipping the spin whenever the if
## statement is true.
def Ising2D():
    plot = []

    for p in product(range(N), range(N)):
        i, j = p
        energy = (lattice[i][(j - 1) % N] + lattice[i][(j + 1) % N] + \
                          lattice[(i - 1) % N][j] + lattice[(i + 1) % N][j]) * \
                          lattice[i][j]
        if energy <= 0 or np.random.random() < np.exp(-2 * Jkbt * energy):
            lattice[i][j] *= -1
            
        plot.append(lattice[i][j])

    return plot

## Animates each sweep of the lattice with each sweep being one frame on
## the plot.    
def IsingPlots(intr):
    ax1.clear()
    plt.title(intr + 1)
    ax1.imshow(iPlot[intr], cmap = cm.Greys, extent = (0, N, 0, N),
               interpolation = 'nearest')

## Each element of iPlot contains one frame of the simulation and it is run
## with FuncAnimation.
iPlot = [0] * mcsmax
for i in range(mcsmax):
    iPlot[i] = Ising2D()
iPlot = np.array(iPlot).reshape(mcsmax, N, N)

anim = animation.FuncAnimation(fig, IsingPlots, frames = mcsmax, 
                          repeat = False, interval = 5)

## Saves the animation as ising.mp4 in the below path.
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Program Files (x86)\\Anaconda3\\Library\\bin\\ffmpeg.exe'
mywriter = animation.FFMpegWriter(fps=20)
anim.save('ising.mp4', writer = mywriter)