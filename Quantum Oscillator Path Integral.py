'''
Calculates the correlation between slices of a path for a specific potential
and correlation function (G). Also, the excitation energy for the quantum harmonic
oscillator is calculated from E0 to E1 due to assumptions that can be made on
the correlation function for the QHO.
'''
import numpy as np
import matplotlib.pyplot as plt

## Set parameters:
N = 20                      ## Number of slices
N_cor = 20                  ## Paths thrown away between measurements
N_cf = 10000                ## Number of paths
a = 0.5                     ## Spacing of slices
eps = 1.4
binNum = 100                ## Number of bins
bsNum = 100                 ## Number of bootstraps
lenG = int(N_cf / binNum)

## Create arrays:
x = np.zeros(N)
GCorr = np.zeros((N_cf, N))
avgG = []
dE = []
bootE = []
bootError = []
error1 = []
error2 = []

## Calculate the change in energy for some change in the path.
def dS(j, alpha):
    return (alpha / a) * ((a**2 + 2)*x[j] + (a**2/2 + 1)*alpha - x[(j + 1) % N] - x[(j - 1) % N])

## Updates the path by adding a small amount (eps) to one point on the path and
## checking if it should stay via Metropolis algorithm for every point.
def update(x):
    for j in range(N):
        alpha = np.random.uniform(-eps, eps)
        deltaS = dS(j, alpha)        
        if deltaS < 0 or np.random.random() < np.exp(-deltaS):    
            x[j] += alpha

## Computes the correlation function for some path configuration by multiyplying
## together path points since the weighted average can be approximated by
## the unweighted average.
def computeG(x, k):
    g = 0
    for j in range(N):
        g += x[j] * x[(j + k) % N]
    return g / N

## Generates N_cf paths at equilibrium that are separated by N_cor Metropolis
## updates to get rid of any correlation.
def MCPaths(G):
    for i in range(N_cf):
        ## Updates the path to remove any correlation of paths.
        for j in range(N_cor):
            update(x)
        
        ## For the ith configuration, the values for the correlation function are
        ## stored in the ith element of list G.
        for k in range(N):
            G[i][k] = computeG(x, k)

## Creates a list of randomnly chosen paths from the original paths.
def bootstrap(G):
    BSList = []
    for i in range(lenG):
        BSList.append(G[np.random.randint(0, lenG)])
    return BSList

## Averages every binNum paths together.
def binning(G):
    binG = []
    for i in range(0, len(G), binNum):
        avgG = 0
        for j in range(binNum):
            avgG += G[i + j]
        binG.append(avgG / binNum)
    return binG

## Calculates the average for each N segments of every lenG paths.
def avgG(G):
    avgGList = []
    for i in range(N):
        GAvg = 0
        for j in range(lenG):
            GAvg += G[j][i]
        GAvg /= lenG
        avgGList.append(GAvg)
    return avgGList

## Thermalizes the lattice.
for _ in range(5 * N_cor):
    update(x)

## Calculates the correlation function for N_cf paths and bins the results.
MCPaths(GCorr)
binG = binning(GCorr)

## Creates a list of bsNum bootstraps for error calculation.
for i in range(bsNum):
    bootE.append(avgG(bootstrap(binG)))

## Averages all paths of G for the binned results.
avgG = avgG(binG)

## Calulcates the data for the excitation energy, correlation function and
## their respective errors.
for i in range(N):
    pathPart1 = []
    pathPart2 = []
    dE.append(np.log(np.abs(avgG[i] / avgG[(i + 1) % N])) / a)
    for j in range(bsNum):
        pathPart1.append(np.log(np.abs(bootE[j][i] / bootE[j][(i + 1) % N])) / a)
        pathPart2.append(bootE[j][i])
    error1.append(np.std(pathPart1))
    error2.append(np.std(pathPart2))


## Plots the data.
xRange = []
for i in range(N):
   xRange.append(i / 2 + 0.5)  

plt.figure()
plt.ylim(0.5, 1.5)
plt.xlim(0, 3.3)
plt.ylabel('$E_1-E_0$')
plt.xlabel('t')
plt.errorbar(xRange, dE, yerr = error1, fmt = 'o')
plt.plot(xRange, dE, 'o')
plt.plot((-0.1, N), (1, 1), color = 'b')

plt.figure()
plt.ylim(0, 0.5)
plt.xlim(0, N / 2)
plt.ylabel('$G(t)$')
plt.xlabel('t')
plt.errorbar(xRange, avgG, yerr = error2)
plt.plot(xRange, avgG, 'o')