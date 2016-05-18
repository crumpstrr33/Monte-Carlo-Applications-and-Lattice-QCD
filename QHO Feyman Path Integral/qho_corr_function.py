'''
Calculates the correlation between slices of a path for a specific potential
and correlation function (G). Also, the excitation energy for the quantum
harmonic oscillator is calculated from E0 to E1 due to assumptions that can be 
made on the correlation function for the QHO. N_CF are generated and binned
into BIN_NUM paths to save on computing time.
'''
import numpy as np
import matplotlib.pyplot as plt

N = 20  ## Number of slices
N_COR = 20  ## Paths thrown away between measurements
N_CF = 100  ## Number of paths
A = 0.5  ## Spacing of slices
EPS = 1.4 ## Upper bound for change in path site
BIN_NUM = 10  ## Number of bins
BS_NUM = 100  ## Number of bootstraps
LEN_G = int(N_CF / BIN_NUM)


'''
Calculate the change in energy for some change in the path.
'''
def dS(x, j, alpha):
    return (alpha / A) * ((A**2 + 2)*x[j] + (A**2/2 + 1)*alpha - 
                          x[(j + 1) % N] - x[(j - 1) % N])


'''
Updates the path by adding a small amount (EPS) to one point on the path and
checking if it should stay via Metropolis algorithm for every point.
'''
def update(x):
    for j in range(N):
        alpha = np.random.uniform(-EPS, EPS)
        delta_S = dS(x, j, alpha)        
        if delta_S < 0 or np.random.random() < np.exp(-delta_S):
            x[j] += alpha


'''
Computes the correlation function for some path configuration by multiyplying
together path points since the weighted average can be approximated by
the unweighted average.
'''
def compute_G(x, k):
    g = 0

    for j in range(N):
        g += x[j] * x[(j + k) % N]

    return g / N


'''
Generates N_cf paths at equilibrium that are separated by N_cor Metropolis
updates to get rid of any correlation.
'''
def mc_paths(G, x):
    for i in range(N_CF):
        ## Updates the path to remove any correlation of paths.
        for j in range(N_COR):
            update(x)
        
        ## For the ith configuration, the values for the correlation function are
        ## stored in the ith element of list G.
        for k in range(N):
            G[i][k] = compute_G(x, k)


'''
Averages every binNum paths together.
'''
def binning(G):
    bin_G = []

    for i in range(0, len(G), BIN_NUM):
        sum_bin_G = 0
        for j in range(BIN_NUM):
            sum_bin_G += G[i + j]
        bin_G.append(sum_bin_G / BIN_NUM)

    return bin_G


'''
Creates a list of randomnly chosen paths from the original paths of the same
length as the original list.
'''
def boot_strap(G):
    bs_list = []

    for i in range(LEN_G):
        bs_list.append(G[np.random.randint(0, LEN_G)])

    return bs_list


'''
Calculates the average for each N segments of every lenG paths.
'''
def avg_G(G):
    avg_G_list = []

    for i in range(N):
        avg_G = 0
        for j in range(LEN_G):
            avg_G += G[j][i]
        avg_G /= LEN_G
        avg_G_list.append(avg_G)

    return avg_G_list


'''
Calulcates the data for the excitation energy, correlation function and
their respective errors
'''
def calc_data_error(avg_tot_G, bs_paths):
    error1 = []
    error2 = []
    delta_E = []

    for i in range(N):
        path_part_1 = []
        path_part_2 = []
        delta_E.append(np.log(np.abs(avg_tot_G[i] / 
                       avg_tot_G[(i + 1) % N])) / A)
    
        for j in range(BS_NUM):
            path_part_1.append(np.log(np.abs(bs_paths[j][i] / 
                             bs_paths[j][(i + 1) % N])) / A)
            path_part_2.append(bs_paths[j][i])
    
        error1.append(np.std(path_part_1))
        error2.append(np.std(path_part_2))

    return error1, error2, delta_E


'''
Main function
'''
def main():
    x = np.zeros(N)
    G_corr = np.zeros((N_CF, N))
    avg_tot_G = []

    bs_paths = []

    ## Thermalizes the lattice
    for _ in range(5 * N_COR):
        update(x)
    
    ## Calculates the correlation function for N_cf paths and bins the results
    mc_paths(G_corr, x)
    bin_G = binning(G_corr)
    
    ## Creates a list of bsNum bootstraps for error calculation
    for _ in range(BS_NUM):
        bs_paths.append(avg_G(boot_strap(bin_G)))
    
    ## Averages all paths of G for the binned results
    avg_tot_G = avg_G(bin_G)

    ## Calculates the exictation energy data and the error sets
    excite_error, corr_error, delta_E = calc_data_error(avg_tot_G, bs_paths)
    
    ## Plots the data.
    x_range = []
    for i in range(N):
       x_range.append(i / 2 + 0.5)  

    plt.figure()
    plt.ylim(0.5, 1.5)
    plt.xlim(0, 3.3)
    plt.ylabel('$E_1-E_0$')
    plt.xlabel('t')
    plt.errorbar(x_range, delta_E, yerr = excite_error, fmt = 'o')
    plt.plot(x_range, delta_E, 'o')
    plt.plot((-0.1, N), (1, 1), color = 'b')
    
    plt.figure()
    plt.ylim(0, 0.5)
    plt.xlim(0, N / 2)
    plt.ylabel('$G(t)$')
    plt.xlabel('t')
    plt.errorbar(x_range, avg_tot_G, yerr = corr_error)
    plt.plot(x_range, avg_tot_G, 'o')


if __name__ == "__main__":
    main()