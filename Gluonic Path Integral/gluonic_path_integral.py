'''
Three different sized Wilson loops (AxA, Ax2A and Ax3A) are calculated on a 
spacetime lattice using an improved lattice involving rectangle operators
for a second-order correction. The program can also be run for the unimproved 
action with the boolean IMP_ACT.
'''
import numpy as np
import cmath
import datetime as dt
from itertools import product
import matrix_functions

IMP_ACT = False
N = 4 ## Lattice size
N_COR = 1 ## Number of sweeps before measurement
N_CF = 25 ## Number of measurements
N_MA = 10 ## Metropolis steps per link
EPS = 0.24
NUM_MAT = 100 ## Number of SU3 matrices and their inverses used
if IMP_ACT: ## Extra constants for the improved action
    BETA = 1.719
    MU0 = 0.797
    C1 = 5 / (3 * MU0**4)
    C2 = 1 / (12 * MU0**6)
else:
    BETA = 5.5

mat3 = matrix_functions.Matrix_NxN(3, EPS)

j = cmath.sqrt(-1) ## Imaginary unit
I = mat3.I() ## Identity matrix 
SU3_LIST = []
for i in range(NUM_MAT): ## Generate SU3 matrices and their inverses
    SU3_LIST.append(mat3.SUn())
    SU3_LIST.append(np.linalg.inv(SU3_LIST[2 * i]))


'''
Returns modulus of x of n (used for periodic boundary conditions)
'''
def m(x, n=N):
    return x % n


'''
Returns a matrix's conjugate transpose (i.e. inverse of an SU3 matrix).
'''
def ct(matrix):
    return matrix.conj().T


'''
Returns the matrix product of mats.
'''
def mm(*mats):
    product = I
    for i in mats:
        product = np.dot(product, i)

    return product


'''
Calculates the staples for the rectangle operators R_mn and R_nm for the link at
(x, y, z, t, n). This returns the matrix total for the link and the conjugate
of the link.
'''
def n_rect(L, x, y, z, t, n):
    a = [0] * 10
    a[n] = 1
    r_mn = 0
    rdmn = 0
    r_nm = 0
    rdnm = 0

    for i in range(3 - n):
        r_mn += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n],                                       \
                   ct(L[m(x+2*a[0])][m(y+2*a[1]-a[-i])][m(z+2*a[2]-a[1-i])][m(t-a[2-i])][i+n+1]),  \
                   ct(L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][n]),            \
                   ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n]),                              \
                   L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][i+n+1]) +                            \
                mm(ct(L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][i+n+1]),        \
                   ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n]),                              \
                   ct(L[m(x-a[0])][m(y-a[1]-a[-i])][m(z-a[2]-a[1-i])][m(t-a[2-i])][n]),            \
                   L[m(x-a[0])][m(y-a[1]-a[-i])][m(z-a[2]-a[1-i])][m(t-a[2-i])][i+n+1],            \
                   L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][n])
        rdmn += mm(L[x][y][z][t][i+n+1],                                                           \
                   L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n],                                  \
                   L[m(x+a[0])][m(y+a[1]+a[-i])][m(z+a[2]+a[1-i])][m(t+a[2-i])][n],                \
                   ct(L[m(x+2*a[0])][m(y+2*a[1])][m(z+2*a[2])][t][i+n+1]),                         \
                   ct(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n])) +                                 \
                mm(ct(L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][n]),                                   \
                   L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][i+n+1],                                   \
                   L[m(x-a[0])][m(y-a[1]+a[-i])][m(z-a[2]+a[1-i])][m(t+a[2-i])][n],                \
                   L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n],                                  \
                   ct(L[m(x+a[0])][m(y+-a[1])][m(z+a[2])][t][i+n+1]))
        r_nm += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][i+n+1],                                   \
                   L[m(x+a[0])][m(y+a[1]+a[-i])][m(z+a[2]+a[1-i])][m(t+a[2-i])][i+n+1],            \
                   ct(L[x][m(y+2*a[-i])][m(z+2*a[1-i])][m(t+2*a[2-i])][n]),                        \
                   ct(L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][i+n+1]),                          \
                   ct(L[x][y][z][t][i+n+1]))
        rdnm += mm(ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][i+n+1]),                          \
                   ct(L[x][m(y-2*a[-i])][m(z-2*a[1-i])][m(t-2*a[2-i])][i+n+1]),                    \
                   L[x][m(y-2*a[-i])][m(z-2*a[1-i])][m(t-2*a[2-i])][n],                            \
                   L[m(x+a[0])][m(y+a[1]-2*a[-i])][m(z+a[2]-2*a[1-i])][m(t-2*a[2-i])][i+n+1],      \
                   L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][i+n+1])
    for i in range(n):
        r_mn += mm(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i],                                       \
                   L[m(x+a[i+n])][m(y+a[i+n-1]+a[1])][m(z+a[i+n-2]+a[2])][m(t+a[3])][i],           \
                   ct(L[m(x+2*a[i+n])][m(y+2*a[i+n-1])][m(z+2*a[i+n-2])][t][n]),                   \
                   ct(L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][i]),                         \
                   ct(L[x][y][z][t][i]))
        rdmn += mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i]),                         \
                   ct(L[m(x-2*a[i+n])][m(y-2*a[i+n-1])][m(z-2*a[i+n-2])][t][i]),                   \
                   L[m(x-2*a[i+n])][m(y-2*a[i+n-1])][m(z-2*a[i+n-2])][t][n],                       \
                   L[m(x-2*a[i+n])][m(y-2*a[i+n-1]+a[1])][m(z-2*a[i+n-2]+a[2])][m(t+a[3])][i],     \
                   L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i])
        r_nm += mm(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][n],                                       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]+2*a[1])][m(z-a[i+n-2]+2*a[2])][m(t+2*a[3])][i]), \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][n]),       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                         \
                   L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i]) +                           \
                mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i]),       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                         \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]-a[1])][m(z-a[i+n-2]-a[2])][m(t-a[3])][n]),       \
                   L[m(x-a[i+n])][m(y-a[i+n-1]-a[1])][m(z-a[i+n-2]-a[2])][m(t-a[3])][i],           \
                   L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][n])
        rdnm += mm(L[x][y][z][t][i],                                                               \
                   L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                             \
                   L[m(x+a[i+n])][m(y+a[i+n-1]+a[1])][m(z+a[i+n-2]+a[2])][m(t+a[3])][n],           \
                   ct(L[x][m(y+2*a[1])][m(z+2*a[2])][m(t+2*a[3])][i]),                             \
                   ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][n])) +                                 \
                mm(ct(L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][n]),                                   \
                   L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][i],                                       \
                   L[m(x+a[i+n])][m(y+a[i+n-1]-a[1])][m(z+a[i+n-2]-a[2])][m(t-a[3])][n],           \
                   L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                             \
                   ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i]))    

    return r_mn + r_nm, rdmn + rdnm


'''
Calculates the staples for the plaquette for the link at
(x, y, z, t, n). This returns the matrix total for the link and the conjugate
of the link.
'''
def n_plaq(L, x, y, z, t, n):
    a = [0] * 10
    a[n] = 1
    staple = 0
    staple_conj = 0

    for i in range(3 - n):
        staple += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n+i+1],                             \
                  ct(L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n]),                           \
                  ct(L[x][y][z][t][n+i+1]))
        staple_conj += mm(ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n+i+1]),               \
                       L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n],                          \
                       L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][t-a[2-i]][n+i+1])
    for i in range(n):
        staple += mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i]), \
                     ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                   \
                     L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i])
        staple_conj += mm(L[x][y][z][t][i],                                                    \
                       L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                     \
                       ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i]))

    return staple, staple_conj


'''
Calculates the change in the action for the link at (x, y, z, t, n). If 
IMP_ACT = True, then the improved action is calcualted, otherwise the 
unimproved action is calculated. The change in action and the random SU3 matrix
used is returned.
'''
def dS(L, staples, stapleC, rects, rectC, x, y, z, t, n):
    M = SU3_LIST[np.random.randint(0, 2 * NUM_MAT)]
    U = L[x][y][z][t][n]

    ## Calculates change in the action based on which action is used
    if IMP_ACT:
        delta_S = (BETA / 3) * np.real(np.trace(np.dot(U - mm(M, U), C1 * staples - C2 * rects) + \
                                            np.dot(ct(U) - ct(mm(M, U)), C1 * stapleC - C2 * rectC)))
    else:
        delta_S = (BETA / 3) * np.real(np.trace(mm(U - mm(M, U), staples) + \
                                               mm(ct(U) - ct(mm(M, U)), stapleC)))

    return delta_S, M


'''
Completes one sweep of the lattice. The function n_rect is run only if the
improved action is being used, otherwise the variables rect and rect_conj
are set to 0.
'''
def sweep(L):
    for p in product(range(N), range(N), range(N), range(N), range(4)):
        x, y, z, t, n = p
        rect = 0
        rect_conj = 0

        ## Calculates the respective staples
        staple, staple_conj = n_plaq(L, x, y, z, t, n)
        if IMP_ACT:
            rect, rect_conj = n_rect(L, x, y, z, t, n)

        for _ in range(N_MA):
            delta_S, M = dS(L, staple, staple_conj, rect, rect_conj, x, y, z, t, n)
            if delta_S < 0 or np.random.random() < np.exp(-delta_S):
                L[x][y][z][t][n] = np.dot(M, L[x][y][z][t][n])


'''
Function that the function measure uses since the Wilson loops is longer for 
Ax2A and Ax3A loops. It returns the extra part of the loop.
'''
def measureLink(L, x, y, z, t, loop_length):
    if loop_length == 3:
        link = [mm(L[m(x+1)][m(y+1)][z][t][1], L[m(x+1)][m(y+2)][z][t][1], ct(L[x][m(y+3)][z][t][0]), ct(L[x][m(y+2)][z][t][1]), ct(L[x][m(y+1)][z][t][1])), \
                mm(L[m(x+1)][y][m(z+1)][t][2], L[m(x+1)][y][m(z+2)][t][2], ct(L[x][y][m(z+3)][t][0]), ct(L[x][y][m(z+2)][t][2]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[m(x+1)][y][z][m(t+1)][3], L[m(x+1)][y][z][m(t+2)][3], ct(L[x][y][z][m(t+3)][0]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][m(y+1)][m(z+1)][t][2], L[x][m(y+1)][m(z+2)][t][2], ct(L[x][y][m(z+3)][t][1]), ct(L[x][y][m(z+2)][t][2]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[x][m(y+1)][z][m(t+1)][3], L[x][m(y+1)][z][m(t+2)][3], ct(L[x][y][z][m(t+3)][1]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][y][m(z+1)][m(t+1)][3], L[x][y][m(z+1)][m(t+2)][3], ct(L[x][y][z][m(t+3)][2]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3]))]   
    elif loop_length == 2:
        link = [mm(L[m(x+1)][m(y+1)][z][t][1], ct(L[x][m(y+2)][z][t][0]), ct(L[x][m(y+1)][z][t][1])), \
                mm(L[m(x+1)][y][m(z+1)][t][2], ct(L[x][y][m(z+2)][t][0]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[m(x+1)][y][z][m(t+1)][3], ct(L[x][y][z][m(t+2)][0]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][m(y+1)][m(z+1)][t][2], ct(L[x][y][m(z+2)][t][1]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[x][m(y+1)][z][m(t+1)][3], ct(L[x][y][z][m(t+2)][1]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][y][m(z+1)][m(t+1)][3], ct(L[x][y][z][m(t+2)][2]), ct(L[x][y][z][m(t+1)][3]))]
    else:
        link = [ct(L[x][m(y+1)][z][t][0]), ct(L[x][y][m(z+1)][t][0]), ct(L[x][y][z][m(t+1)][0]), \
                ct(L[x][y][m(z+1)][t][1]), ct(L[x][y][z][m(t+1)][1]), ct(L[x][y][z][m(t+1)][2])]

    return link


'''
Measures Wilson loops and normalizes and averages them over the lattice.
'''
def measure(L, loop_length):
    LAvg = 0

    for p in product(range(N), range(N), range(N), range(N)):
        x, y, z, t = p
        link = measureLink(L, x, y, z, t, loop_length)
        LAvg += np.real(np.trace( \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][1], link[0], ct(L[x][y][z][t][1])) + \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][2], link[1], ct(L[x][y][z][t][2])) + \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][3], link[2], ct(L[x][y][z][t][3])) + \
        mm(L[x][y][z][t][1], L[x][m(y+1)][z][t][2], link[3], ct(L[x][y][z][t][2])) + \
        mm(L[x][y][z][t][1], L[x][m(y+1)][z][t][3], link[4], ct(L[x][y][z][t][3])) + \
        mm(L[x][y][z][t][2], L[x][y][m(z+1)][t][3], link[5], ct(L[x][y][z][t][3]))))

    return LAvg / N**4 / 18


'''
Main function
'''
def main():
    lat_tots = [[], [], []]
    L = np.array([I] * N**4 * 4).reshape((N,N,N,N,4,3,3))    

    t0 = dt.datetime.now()

    ## Thermalize the lattice
    for i in range(100):
        sweep(L)

    t1 = dt.datetime.now()

    for i in range(N_CF):
        for _ in range(N_COR):
            sweep(L)

        ## Records the measurement for each Wilson loop
        for j in range(3):
            lat_tots[j].append(measure(L, j + 1))

        ## Print out the data for the current sweep
        space = ' ' * (5 - len(str(i + 1)))
        print('Config %d:%sAxA  - %.5f' % ((i + 1), space, lat_tots[0][i]))
        print('             Ax2A - %.5f' % (lat_tots[1][i]))
        print('             Ax3A - %.5f' % (lat_tots[2][i]))

    t2 = dt.datetime.now()

    ## Print out the time taken for various parts of the program.
    print('\nTime total:             %.3f seconds' %  (t2 - t0).total_seconds())
    print('Time per configuration: %.3f seconds' % ((t2 - t1).total_seconds() / N_CF))
    print('Time per sweep:         %.3f seconds\n' % ((t2 - t1).total_seconds() / N_CF / N_COR))
   
    ## Prints out average and standard deviation of the loops
    for i in range(3):
        print('Loop AxA:   Average - %.8f' % np.average(lat_tots[i]))
        print('                Std - %.8f\n' % np.std(lat_tots[i]))

    ## Saves relevant data to be analyzed inlqcd_analysis.py
    for i in range(1, 4):
        with open('NoImp_Ax%sA_Measurements.txt' % i, 'w') as f:
            f.write(', '.join(str(x) for x in lat_tots[i - 1]))
            f.write('\n\nAverage: %.8f' % np.average(lat_tots[i - 1]))
            f.write('\n\nStd: %.8f' % np.std(lat_tots[i - 1]))

    act = ''
    if not IMP_ACT:
        act = 'un'
    with open('parameters_%simproved.txt' % act, 'w') as f:
        f.write('%simproved action parameters:\n' % act)        
        f.write('N = %d\n' % N)
        f.write('N_cor = %d\n' % N_COR)
        f.write('N_cf = %d\n' % N_CF)
        f.write('N_ma = %d\n' % N_MA)
        f.write('eps = %.2f\n' % EPS)
        f.write('beta = %.3f\n' % BETA)
        if IMP_ACT:
            f.write('mu0 = %.3f\n' % MU0)


if __name__ == "__main__":
    main()