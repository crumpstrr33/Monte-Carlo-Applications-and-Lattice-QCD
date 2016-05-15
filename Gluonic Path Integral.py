import numpy as np
import cmath
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import product
import Matrix3x3Stuff

N = 8
N_cor = 1
N_cf = 100
N_ma = 10
eps = 0.24
beta = 1.719
#beta = 5.5
mu0 = 0.797
c1 = 5 / (3 * mu0**4)
c2 = 1 / (12 * mu0**6)
numMatrices = 100


j = cmath.sqrt(-1)
Mat3 = Matrix3x3Stuff.Matrix_nxn(3, eps)
I = Mat3.I()

latTots0 = []
latTots1 = []
latTots2 =[]
SU3list = []
for i in range(numMatrices):
    SU3list.append(Mat3.SUn())
    SU3list.append(np.linalg.inv(SU3list[2 * i]))
L = np.array([I] * N**4 * 4).reshape((N,N,N,N,4,3,3))

def ct(matrix):
    return matrix.conj().T
    
def m(x, n = N):
    return x % n
    
def mm(*mats):
    product = I
    for i in mats:
        product = np.dot(product, i)
    return product

def nRect(x, y, z, t, n):
    a = [0] * 10
    a[n] = 1
    Rmns = 0
    RmnC = 0
    Rnms = 0
    RnmC = 0

    for i in range(3 - n):
        Rmns += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n],                                       \
                   ct(L[m(x+2*a[0])][m(y+2*a[1]-a[-i])][m(z+2*a[2]-a[1-i])][m(t-a[2-i])][i+n+1]),  \
                   ct(L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][n]),            \
                   ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n]),                              \
                   L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][i+n+1]) +                            \
                mm(ct(L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][i+n+1]),        \
                   ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n]),                              \
                   ct(L[m(x-a[0])][m(y-a[1]-a[-i])][m(z-a[2]-a[1-i])][m(t-a[2-i])][n]),            \
                   L[m(x-a[0])][m(y-a[1]-a[-i])][m(z-a[2]-a[1-i])][m(t-a[2-i])][i+n+1],            \
                   L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][n])
        RmnC += mm(L[x][y][z][t][i+n+1],                                                           \
                   L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n],                                  \
                   L[m(x+a[0])][m(y+a[1]+a[-i])][m(z+a[2]+a[1-i])][m(t+a[2-i])][n],                \
                   ct(L[m(x+2*a[0])][m(y+2*a[1])][m(z+2*a[2])][t][i+n+1]),                         \
                   ct(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n])) +                                 \
                mm(ct(L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][n]),                                   \
                   L[m(x-a[0])][m(y-a[1])][m(z-a[2])][t][i+n+1],                                   \
                   L[m(x-a[0])][m(y-a[1]+a[-i])][m(z-a[2]+a[1-i])][m(t+a[2-i])][n],                \
                   L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n],                                  \
                   ct(L[m(x+a[0])][m(y+-a[1])][m(z+a[2])][t][i+n+1]))
        Rnms += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][i+n+1],                                   \
                   L[m(x+a[0])][m(y+a[1]+a[-i])][m(z+a[2]+a[1-i])][m(t+a[2-i])][i+n+1],            \
                   ct(L[x][m(y+2*a[-i])][m(z+2*a[1-i])][m(t+2*a[2-i])][n]),                        \
                   ct(L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][i+n+1]),                          \
                   ct(L[x][y][z][t][i+n+1]))
        RnmC += mm(ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][i+n+1]),                          \
                   ct(L[x][m(y-2*a[-i])][m(z-2*a[1-i])][m(t-2*a[2-i])][i+n+1]),                    \
                   L[x][m(y-2*a[-i])][m(z-2*a[1-i])][m(t-2*a[2-i])][n],                            \
                   L[m(x+a[0])][m(y+a[1]-2*a[-i])][m(z+a[2]-2*a[1-i])][m(t-2*a[2-i])][i+n+1],      \
                   L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][m(t-a[2-i])][i+n+1])
    for i in range(n):
        Rmns += mm(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i],                                       \
                   L[m(x+a[i+n])][m(y+a[i+n-1]+a[1])][m(z+a[i+n-2]+a[2])][m(t+a[3])][i],           \
                   ct(L[m(x+2*a[i+n])][m(y+2*a[i+n-1])][m(z+2*a[i+n-2])][t][n]),                   \
                   ct(L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][i]),                         \
                   ct(L[x][y][z][t][i]))
        RmnC += mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i]),                         \
                   ct(L[m(x-2*a[i+n])][m(y-2*a[i+n-1])][m(z-2*a[i+n-2])][t][i]),                   \
                   L[m(x-2*a[i+n])][m(y-2*a[i+n-1])][m(z-2*a[i+n-2])][t][n],                       \
                   L[m(x-2*a[i+n])][m(y-2*a[i+n-1]+a[1])][m(z-2*a[i+n-2]+a[2])][m(t+a[3])][i],     \
                   L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i])
        Rnms += mm(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][n],                                       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]+2*a[1])][m(z-a[i+n-2]+2*a[2])][m(t+2*a[3])][i]), \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][n]),       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                         \
                   L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i]) +                           \
                mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i]),       \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                         \
                   ct(L[m(x-a[i+n])][m(y-a[i+n-1]-a[1])][m(z-a[i+n-2]-a[2])][m(t-a[3])][n]),       \
                   L[m(x-a[i+n])][m(y-a[i+n-1]-a[1])][m(z-a[i+n-2]-a[2])][m(t-a[3])][i],           \
                   L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][n])
        RnmC += mm(L[x][y][z][t][i],                                                               \
                   L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                             \
                   L[m(x+a[i+n])][m(y+a[i+n-1]+a[1])][m(z+a[i+n-2]+a[2])][m(t+a[3])][n],           \
                   ct(L[x][m(y+2*a[1])][m(z+2*a[2])][m(t+2*a[3])][i]),                             \
                   ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][n])) +                                 \
                mm(ct(L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][n]),                                   \
                   L[x][m(y-a[1])][m(z-a[2])][m(t-a[3])][i],                                       \
                   L[m(x+a[i+n])][m(y+a[i+n-1]-a[1])][m(z+a[i+n-2]-a[2])][m(t-a[3])][n],           \
                   L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                             \
                   ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i]))    
    return Rmns + Rnms, RmnC + RnmC

def nPlaq(x, y, z, t, n):
    a = [0] * 10
    a[n] = 1
    staples = 0
    stapleC = 0

    for i in range(3 - n):
        staples += mm(L[m(x+a[0])][m(y+a[1])][m(z+a[2])][t][n+i+1],                             \
                      ct(L[x][m(y+a[-i])][m(z+a[1-i])][m(t+a[2-i])][n]),                        \
                      ct(L[x][y][z][t][n+i+1]))
        stapleC += mm(ct(L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n+i+1]),                    \
                      L[x][m(y-a[-i])][m(z-a[1-i])][m(t-a[2-i])][n],                            \
                      L[m(x+a[0])][m(y+a[1]-a[-i])][m(z+a[2]-a[1-i])][t-a[2-i]][n+i+1])
    for i in range(n):
        staples += mm(ct(L[m(x-a[i+n])][m(y-a[i+n-1]+a[1])][m(z-a[i+n-2]+a[2])][m(t+a[3])][i]), \
                      ct(L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][n]),                   \
                      L[m(x-a[i+n])][m(y-a[i+n-1])][m(z-a[i+n-2])][t][i])
        stapleC += mm(L[x][y][z][t][i],                                                         \
                      L[m(x+a[i+n])][m(y+a[i+n-1])][m(z+a[i+n-2])][t][n],                       \
                      ct(L[x][m(y+a[1])][m(z+a[2])][m(t+a[3])][i]))
    return staples, stapleC

def dS(staples, stapleC, rects, rectC, x, y, z, t, n):
    M = SU3list[np.random.randint(0, 2 * numMatrices)]
    U = L[x][y][z][t][n]
    deltaS = (beta / 3) * np.real(np.trace(np.dot(U - mm(M, U), c1 * staples - c2 * rects) + \
                                            np.dot(ct(U) - ct(mm(M, U)), c1 * stapleC - c2 * rectC)))
    #deltaS = (beta / 3) * np.real(np.trace(mm(U - mm(M, U), staples) + \
    #                                       mm(ct(U) - ct(mm(M, U)), stapleC)))                
    return deltaS, M

def sweep():
    for p in product(range(N), range(N), range(N), range(N), range(4)):
        x, y, z, t, n = p
        staples, stapleC = nPlaq(x, y, z, t, n)
        rects, rectC = nRect(x,y,z,t,n)
        for _ in range(N_ma):
            deltaS, M = dS(staples, stapleC, rects, rectC, x, y, z, t, n)
            if deltaS < 0 or np.random.random() < np.exp(-deltaS):
                L[x][y][z][t][n] = np.dot(M, L[x][y][z][t][n])

def measureLink(x, y, z, t, loopLength):
    if loopLength == 3:
        link = [mm(L[m(x+1)][m(y+1)][z][t][1], L[m(x+1)][m(y+2)][z][t][1], ct(L[x][m(y+3)][z][t][0]), ct(L[x][m(y+2)][z][t][1]), ct(L[x][m(y+1)][z][t][1])), \
                mm(L[m(x+1)][y][m(z+1)][t][2], L[m(x+1)][y][m(z+2)][t][2], ct(L[x][y][m(z+3)][t][0]), ct(L[x][y][m(z+2)][t][2]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[m(x+1)][y][z][m(t+1)][3], L[m(x+1)][y][z][m(t+2)][3], ct(L[x][y][z][m(t+3)][0]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][m(y+1)][m(z+1)][t][2], L[x][m(y+1)][m(z+2)][t][2], ct(L[x][y][m(z+3)][t][1]), ct(L[x][y][m(z+2)][t][2]), ct(L[x][y][m(z+1)][t][2])), \
                mm(L[x][m(y+1)][z][m(t+1)][3], L[x][m(y+1)][z][m(t+2)][3], ct(L[x][y][z][m(t+3)][1]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3])), \
                mm(L[x][y][m(z+1)][m(t+1)][3], L[x][y][m(z+1)][m(t+2)][3], ct(L[x][y][z][m(t+3)][2]), ct(L[x][y][z][m(t+2)][3]), ct(L[x][y][z][m(t+1)][3]))]   
    elif loopLength == 2:
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

def measure(loopLength):
    LAvg = 0
    for p in product(range(N), range(N), range(N), range(N)):
        x, y, z, t = p
        link = measureLink(x, y, z, t, loopLength)
        LAvg += np.real(np.trace( \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][1], link[0], ct(L[x][y][z][t][1])) + \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][2], link[1], ct(L[x][y][z][t][2])) + \
        mm(L[x][y][z][t][0], L[m(x+1)][y][z][t][3], link[2], ct(L[x][y][z][t][3])) + \
        mm(L[x][y][z][t][1], L[x][m(y+1)][z][t][2], link[3], ct(L[x][y][z][t][2])) + \
        mm(L[x][y][z][t][1], L[x][m(y+1)][z][t][3], link[4], ct(L[x][y][z][t][3])) + \
        mm(L[x][y][z][t][2], L[x][y][m(z+1)][t][3], link[5], ct(L[x][y][z][t][3]))))
    return LAvg / N**4 / 18


t0 = dt.datetime.now()

for i in range(100):
    sweep()

t1 = dt.datetime.now()

for i in range(N_cf):
    for _ in range(N_cor):
        sweep()

    latTots0.append(measure(1))
    latTots1.append(measure(2))
    latTots2.append(measure(3))

    if i < 9:
        print('Config %d:     AxA - %.5f' % ((i+1), latTots0[i]))
    elif i < 99:
        print('Config %d:    AxA - %.5f' % ((i+1), latTots0[i]))
    else:
        print('Config %d:   AxA - %.5f' % ((i+1), latTots0[i]))
    print('             Ax2A - %.5f' % (latTots1[i]))
    print('             Ax3A - %.5f' % (latTots2[i]))

t2 = dt.datetime.now()

plt.figure()
plt.ylim(0, 1)
plt.xlim(0,25)
plt.ylabel('Monte Carlo Measurement')
plt.xlabel('Sweeps')
plt.plot(latTots0, 'b', lw = '1.5')
plt.plot(latTots1, 'g', lw = '1.5')
plt.plot(latTots2, 'r', lw = '1.5')

Ax1A = mpatches.Patch(color = 'blue', label = '$a$ x $a$')
Ax2A = mpatches.Patch(color = 'green', label = '$a$ x $2a$')
Ax3A = mpatches.Patch(color = 'red', label = '$a$ x $3a$')
plt.legend(handles = [Ax1A, Ax2A, Ax3A])


print('\nTime total:             %.3f seconds' %  (t2 - t0).total_seconds())
print('Time per configuration: %.3f seconds' % ((t2 - t1).total_seconds() / N_cf))
print('Time per sweep:         %.3f seconds\n' % ((t2 - t1).total_seconds() / N_cf / N_cor))
print('Loop AxA:   Average - %.8f' % np.average(latTots0))
print('                Std - %.8f\n' % np.std(latTots0))
print('Loop Ax2A:  Average - %.8f' % np.average(latTots1))
print('                Std - %.8f\n' % np.std(latTots1))
print('Loop Ax3A:  Average - %.8f' % np.average(latTots2))
print('                Std - %.8f' % np.std(latTots2))


with open('AxA_Measurements.txt', 'w') as f:
    f.write(', '.join(str(x) for x in latTots0))
    f.write('\n\nAverage: %.8f' % np.average(latTots0))
    f.write('\n\nStd: %.8f' % np.std(latTots0))
with open('Ax2A_Measurements.txt', 'w') as f:
    f.write(", ".join(str(x) for x in latTots1))
    f.write('\n\nAverage: %.8f' % np.average(latTots1))
    f.write('\n\nStd: %.8f' % np.std(latTots1))
with open('Ax3A_Measurements.txt', 'w') as f:
    f.write(", ".join(str(x) for x in latTots2))
    f.write('\n\nAverage: %.8f' % np.average(latTots2))
    f.write('\n\nStd: %.8f' % np.std(latTots2))