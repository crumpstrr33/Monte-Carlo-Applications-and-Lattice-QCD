'''
In this file, the class Matrix_nxn allows for certain operations and matrix 
creation for a specific n. The class Matrix_Const contains some common matrix
constants that may be needed.
'''
import numpy as np
import cmath
import math

j = cmath.sqrt(-1)

## Run the main program?
run = False

## Various functions for matrix creation/manipulation for dxd matrices with an
## optional factor eps added as a multiplicative constant.
## (defaults to 1 if not wanted)
class Matrix_nxn(object):
    def __init__(self, d, eps = 1):
        self.eps = eps
        self.d = d

    ## Creates an nxn Hermitian matrix with real/imaginary values between -1 and 1.
    def Hermitian(self):
        c = 0
        x = np.random.uniform(-1, 1, self.d**2)
        H = np.array([0 + 0.j] * self.d**2).reshape((self.d, self.d))
        for i in range(self.d):
            H[i][i] = x[-(i + 1)] + 0.j
        for a in range(self.d - 1):
            for b in range(self.d - 1 - a):
                H[a][b + a + 1] = x[2*b + 2*c] + x[2*b + 2*c + 1]*j
                H[b + a + 1][a] = x[2*b + 2*c] - x[2*b + 2*c + 1]*j
            c += self.d - 1 - a
        return H
    
    ## Creates an nxn Unitary matrix by estimating e^(iH) to 25 terms.
    def Unitary(self):
        H = self.Hermitian()
        U = self.expEst(H, 25)
        return U

    ## Creates an nxn SU(n) matrix by normalizing.
    def SUn(self):
        H = self.Hermitian()
        U = self.expEst(H, 25)
        SU3 = U / np.linalg.det(U)**(1/self.d)
        return SU3

    ## Estimates e^(i*a*M) for some matrix M and some constant a to n terms.
    def expEst(self, H, n):
        expTot = 0
        for k in range(n):
            expTot += (j * self.eps)**k / math.factorial(k) * np.linalg.matrix_power(H, k)
        return expTot

    ## Gives the nth column of matrix.
    def col(self, n, matrix):
        col = []
        for i in range(self.d):
            col.append(matrix[i][n])
        return col

    ## (WIP: Does not work yet!) Orthonormalizes the columns of a dxd matrix.
    def Gram_Schmidt(self, mat):
        n = self.d
        nmat = np.array([0. + 0.j] * n**2).reshape((n,n))
        for i in range(n):
            nmat[i][0] = self.col(0, mat)[i]
        for i in range(n):
            for j in range(n):
                nmat[j][i] = self.col(i, mat)[j]
                for k in range(i):
                    nmat[k][i] += np.vdot(self.col(i, nmat), self.col(k, nmat)) \
                                / np.vdot(self.col(k, nmat), self.col(k, nmat)) \
                                * self.col(k, nmat)[k]
        norms = []
        for i in range(n):
            norms.append(np.linalg.norm(self.col(i, nmat)))
        for i in range(n):
            for j in range(n):
                nmat[i][j] /= norms[i] 
        return nmat

    ## Returns a dxd identity matrix.  
    def I(self):
        Id = self.zero()
        for i in range(self.d):
            Id[i][i] = 1 + 0.j
        return Id

    ## Creates a dxd zero matrix.
    def zero(self):
        return np.array([0 + 0.j] * self.d**2).reshape((self.d, self.d))

## Common matrix constants.
class Matrix_Const(object):
    def __init__(self):
        ## Pauli Matrices.
        self.pm1 = np.array([0 + 0.j, 1 + 0.j,
                             1 + 0.j, 0 + 0.j])
        self.pm2 = np.array([0 + 0.j, 0 - 1.j,
                             0 + 1.j, 0 + 0.j])
        self.pm3 = np.array([1 + 0.j, 0 + 0.j,
                             0 + 0.j, -1 + 0.j])
        self.pmList = np.array([self.pm1, self.pm2, self.pm3])

        ## Gell-Mann Matrices.
        self.gm1 = np.array([0 + 0.j,  1 + 0.j,  0 + 0.j,
                             1 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j,  0 + 0.j]).reshape((3,3))
        self.gm2 = np.array([0 + 0.j,  0 - 1.j,  0 + 0.j,
                             0 + 1.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j,  0 + 0.j]).reshape((3,3))                         
        self.gm3 = np.array([1 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j, -1 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j,  0 + 0.j]).reshape((3,3))
        self.gm4 = np.array([0 + 0.j,  0 + 0.j,  1 + 0.j,
                             0 + 0.j,  0 + 0.j,  0 + 0.j,
                             1 + 0.j,  0 + 0.j,  0 + 0.j]).reshape((3,3))
        self.gm5 = np.array([0 + 0.j,  0 + 0.j,  0 - 1.j,
                             0 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 1.j,  0 + 0.j,  0 + 0.j]).reshape((3,3))
        self.gm6 = np.array([0 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j,  1 + 0.j,
                             0 + 0.j,  1 + 0.j,  0 + 0.j]).reshape((3,3))
        self.gm7 = np.array([0 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j,  0 - 1.j,
                             0 + 0.j,  0 + 1.j,  0 + 0.j]).reshape((3,3))
        self.gm8 = np.array([1 + 0.j,  0 + 0.j,  0 + 0.j,
                             0 + 0.j,  1 + 0.j,  0 + 0.j,
                             0 + 0.j,  0 + 0.j, -2 + 0.j]).reshape((3,3)) / np.sqrt(3)
        self.gmList = np.array([self.gm1, self.gm2, self.gm3,
                                self.gm4, self.gm5, self.gm6,
                                self.gm7, self.gm8])

## Main program, runs if run = True.         
if __name__ == "__main__" and run:
    numMatrices = 100000
    M = Matrix_nxn(3) ## What size matrix you want?
    SU3list = []
    detsRe = []
    detsIm = []

    ## Creates numMatrices matrices.
    for i in range(numMatrices):
        SU3list.append(M.SUn())

    ##Creates lists for the real and imaginary parts of the determinants of
    ## each matrix.
    for i in SU3list:
        detsRe.append(np.real(np.linalg.det(i)))
        detsIm.append(np.imag(np.linalg.det(i)))

    ## Prints out the averages/standard deviations of these lists, used
    ## primarily for troubleshooting.
    print('Imaginary')
    print('Average:', np.average(detsIm))
    print('Standard Deviation:', np.std(detsIm), '\n')

    print('Real:')
    print('Average:', np.average(detsRe))
    print('Standard Deviation:', np.std(detsRe), '\n')

    ## Plots the real part of the trace of each matrix versus the imaginary
    ## part. Makes a cool concave triangle thing for d=3. Generally, creates
    ## a concave shape with d points but is too fuzzy past d = 5 for anything
    ## cool. :(    
    imTrace = []
    reTrace = []
    for i in range(numMatrices):
        imTrace.append(np.imag(np.trace(SU3list[i])))
        reTrace.append(np.real(np.trace(SU3list[i])))

    import matplotlib.pyplot as plt
    plt.figure()
    plt.xlabel('ReTr(M)')
    plt.ylabel('ImTr(M)')
    plt.xlim(-2, 3)
    plt.ylim(-3, 3)
    plt.scatter(reTrace, imTrace, s = 1, color = 'b')

    ## For print matrices without a million decimal points...
    '''
    matrix = np.array([x.__format__('.3f') for x in matrix.flatten()]).reshape((3,3))
    '''

## Some crap. (Grahm-Schmidt process for 3x3, if you must know)
'''
def orthonormalize(self, matrix):
    newmatrix = ((np.zeros(9) + 0 * j).reshape((3,3)))

    ## Replaces the first column
    for i in range(3):
        newmatrix[i][0] = self.col(0, matrix)[i]

    ## Replaces the second column
    for i in range(3):
        newmatrix[i][1] = self.col(1, matrix)[i] - np.vdot(self.col(0, newmatrix), self.col(1, matrix))  \
                                                 / np.vdot(self.col(0, newmatrix), self.col(0, newmatrix))  \
                                                         * self.col(0, newmatrix)[i]

    ## Replaces the third column
    for i in range(3):
        newmatrix[i][2] = self.col(2, matrix)[i] - np.vdot(self.col(0, newmatrix), self.col(2, matrix))  \
                                                 / np.vdot(self.col(0, newmatrix), self.col(0, newmatrix))  \
                                                         * self.col(0, newmatrix)[i]  \
                                                 - np.vdot(self.col(1, newmatrix), self.col(2, matrix))  \
                                                 / np.vdot(self.col(1, newmatrix), self.col(1, newmatrix))  \
                                                         * self.col(1, newmatrix)[i]
    
    ## Normalizes the columns   
    normCol1 = np.linalg.norm(self.col(0, newmatrix))
    normCol2 = np.linalg.norm(self.col(1, newmatrix))
    normCol3 = np.linalg.norm(self.col(2, newmatrix))

    for i in range(3):
        newmatrix[i][0] /= normCol1
        newmatrix[i][1] /= normCol2
        newmatrix[i][2] /= normCol3

    return newmatrix
'''
