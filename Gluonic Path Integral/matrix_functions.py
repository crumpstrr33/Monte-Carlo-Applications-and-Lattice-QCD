'''
In this file, the class MatrixNxN allows for certain operations and matrix 
creation for a specific d. The class MatrixConst contains some common matrix
constants that may be needed.
'''
import numpy as np
import cmath
import math

j = cmath.sqrt(-1)

## Run the main program?
RUN = False


class MatrixNxN(object):
    '''
    Various functions for matrix creation/manipulation for dxd matrices.
    '''
    def __init__(self, d, eps=1):
        '''
        Parameters:
        d - the dimension of the matrix (dxd)
        eps - (optional) a multiplicative factor when estimating the
              exponential of a matrix, defaults to 1.
        '''
        self.eps = eps
        self.d = d


    def hermitian(self, min_val=-1, max_val=1):   
        '''
        Parameters:
        min_val - (optional) the minimum value allowed for the real and
                  imaginary part (is not the min value for the norm)
        max_val - (optional) the maximum value allowed for the real and
                  imaginary part (is not the max value for the norm)
        Returns:
        H - a complex dxd Hermitian matrix

        Creates an nxn Hermitian matrix with real/imaginary values uniformly
        distributed between min_val and max_val.
        '''
        H = np.random.uniform(min_val, max_val, (self.d, self.d)) + \
            np.random.uniform(min_val, max_val, (self.d, self.d)) * j
        H = (H + np.conj(H.T)) / 2

        return H


    def unitary(self):
        '''
        Returns:
        U - a complex dxd Unitary matrix

        Creates an dxd Unitary matrix by estimating e^(iH) where H is a
        Hermitian matrix generated by the function self.hermitian.
        '''
        H = self.hermitian()
        U = self._exp_est(H, 25)

        return U


    def sun(self):
        '''
        Returns:
        SUd - a complex dxd SU(d) matrix

        Creates an SU(n) (for n=d) matrix by normalizing a unitary matrix
        with the dth root of the unitary matrix's determinant where d is the
        size of the matrix.
        '''
        H = self.hermitian()
        U = self._exp_est(H)
        SUd = U / np.linalg.det(U)**(1 / self.d)

        return SUd


    def _exp_est(self, H, n=25):
        '''
        Parameters:
        H - a dxd matrix. The variable H is used because this function is used
            for calculating a Unitary matrix by exponentiating a Hermitian
            matrix
        n - (optional) the number of terms to estimate the sum to
        Returns:
        exp_tot - the estimation of the exponential

        Estimates e^(i*a*M) by expanding it's Taylor seris for some matrix M
        and some constant a to n terms. With n = 25, it's a good balance, as
        the error and the computational time both very small, but it can be
        changed as n is an input variable.
        '''
        exp_tot = 0

        for k in range(n):
            exp_tot += (j * self.eps)**k / math.factorial(k) * \
                                           np.linalg.matrix_power(H, k)
        return exp_tot


    def I(self, dtype=complex):
        '''
        Parameters:
        dtype - (optional) data type for the matrix

        Returns a dxd complex identity matrix.  
        '''
        return np.eye(self.d, dtype=dtype)

  
    def zero(self, dtype=complex):
        '''
        Parameters:
        dtype - (optional) data type for the matrix
        
        Creates a dxd complex zero matrix.
        '''  
        return np.zeros((self.d, self.d), dtype=dtype)


'''
Main function: Shows the distrubtion of the determinants of SU(n) matrices
               on a complex plot.
'''
def main():
    num_mat = 10000
    su3_list = []
    dets_re = []
    dets_im = []
    
    ## What size matrix you want?
    M = MatrixNxN(3) 

    ## Creates numMatrices matrices.
    for i in range(num_mat):
        su3_list.append(M.sun())

    ##Creates lists for the real and imaginary parts of the determinants of
    ## each matrix.
    for i in su3_list:
        dets_re.append(np.real(np.linalg.det(i)))
        dets_im.append(np.imag(np.linalg.det(i)))

    ## Prints out the averages/standard deviations of these lists, used
    ## primarily for troubleshooting.
    print('Imaginary')
    print('Average:', np.average(dets_im))
    print('Standard Deviation:', np.std(dets_im), '\n')

    print('Real:')
    print('Average:', np.average(dets_re))
    print('Standard Deviation:', np.std(dets_re), '\n')

    ## Plots the real part of the trace of each matrix versus the imaginary
    ## part. Makes a cool concave triangle thing for d=3. Generally, creates
    ## a concave shape with d vertices.
    im_trace = []
    re_trace = []
    for i in range(num_mat):
        im_trace.append(np.imag(np.trace(su3_list[i])))
        re_trace.append(np.real(np.trace(su3_list[i])))

    import matplotlib.pyplot as plt
    plt.figure()
    plt.xlabel('ReTr(M)')
    plt.ylabel('ImTr(M)')
    plt.xlim(-2, 3)
    plt.ylim(-3, 3)
    plt.scatter(re_trace, im_trace, s = 1, color = 'b')

    ## To print 'matrix' with a limited number of decimal points
    # np.array([x.__format__('.3f') for x in matrix.flatten()]).reshape((3,3))

  
if __name__ == "__main__":
    ## Runs the main function based on value of RUN
    if RUN:
        main()