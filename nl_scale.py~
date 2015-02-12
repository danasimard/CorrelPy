def __init__():
    return

from plot_eigen import *
import numpy as np 
import matplotlib.pyplot as plt
import string

# checkpoints needs to be the same list of checkpoints as used in the 
# parameter file in diag_cross_spectrum. It would be ideal to later 
# change diag_cross_spectrum to print the parameter file to the 
# data directory so that it can be read in by read_eigenvalues and
# read_eigenvectors
def eigenvalue_nonlinear_scale( dirname ):
    numcheckpoints,numpoints,k,eigenvalues = read_eigenvalues( dirname )
    eigen_t = numcheckpoints*0.5
    index = np.abs( eigenvalues[-1,:] - eigen_t ).argmin()
    k_t = k[index] 
    return k_t, eigenvalues[:,0]


def eigenvalue_old( dirname, checkpoints = [] ):
    output = dirname + '/eigenvalues_nonlinear.ps'
    numcheckpoints, numpoints, k, eigenvalues = read_eigenvalues( dirname )
    if len(checkpoints)<1:
        checkpoints = np.arange( numcheckpoints )
    eigen_t = numcheckpoints * 0.5
    k_t = np.zeros( numcheckpoints )
    for i in range(numcheckpoints):
        eigenvalues_now = eigenvalues[:,i]
        index = np.abs( eigenvalues_now - eigen_t ).argmin()
        k_t[i] = k[index]
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.plot( checkpoints, k_t )
    plt.ylabel( 'k_NL' )
    plt.xlabel( 'redshift' )
    plt.title( dirname )
    plt.savefig( output )
    fig.show()
    return checkpoints, k_t
