import numpy as np
import matplotlib.pyplot as plt
import string
import matplotlib as mpl

def __init__():
    return

def read_parameters( dirname ):
    infile = dirname + '/correls_parameters.txt'
    f = open( infile )
    lines = f.readlines()[1:]
    checkpoints = []
    for line in lines:
        items = line.split()
        checkpoints.append( float(items[0]) )
    f.close()
    return checkpoints

def plot_eigenvalues( dirname ):
    outpath = dirname + '/coefficients_eigenvalues.ps'
    numcheckpoints, numpoints, k, eigenvalues = read_eigenvalues( dirname )
    eigenvalues = np.array( eigenvalues )
    eigenvalues_n = eigenvalues / numcheckpoints
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(numcheckpoints):
        plt.plot( k, eigenvalues_n[i] )
    plt.plot( k, np.zeros( numpoints), 'k-' )
    ax.set_xscale('log')
    plt.xlabel( 'k' )
    plt.title( dirname )
#    mpl.rcParams['text.usetex']=True
#    plt.ylabel( r"$\lambda$" )
    plt.ylabel( 'lambda' )
    plt.savefig( outpath )
    fig.show()
    return


# Use this routine - this is the plot we want to see when we 
# plot the eigenvectors - basically comparing the ith component of 
# each of the eigenvectors
def plot_eigenvectors(dirname, k_value ):
    numcheckpoints, numpoints, k, eigenvectors = read_eigenvectors(dirname)
    checkpoints = read_parameters( dirname )
    eigenvectors = np.array( eigenvectors )
    eigenvectors_n = eigenvectors * np.sqrt( numcheckpoints )
    k_index = np.abs( k - k_value ).argmin()
    output = dirname + '/coefficients_eigenvectors_' + str(k_value) + '.ps'
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(numcheckpoints):
        plt.plot( checkpoints, eigenvectors_n[k_index,i] )
    plt.ylabel( 'eigenvectors * sqrt[n]' )
    plt.xlabel( 'redshift' )
    plt.title( dirname + '   k = ' + str(k[k_index] ) )
    #ax.set_xscale('log')
    plt.gca().invert_xaxis()
    plt.savefig( output )
    fig.show()
    return
    
    

def read_eigenvalues( dirname ):
    datapath = dirname + '/coefficients_eigenvalues.dat'
    f = open(datapath)
    lines = f.readlines()
    line = lines[0]
    items = line.split()
    numpoints = int(items[0])
    numcheckpoints = int(items[1])
    k = np.zeros( numpoints )
    eigenvalues = np.zeros( (numcheckpoints, numpoints) )
    for i in range( numpoints):
        items = lines[i+1].split()
        k[i] = items[0]
        for j in range( numcheckpoints ) :
            eigenvalues[j,i] = items[j+1]
    f.close()
    return numcheckpoints, numpoints, k, eigenvalues

def plot_eigenvectors2(dirname, k_value):
    numcheckpoints, numpoints, k, eigenvectors = read_eigenvectors2(dirname)
    k_index = np.abs( k - k_value ).argmin()
    output = dirname + '/coefficients_eigenvectors_' + str(k_value) + '_2.ps'
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(numcheckpoints):
        plt.plot( range(numcheckpoints), eigenvectors[k_index,i] )
    plt.ylabel( 'eigenvectors' )
    plt.xlabel( 'checkpoint' )
    plt.title( dirname + '   k = ' + str(k[k_index] ) )
    plt.savefig( output )
    fig.show()
    return

def read_eigenvectors( dirname ):
    datapath = dirname + '/coefficients_eigenvectors.dat' 
    f = open(datapath)
    lines = f.readlines()
    items = lines[0].split()
    numpoints = int(items[0])
    numcheckpoints = int(items[1])
    k = np.zeros(numpoints)
    eigenvectors = np.zeros( (numpoints, numcheckpoints, numcheckpoints) )
    for l in range( numpoints ):
        items = lines[l+1].split()
        k[l] = items[0]
        for i in range( numcheckpoints ):
            for j in range( numcheckpoints ):
                eigenvectors[l,i,j] = items[i*numcheckpoints+j+1]
    f.close()
    return numcheckpoints, numpoints, k, eigenvectors

def read_eigenvectors2( dirname ):
    datapath = dirname + '/coefficients_eigenvectors.dat' 
    f = open(datapath)
    lines = f.readlines()
    items = lines[0].split()
    numpoints = int(items[0])
    numcheckpoints = int(items[1])
    k = np.zeros(numpoints)
    eigenvectors = np.zeros( (numpoints, numcheckpoints, numcheckpoints) )
    for l in range( numpoints ):
        items = lines[l+1].split()
        k[l] = items[0]
        for i in range( numcheckpoints ):
            for j in range( numcheckpoints ):
                eigenvectors[l,j,i] = items[i*numcheckpoints+j+1]
    f.close()
    return numcheckpoints, numpoints, k, eigenvectors
    
