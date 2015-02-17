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


def plot_contribution( dirname, k_value ):
    outpath = dirname + '/coefficients_contribution_' + str( k_value ) + '.ps'
    numcheckpoints, numpoints, k, eigenvalues = read_eigenvalues( dirname )
    numcheckpoints, numpoints, k, eigenvectors = read_eigenvectors( dirname, 0 )
    checkpoints = read_parameters( dirname )
    eigenvectors = np.array( eigenvectors )
    eigenvectors_n = np.zeros( (numpoints, numcheckpoints, numcheckpoints) )
    for l in range( numpoints ):
        for i in range( numcheckpoints):
            eigenvectors_n[l,i,:] = eigenvectors[l,i,:] *np.sqrt( abs(eigenvalues[i,l]) )
    k_index = np.abs( k - k_value ).argmin()
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(numcheckpoints):
        plt.plot( checkpoints, eigenvectors_n[k_index,i] )
    plt.ylabel( 'v * sqrt[|lambda|]' )
    plt.xlabel( 'redshift' )
    plt.title( dirname + '   k = ' + str(k[k_index] ) )
    #ax.set_xscale('log')
    plt.gca().invert_xaxis()
    plt.savefig( outpath )
    fig.show()
    return


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
# Note that when tauflag is 1, the code uses H0=67.  If you want use a 
# different H0, you will have to change this later
def plot_eigenvectors(dirname, k_value, wflag, tauflag ):
    numcheckpoints, numpoints, k, eigenvectors = read_eigenvectors(dirname, wflag)
    checkpoints = read_parameters( dirname )
    if (tauflag):
        checkpoints_old = np.array( checkpoints )
        checkpoints_old = 1.0 / ( 1.0 + checkpoints_old)
        checkpoints = -2./67. * ( checkpoints_old)**(-0.5)
    eigenvectors = np.array( eigenvectors )
    eigenvectors_n = eigenvectors * np.sqrt( numcheckpoints )
    k_index = np.abs( k - k_value ).argmin()
    if wflag:
        output = dirname + '/coefficients_wvectors_' + str(k_value) + '.ps'
    else:
        output = dirname + '/coefficients_eigenvectors_' + str(k_value) + '.ps'
    fig = plt.figure()
    ax = plt.subplot(111)
    eigenvectors_sum = np.zeros( (numpoints, numcheckpoints) )
    if wflag:
        for i in range(numcheckpoints):
            plt.plot( checkpoints, eigenvectors[k_index,i] )
            eigenvectors_sum = eigenvectors_sum + eigenvectors[:,i,:]
        plt.plot( checkpoints, eigenvectors_sum[k_index ], 'k-' )
        plt.ylabel( 'wvectors' )
        if tauflag:
            plt.xlabel( 'tau' )
        else:
            plt.xlabel( 'redshift' )
        plt.title( dirname + ' wvectors k = ' + str(k[k_index]))
    else:
        for i in range(numcheckpoints):
            plt.plot( checkpoints, eigenvectors_n[k_index,i] )
        plt.ylabel( 'eigenvectors * sqrt[n]' )
        if tauflag:
            plt.xlabel( 'tau' )
        else:
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

def read_eigenvectors( dirname, wflag ):
    if wflag:
        datapath = dirname +'/coefficients_wvectors.dat'
    else:
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
    
def read_delta( dirname ):
    datapath = dirname + '/delta_cs0.dat'
    f = open( datapath )
    lines = f.readlines()
    items = lines[0].split()
    numpoints = int( items[0] )
    numcheckpoints = int( items[1] )
    checkpoints = np.zeros( numcheckpoints )
    k = np.zeros( numpoints )
    delta = np.zeros( (numpoints, numcheckpoints ) )
    items = lines[1].split()
    for i in range( numcheckpoints ):
        checkpoints[i] = float( items[i] )
    for l in range(numpoints):
        items = lines[l+2].split()
        k[l] = float( items[0] )
        for i in range( numcheckpoints):
            delta[l][i] = float( items[i+1] )
    return k, checkpoints, delta, numpoints, numcheckpoints

def plot_delta( dirname ):
    k, checkpoints, delta, numpoints, numcheckpoints = read_delta( dirname )
    outpath = dirname + "/delta_cs0.ps"
    k = np.array( k )
    delta = np.array( delta )
    delta_true = np.zeros( (numpoints, numcheckpoints) )
    for i in range(numcheckpoints):
        delta_true[:,i] = ( 2 * np.pi**2. / ( k**3.))**0.5 * delta[:,i]
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(numcheckpoints):
        plt.plot( k, delta[:,i], label=checkpoints[i])
    plt.xlabel( 'k [h/Mpc]' )
    plt.ylabel( 'delta [(Mpc/h)^3/2]' )
    plt.title( dirname )
    plt.legend()
    plt.savefig( outpath )
    fig.show()
    return
    
