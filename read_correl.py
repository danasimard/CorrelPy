# read_correl.py
# Dana Simard January 2015
#
# Module to read in correlation files and return k, delta and the 
# uncertainty to delta.  read_correl is the main function that 
# does this.  Other functions are called from read_correl.
# Note that for the equal time correlators, the poisson noise is 
# subtracted from delta2, so delta = delta2 - delta2poisson.
# For the unequal time correlators, this does not happen, so 
# delta = delta2

import numpy as np

def __init__():
    return

# Main function
def read_correl(dirname, checkpoint1, checkpoint2):
    if checkpoint1 == checkpoint2:
        k, delta, edelta = read_auto( dirname, checkpoint1)
    else:
        k, delta, edelta = read_cross( dirname, checkpoint1, checkpoint2 )
    return k, delta, edelta

# Function for reading in delta, edelta and k from autocorrelators
# To find delta, this function subtracts the poisson noise from delta2
def read_auto(dirname, checkpoint):
    datapath = dirname + '/' + "{0:.3f}".format(checkpoint) + 'ngpps.dat'
    dt = [('k', float),('delta2',float),('edelta2',float),('delta2poisson',float),('edelta2poisson',float)]
    d = np.loadtxt( datapath, dtype=np.dtype(dt) )
    delta = np.array( d['delta2'] ) - np.array( d['delta2poisson'] )
    edelta = (np.array( d['edelta2'])**2. + np.array( d['edelta2poisson'] )**2.)**0.5
    return np.array(d['k']), delta, edelta

# Function for reading in delta, edelta and k from cross correlators
# Doesn't make any effort to subtract noise - assume noise cancels when
# taking cross correlation spectrum
def read_cross(dirname, checkpoint1, checkpoint2):
    datapath1 = dirname + '/' + "{0:.3f}".format(checkpoint1) + 'x' + "{0:.3f}".format(checkpoint2) + 'ngpps.dat'
    datapath2 = dirname + '/' + "{0:.3f}".format(checkpoint2) + 'x' + "{0:.3f}".format(checkpoint1) + 'ngpps.dat'
    try:
        f = open( datapath1 )
    except IOError:
        f = open( datapath2 )
    lines = f.readlines()
    k = []
    delta = []
    edelta = []
    for line in lines:
        items = line.split()
        k.append( float( items[0] ) )
        delta.append( float( items[1] ) )
        edelta.append( float( items[2] ) )
    f.close()
    return np.array(k), np.array(delta), np.array(edelta)


        
