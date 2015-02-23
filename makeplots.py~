# makeplots.py
# Dana Simard January 2015
# 
# This file contains all of the plotting routines that have to do with the 
# cross and auto correlators.
#

from CorrelPy.read_correl import *
from CorrelPy.correl_factors import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import pylab
from matplotlib.patches import Rectangle
import matplotlib.cm as cmx
import matplotlib

def __init__():
    return

def plot_correl_factors( dirname, checkpoints, log_switch):
    # dirname = path to the directory with all of the correlation data files
    # checkpoints = array with all of the checkpoints for which we want 
    #               to compute the cross correlation factors
    # log_switch = 0 for linear y axis, 1 for log y axis
    # eg. > plot_correl_factors( 'cubep3m8', [0., 1., 10.], 0 )
    fig1 = plt.figure()
    ax1 = plt.subplot(111)
    ax1.set_xscale('log')
    if log_switch:
        ax1.set_yscale('log')
    ax1.set_xlabel('k')
    ax1.set_ylabel('r')
    if len(checkpoints) == 2:
        ax1.set_title( dirname + ' ' + str(checkpoints[0]) + 'x' + str(checkpoints[1]) )
    nplots = 0
    for i in range(len(checkpoints) ):
        nplots = nplots + i
    if nplots > 4:
        fig2 = plt.figure()
        ax2 = plt.subplot(111)
        ax2.set_xscale('log')
        ax2.set_xlabel('k')
        ax2.set_ylabel('r')
        fig3 = plt.figure()
        ax3 = plt.subplot(111)
        ax3.set_xscale('log')
        ax3.set_xlabel('k')
        ax3.set_ylabel('r')
        fig4 = plt.figure()
        ax4 = plt.subplot(111)
        ax4.set_xscale('log')
        ax4.set_xlabel('k')
        ax4.set_ylabel('r')
        if log_switch:
            ax2.set_yscale('log')
            ax3.set_yscale('log')
            ax4.set_yscale('log')
        nperplot = int(round(nplots/4.0))
    else:
        nperplot = nplots
    formats = ['bo','gv','r^','cs','mp','y*','kd']
    m = 0
    for i in range(len(checkpoints)):
        for j in range(i+1, len(checkpoints)):
            k, r, er = correl_factors( dirname, checkpoints[i], checkpoints[j])
            if m < nperplot:
                ax1.errorbar( k, r, yerr=er, fmt=formats[m])
            elif m < nperplot*2:
                ax2.errorbar( k, r, yerr=er, fmt=formats[m-nperplot])
            elif m < nperplot*3:
                ax3.errorbar( k, r, yerr=er, fmt=formats[m-nperplot*2])
            else:
                ax4.errorbar( k, r, yerr=er, fmt=formats[m-nperplot*3])
            m = m + 1
    if nplots > 4:
        fig1.subplots_adjust(wspace=0,hspace=0)
        fig1.show()
        fig2.subplots_adjust(wspace=0,hspace=0)
        fig2.show()
        fig3.subplots_adjust(wspace=0,hspace=0)
        fig3.show()
        fig4.subplots_adjust(wspace=0,hspace=0)
        fig4.show()
    else:
        fig1.subplots_adjust(wspace=0,hspace=0)
        fig1.show()
    return

def plot_correls( dirname, checkpoint, log_switch):
    k, delta, edelta = read_correl( dirname, checkpoint, checkpoint )
    delta = abs(delta)
    fig1 = plt.figure()
    ax1 = plt.subplot(111)
    if log_switch:
        ax1.set_yscale( 'log' )
        ax1.set_xscale('log')
    ax1.set_title( dirname + ', redshift ' + str(checkpoint) )
    ax1.set_xlabel( 'k' )
    ax1.set_ylabel('P(k)' )
    ax1.errorbar( k[:-1], delta[:-1], yerr = edelta[:-1], fmt = 'bo' )
    fig1.show()
    return

def plot_correl_factors_3d( dirname, k_value, checkpoints ):
    numcheckpoints = len(checkpoints)
    r = np.empty( ( numcheckpoints, numcheckpoints), np.float)
    er = np.empty( ( numcheckpoints, numcheckpoints), np.float )
    k, rtemp, ertemp = correl_factors(dirname, checkpoints[0], checkpoints[0] )
    k_index = np.abs( k - k_value ).argmin()
    for i in range( len(checkpoints) ):
        for j in range( len (checkpoints ) ):
            ktemp, rtemp, ertemp = correl_factors( dirname, checkpoints[i], checkpoints[j])
            r[i][j] = rtemp[k_index]
            er[i][j] = ertemp[k_index]
    rplot = r.flatten()
    erplot = er.flatten()
    x = []
    y = []
    for i in range(numcheckpoints):
        for j in range(numcheckpoints):
            x.append( checkpoints[i] )
            y.append( checkpoints[j] )
    fig = pylab.figure()
    ax = Axes3D(fig)
    ax.set_title( 'k = ' + str(k[k_index]) )
    ax.set_xlabel( 'Checkpoint a' )
    ax.set_ylabel( 'Checkpoint b' )
    ax.set_zlabel( 'R_ab' )
    #    ax = fig.add_subplot(111, projection='3d' )
    colours = [u'k',u'b',u'r',u'g',u'm',u'c',u'y']
    for i in range(numcheckpoints):
        p = ax.scatter( x[numcheckpoints*i:numcheckpoints*(i+1)], y[numcheckpoints*i:numcheckpoints*(i+1)], rplot[numcheckpoints*i:numcheckpoints*(i+1)], zdir=u'z', s=40, c=colours[i%len(colours)], marker='^')
        #for i in range( len( x ) ):
        #x2, y2, _ = proj3d.proj_transform( x[i], y[i], rplot[i], ax.get_proj() )
        #label = pylab.annotate( str(x[i]) + ", " + str(y[i]), xy = (x2, y2), xytext = (-5,5), textcoords = 'offset points', ha = 'right', va = 'bottom')
        #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0') )
    fig.show()
    return

def plot_correl_factors_2d( dirname, k_value, checkpoints ):
    numcheckpoints = len(checkpoints)
    r = np.empty( ( numcheckpoints, numcheckpoints), np.float)
    er = np.empty( ( numcheckpoints, numcheckpoints), np.float )
    k, rtemp, ertemp = correl_factors(dirname, checkpoints[0], checkpoints[0] )
    k_index = np.abs( k - k_value ).argmin()
    for i in range( len(checkpoints) ):
        for j in range( len (checkpoints ) ):
            ktemp, rtemp, ertemp = correl_factors( dirname, checkpoints[i], checkpoints[j])
            r[i][j] = rtemp[k_index]
            er[i][j] = ertemp[k_index]
    width = []
    width.append( 0. )
    for i in range( 1, len(checkpoints) ):
        width.append( (checkpoints[i] - checkpoints[i-1])/2 )
    width.append( width[ - 1 ] )
    print width
    cm = plt.get_cmap('Greys')
    cNorm = matplotlib.colors.Normalize( vmin = min( r.flatten() ), vmax = max( r.flatten() ))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap = cm)
    pylab.plot([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], color="black" )
    for i in range(  len(checkpoints ) ):
        for j in range( len(checkpoints) ):
            pylab.gca().add_patch( Rectangle((checkpoints[i]-width[i],checkpoints[j]-width[j]),width[i]+width[i+1], width[j]+width[j+1], facecolor=scalarMap.to_rgba(r[i][j])))
    scalarMap.set_array(r.flatten() )
    pylab.colorbar( scalarMap )
    pylab.title( 'k = ' + str(k[k_index]) )
    pylab.xlabel( "z_a" )
    pylab.ylabel( "z_b" )
    pylab.show()
    return

            
