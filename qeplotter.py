#!/usr/bin/env python
'''
This script reads data from the output band file produced by QE. It can be 
used to plot electronic and phononic band-strcuture.
Initially written on 23 Jul 2013
Functionized on 26 Jan 2015
Added an option to compare between different dispersions on 2 Jun 2015
Author: Hongze Xia
Affiliation: University of New South Wales
'''
import numpy as np
import pylab as plt

class ReadBand(object):
    '''
    This class reads data from the band files produced by Quantum espresso. It basically
    contains the dispersions and the k-points associated.
    nbnd: int, number of bands for one k-point
    nks: total number of k-points
    kvec: matrix of kvectors
    band: matrix of dispersions
    '''
    def __init__(self,filename):
        # validation of the file's name
        assert type(filename) is str
        f = open(filename,'r')
        firstline = f.readline()
        tmp =  firstline.split()
        self.nbnd = int(tmp[2][0:len(tmp[2])-1]) # tmp[2] is '12,'
        self.nks = int(tmp[-2])
        # print nbnd, nks
        self.kvec = np.zeros(shape=(self.nks,3))
        self.band = np.zeros(shape=(self.nks,self.nbnd))
        n = 0
        tmpe = []
        for line in f:
            # if the first few characters are blank, then it is a k vector
            if '        ' in line:
                self.kvec[n] = np.fromstring(line,dtype=float,sep=' ')
                n += 1
            # else it is the eigenvalues
            else:
                e = np.fromstring(line,dtype=float,sep=' ')
                if tmpe == []:
                    tmpe = e
                else:
                    tmpe = np.append(tmpe,e)
        # print tmpe.shape
        self.band = np.reshape(tmpe,self.band.shape)

def plot(infile,Q,point_names,efermi=0.0,edos=None,ymin=None,\
        ymax=None,dosarray=None,phonon=False,phdos=None,\
        phunit="cm-1",comp=None,outfile="band.pdf",flip=False):
    """
    infile: str, QE output band file
    Q: list of integer, indices of spacial points starting
        from zero (C or Python order)
    point_names: list of str, the corresponding names of Q
    efermi: float, Fermi level for electron band structure
    edos: str, QE electron DOS output file (three columns,
        namely energy, DOS, int DOS. The first two are used)
    ymin,ymax: float, limits of y-axis
    dosarray: 2D numpy array, one could supply DOS as array
        if it is desirable
    phonon: boolean, False for electron plot; otherwise,
        it will be a phonon plot (affect y-axis)
    phdos: str, QE phonon DOS output file (two columns)
    phunit: str, can be "meV" or "mev", "cm-1" or "cm",
        "THz" or "thz".
    comp: str, compare with other models. first column is
        k-path in unit of 2pi/alat and the rest is eigenvalues.
        It assumes csv file is comma delimited.
    outfile: str, output band structure file
    flip: boolean
        DFT in solid line and external in dashed line.
    """
    # some checks
    if len(Q) != len(point_names):
        raise ValueError ("Check your special points!")
    if edos != None:
        print "Reading electron DOS from file."
        assert type(edos) is str
        dosarray = np.genfromtxt(edos)[:,0:2]
    if phdos != None:
        print "Reading phonon DOS from file."
        assert type(phdos) is str
        dosarray = np.genfromtxt(phdos)
    if dosarray != None:
        print "Plotting dispersions along with DOS."
    if comp != None:
        print "Plot will overlap with "+comp
        import os.path
        extension = os.path.splitext(comp)[1]
        if extension.lower() != ".csv":
            dt = np.genfromtxt(comp)
        else:
            dt = np.genfromtxt(comp,delimiter=",")
        q1 = dt[:,0]; b1 = dt[:,1:]
        del os
    if not flip:
        dftLine = 'k-'
        extLine = 'k--'
    else:
        dftLine = 'k--'
        extLine = 'k-'
        
    b0 = ReadBand(filename=infile)
    b0.band = np.sort(b0.band)
    # construct the reciprocal path
    sp = b0.kvec[Q] # spacial points
    q = np.zeros(b0.nks)
    start = 0.0
    for i in range(len(Q)-1):
        nks = Q[i+1] - Q[i] + 1
        length = np.linalg.norm(sp[i+1] - sp[i])
        q[Q[i]:Q[i+1]+1] = np.linspace(start,start+length,nks,endpoint=True)
        start += length
    
    plt.figure(0,(10, 6))
    if dosarray != None: plt.axes([.11, .07, .64, .85])
    plt.xlabel("Reduced wave number", fontsize=18)
    # some conversions
    if phonon == False:
        print "The electronic band structure will be plotted."
        print "All energy will be normalised to E_F = %8.3f eV." % efermi
        b0.band -= efermi
        if dosarray != None: dosarray[:,0] -= efermi
        if ymax == None:
            ymax = np.ceil(b0.band.max()/5.0)*5
        if ymin == None:
            ymin = np.floor(b0.band.min()/5.0)*5
        plt.plot(q,[0.0]*len(q),'r--',lw=1) # Fermi level
        plt.ylabel("Electron Energy ($\mathrm{eV}$)", fontsize=18)
    else:
        print "The phonon dispersions will be plotted."
        if phunit.lower() == "mev":
            b0.band /= 8.0532
            if dosarray != None: dosarray[:,0] /= 8.0532
            plt.ylabel("Phonon energy (meV)", fontsize=18)
            if ymax == None:
                ymax = np.ceil(b0.band.max()/10.0)*10
            if ymin == None:
                ymin = np.floor(b0.band.min()/10.0)*10
        if phunit.lower() == "thz":
            b0.band /= 33.33333
            if dosarray != None: dosarray[:,0] /= 33.33333
            plt.ylabel("Phonon frequency (THz)", fontsize=18)
            if ymax == None:
                ymax = np.ceil(b0.band.max()/5.0)*5
            if ymin == None:
                ymin = np.floor(b0.band.min()/5.0)*5
        if phunit == "cm-1" or phunit == "cm":
            plt.ylabel("Phonon frequency ($\mathrm{cm^{-1}}$)", fontsize=18)
            if ymax == None:
                ymax = np.ceil(b0.band.max()/100.)*100
            if ymin == None:
                ymin = np.floor(b0.band.min()/100.)*100
    # plot the figure and some tweeks
    for n in range(b0.nbnd):
        plt.plot(q, b0.band[:,n], dftLine, lw=1)
    if comp != None:
        for n in range(len(b1[0,:])):
            plt.plot(q1,b1[:,n], extLine, lw=1)
    plt.xticks(q[Q], point_names)
    plt.tick_params(axis='x', labeltop='on',labelsize=15,labelbottom='off')
    plt.yticks(fontsize=15)
    plt.xlim(q[0], q[-1])
    plt.ylim(ymin,ymax)
    plt.grid('on')
    # plot DOS if enabled
    if dosarray != None:
        # crop the dosarray
        mask = (dosarray[:,1]>=ymin)*(dosarray[:,1]<=ymax)
        dosarray = dosarray[mask]
        # such that it fills solid color
        dosarray[0,1] =0.0; dosarray[-1,1] =0.0
        plt.axes([.78, .07, 0.17, .85])
        plt.plot(dosarray[:,1],dosarray[:,0],'k-',lw=1)
        plt.fill(dosarray[:,1],dosarray[:,0],color='lightgrey')
        plt.ylim(ymin,ymax)
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel("DOS",fontsize=18)
    
    plt.savefig(outfile,dpi=300)