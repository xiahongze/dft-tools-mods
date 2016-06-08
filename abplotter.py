#!/usr/bin/env python
# This is a python script for plotting the band structure and density of states
# from ABINIT. It reads the file ends with EIG and DOS. The user should provide
# the file names, the lattice vectors and special points' names and positions.
# Created by Hongze Xia, 30 OCT 2015.

import numpy as np
import matplotlib.pylab as plt
ha2ev = 2.7210E1
def eig_reader (fn):
    """
    fn: str, file name
    return:
        (eig,kpts)
    print the number of bands and points
    """
    if type(fn) is not str:
        print "This is not a valid input variable for fn."
    with open(fn,"r") as fp:
        # read the first line
        line1 = fp.readline()
        # print line1.split()
        nkpt = int(line1.split()[-3])
        kpts = []; eig = []
        for line in fp:
            split = line.split()
            if split[0] == "kpt#":
                kpts.extend(map(float,split[7:10]))
            else:
                eig.extend(map(float,split))
        eig = np.array(eig).reshape(nkpt,-1)
        eig = np.sort(eig,axis=1)
        kpts = np.array(kpts).reshape(nkpt,3)
    print "===> %d k-points and %d bands in File \"%s\"" % (eig.shape+(fn,))
    return (eig,kpts)

def dos_reader(fn,pdos=False):
    """docstring for fname"""
    if type(fn) is not str:
        print "This is not a valid input variable for fn."
    with open(fn,"r") as fp:
        for i in range(7): # skip the first few lines
            fp.readline()
        line = fp.readline()
        split = line.split()
        if split[1] == "Fermi":
            ef = float(split[-1])*ha2ev
    dt = np.genfromtxt(fn)
    if not pdos:
        print "===> Fermi energy : %10.5f eV from File \"%s\"" % (ef,fn)
        e = dt[:,0]*ha2ev
        dos = dt[:,1]/ha2ev
        return (ef,e,dos)
    else:
        e = dt[:,0]*ha2ev - ef
        dos = dt[:,1:6]/ha2ev
        return (e,dos)

def kline(kpts,a):
    """return a k-line for the kpts
    kpts: crystal coordinates
    a: lattice vectors
    return: kline
    """
    b = np.linalg.inv(a).T
    kconv = kpts.dot(b)
    segments = np.linalg.norm(kconv[1:] - kconv[:-1],axis=1)
    return np.cumsum(np.hstack((0,segments)))

def plot(a,fileig,Q,point_names,ef=0,fildos=None,ymin=None,ymax=None
        ,pdos_pref=None,atoms=None,pdos_max=None):
    """docstring for plot"""
    fig = plt.figure(0,(12, 8))
    assert len(Q) == len(point_names), "Length of Q and point_names should be the same!"
    if fildos != None:
        ax = plt.axes([.06, .05, .7, .85])
        ef,e,dos = dos_reader(fildos)
        e -= ef
    elif pdos_pref:
        ax = plt.axes([.06, .05, .7, .85])
    else:
        ax = fig.add_subplot(111)
        
    eig,kpts = eig_reader(fileig)
    eig -= ef
    q = kline(kpts,a)
    for i in xrange(eig.shape[1]):
        plt.plot(q,eig[:,i],'k-',lw=1)
    plt.xticks(q[Q], point_names)
    plt.tick_params(axis='x', labeltop='on',labelsize=15,labelbottom='off')
    plt.yticks(fontsize=15)
    plt.xlim(q[0], q[-1])
    plt.plot(q,[0]*len(q),'r--')
    plt.xlabel("Reduced wave number", fontsize=18)
    plt.ylabel("Electron Energy (eV)", fontsize=18)
    plt.grid('on')
    plt.ylim(ymin,ymax)
    ymin,ymax = plt.ylim()
    # add an extra ef on the right of the axis.... so troublesome
#    ax1 = plt.axes([.11, .05, .64, .85],frameon=False)
    ax1 = plt.axes(ax.get_position(),frameon=False)
    ax1.yaxis.tick_right()
    # plt.tick_params(axis='y', labelleft='on',labelright='on',labelsize=15)    
    ax1.xaxis.set_ticklabels([])
    ax1.xaxis.set_ticks_position('none')
    plt.yticks([0],['$\epsilon_{\mathrm{F}}$'],fontsize=15)
    plt.ylim(ymin,ymax)

    if fildos:
        plt.axes([.79, .05, 0.20, .85])
        plt.plot(dos,e,'k-',lw=1)
        plt.ylim(ymin,ymax)
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel("DOS",fontsize=18)
    
    if pdos_pref:
        ax2 = plt.axes([.79, .05, 0.20, .85])
        plot_pdos(pdos_pref,atoms,ax=ax2)
        ax2.set_xticks([],[])
        ax2.set_yticks([],[])
        ax2.set_ylim(ymin,ymax)
        ax2.set_xlim(0,pdos_max)
    plt.show()
    plt.close()

def plot_pdos(prefix,atoms,ax=None,xmin=None,xmax=None,ymin=None,ymax=None):
    from glob import glob
    assert type(prefix) == str, "What is the prefix of your PDOS files?"
    filelist = glob(prefix+"*"); filelist.sort()
    del glob
    assert len(filelist) == len(atoms), "Atom number should match with the number of PDOS files!"
    print "===> Please check if the file list matches the atom list:"
    for i in xrange(len(atoms)):
        print "     %d. %s <--> %s" % (i,filelist[i],atoms[i])
    dt = map(lambda x: dos_reader(x,pdos=True),filelist)
    pdos = np.asarray(map(lambda x: x[1],dt)); e = dt[0][0]
    atomset = []
    for i in atoms:
        if i not in atomset:
            atomset.append(i)
    atoms = np.array(atoms)
    tdos = pdos.sum(axis=0).sum(axis=1)
    pdos_by_atom = np.array([pdos[atoms==atom].sum(axis=0) for atom in atomset])
    orbitals = ['$s$','$p$','$d$','$f$','$g$']
    # for i in xrange(1,len(dt)):
    #     tdos += dt[i][1]
    # atomset = set(atoms)
    if ax == None:
        fig = plt.figure(0,(10, 6))
        ax = plt.gca()
        ax.plot(e,tdos,label="Total DOS")
        for i in xrange(len(atomset)):
            ax.plot(e,pdos_by_atom[i].sum(axis=1),label=atomset[i])
        show = True
        ax.set_xlabel("Electron energy (eV)",fontsize=18)
        ax.set_ylabel("DOS (a.u.)",fontsize=18)
        ax.legend(ncol=1)
    else:
        ax.plot(tdos,e,label="Total DOS")
        for i in xrange(len(atomset)):
            ax.plot(pdos_by_atom[i].sum(axis=1),e,label=atomset[i])
        show = False
        ax.set_xlabel("DOS (a.u.)",fontsize=18)
        ax.legend(bbox_to_anchor=(0., 1.01, 1., .10), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
        
        # for j in xrange(5):
            # plt.plot(e,pdos_by_atom[i,:,j],label="%s %s" %(atomset[i],orbitals[j]))
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    if show:
        plt.show()
