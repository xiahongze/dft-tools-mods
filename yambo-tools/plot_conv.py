#!/usr/bin/python

import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np
mpl.style.use('ggplot')


def plot(ax, taskids, conv_param):
    """
    conv_param: str
        EXX, BndsRn, NGsBlk, Gbnd
    """
    
    ax.set_title("Convergence test for %s" % conv_param,
                fontname='fantasy')
    param_split = taskids[0].split('_')
    indx = 0
    for i in range(len(param_split)):
        if conv_param in param_split[i]:
            indx = i
            break
    labels = [task.split('_')[indx] for task in taskids]
    
    for i in range(len(taskids)):
        dt = np.genfromtxt("o-"+taskids[i]+'.qp')
        # scatter plot (Note: `plt.scatter` doesn't use default colors)
        ax.plot(dt[:,2], dt[:,3], 'o', label=labels[i])
    
    ax.legend()
    ax.set_xlabel(r"$\varepsilon_\mathrm{DFT}$ (eV)",
                fontsize=16)
    ax.set_ylabel(r"$\varepsilon_\mathrm{QP}-\varepsilon_\mathrm{DFT}$ (eV)",
                fontsize=16)


if __name__ == "__main__":
    fig, ax = plt.subplots()
    plot(ax, ['HF5RY_X5Ry_nb_50_ngb_50_paral'], 'HF')
    plt.show()