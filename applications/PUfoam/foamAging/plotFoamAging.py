#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:44:16 2015

@author: Pavel Ferkl
"""
from __future__ import division, print_function
import os
import json
from numpy import arange
import matplotlib.pyplot as plt
import pandas as pd


def main():
    """Main function."""
    line_weight = 2
    directory_in = './inputs/'
    directory_out = './results/foamAging/'
    infile = open(directory_in + 'foamAging.json', 'r')
    inp = json.load(infile)
    titles = inp['foamCondition']['initialComposition'].keys()
    nplots = len(titles)
    if inp['numerics']['progressTime'] == 'linear':
        outs = inp['numerics']['numberOfOutputs']
    elif inp['numerics']['progressTime'] == 'logarithmic':
        outs = inp['numerics']['outputsPerOrder'] * \
            inp['numerics']['numberOfOrders']
    infile.close()
    for i in range(nplots):
        plt.figure(i).clf()
    for i in arange(0, outs + 1, int(outs / 6)):
        j = '{0:04d}'.format(i)
        with open(os.path.join(directory_out,
                               'pres_' + j + '.csv'), 'r') as infile:
            dtf = pd.read_csv(infile)
            time = dtf['time'.rjust(24)]
            x = dtf['position'.rjust(23)]
            conc = []
            for gas in titles:
                conc.append(dtf[gas.ljust(23)])
        for j in range(nplots):
            plt.figure(j)
            plt.plot(
                x, conc[j], lw=line_weight,
                label='t={0:.1f} days'.format(time[0]))
    for i in range(nplots):
        plt.figure(i)
        plt.legend()
        plt.title(titles[i])
        plt.xlabel('Position')
        plt.ylabel('Partial pressure / Pa')
    plt.show()


if __name__ == '__main__':
    main()
