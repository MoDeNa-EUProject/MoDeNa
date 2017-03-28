#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:44:16 2015

@author: Pavel Ferkl
"""
from __future__ import division
from numpy import loadtxt,arange
import matplotlib.pyplot as plt
import json
import os
import pandas as pd
def main():
    nplots=4
    titles=['Oxygen','Nitrogen','Carbon dioxide','Cyclopentane']
    directory_in='./inputs/'
    directory_out='./results/foamAging/'
    infile=open(directory_in+'foamAging.json','r')
    decoded=json.load(infile)
    l=decoded['numerics']['numberOfOutputs']
    infile.close()
    for i in range(nplots):
        plt.figure(i).clf()
    line=[]
    lw=2
    for i in arange(0,l+1,int((l-1)/6)):
        j='{0:04d}'.format(i)
        with open(os.path.join(directory_out,'pres_'+j+'.csv'),'r') as infile:
            df = pd.read_csv(infile)
            time = df['time']
            x = df['position']
            co2 = df['O2']
            cn2 = df['N2']
            ccd = df['CO2']
            ccyp = df['CP']
        conc=[co2,cn2,ccd,ccyp]
        for j in range(nplots):
            plt.figure(j)
            line1,=plt.plot(
                x,conc[j],lw=lw,label='t={0:.1f} days'.format(time[0]))
            line.append(line1)
    for i in range(nplots):
        plt.figure(i)
        plt.legend()
        plt.title(titles[i])
        plt.xlabel('Position')
        plt.ylabel('Partial pressure / Pa')
    plt.show()

if __name__ == '__main__':
    main()
