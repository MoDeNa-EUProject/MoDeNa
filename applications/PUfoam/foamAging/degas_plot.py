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
def main():
    directory_in='./'
    directory_out='results/'
    infile=open(directory_in+'foamAging.json','r')
    decoded=json.load(infile)
    l=decoded['numerics']['numberOfOutputs']
    infile.close()
    plt.figure(1).clf()
    plt.figure(2).clf()
    plt.figure(3).clf()
    line=[]
    lw=2
    for i in arange(0,l+1,int((l-1)/6)):
        j='{0:04d}'.format(i)
        # infile=open(directory_out+'H2perm_'+j+'.dat','r')
        infile=open(directory_out+'ppar_'+j+'.dat','r')
        time,x,cair,ccd,ccyp=loadtxt(infile,skiprows=0,unpack=True)
        plt.figure(1)
        line1,=plt.plot(x,ccd,lw=lw,label='t={0:.1f} days'.format(time[0]))
        line.append(line1)
        plt.figure(2)
        line1,=plt.plot(x,cair,lw=lw,label='t={0:.1f} days'.format(time[0]))
        line.append(line1)
        plt.figure(3)
        line1,=plt.plot(x,ccyp,lw=lw,label='t={0:.1f} days'.format(time[0]))
        line.append(line1)
        infile.close()
    plt.figure(1)
    plt.legend()
    plt.title('Carbon dioxide')
    plt.xlabel('Position')
    plt.ylabel('Partial pressure / Pa')
    plt.figure(2)
    plt.legend()
    plt.title('Air')
    plt.xlabel('Position')
    plt.ylabel('Partial pressure / Pa')
    plt.figure(3)
    plt.legend()
    plt.title('Cyclopentane')
    plt.xlabel('Position')
    plt.ylabel('Partial pressure / Pa')
    plt.show()

if __name__ == '__main__':
    main()
