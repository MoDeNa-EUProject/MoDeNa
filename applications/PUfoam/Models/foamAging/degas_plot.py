# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:44:16 2015

@author: Pavel Ferkl
"""
from __future__ import division
from numpy import loadtxt
import matplotlib.pyplot as plt
directory='/home/me/workspace/degas/'
infile=open(directory+'input.par','r')
l=int(infile.readline().split()[0])
infile.close()
plt.figure(1).clf()
plt.figure(2).clf()
line=[]
for i in range(1,l+1,l-1):
    j='{0:04d}'.format(i)
    infile=open(directory+'H2perm_'+j+'.dat','r')
    time,x,cair,ccd=loadtxt(infile,skiprows=0,unpack=True)
    plt.figure(1)
    line1,=plt.plot(x,ccd,'o',label='t={0:.1f} days'.format(time[0]))
    line.append(line1)
    plt.show()
    plt.figure(2)
    line1,=plt.plot(x,cair,'o',label='t={0:.1f} days'.format(time[0]))
    line.append(line1)
    plt.show()
    infile.close()

plt.figure(1)
plt.legend()
plt.title('Carbon dioxide')
plt.xlabel('Position')
plt.ylabel('Concentration')
plt.figure(2)
plt.legend()
plt.title('Air')
plt.xlabel('Position')
plt.ylabel('Concentration')