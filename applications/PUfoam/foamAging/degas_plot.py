# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:44:16 2015

@author: Pavel Ferkl
"""
from __future__ import division
from numpy import loadtxt
import matplotlib.pyplot as plt
#directory_in='/home/me/workspace/degas/'
#directory_out='/home/me/workspace/degas/'
directory_in='./'
directory_out='results/'
infile=open(directory_in+'input.in','r')
l=int(infile.readline().split()[0])
infile.close()
plt.figure(1).clf()
plt.figure(2).clf()
plt.figure(3).clf()
line=[]
for i in range(1,l+1,l-1):
    j='{0:04d}'.format(i)
    infile=open(directory_out+'H2perm_'+j+'.dat','r')
    time,x,cair,ccd,ccyp=loadtxt(infile,skiprows=0,unpack=True)
    plt.figure(1)
    line1,=plt.plot(x,ccd,'o',label='t={0:.1f} days'.format(time[0]))
    line.append(line1)
    plt.show()
    plt.figure(2)
    line1,=plt.plot(x,cair,'o',label='t={0:.1f} days'.format(time[0]))
    line.append(line1)
    plt.show()
    plt.figure(3)
    line1,=plt.plot(x,ccyp,'o',label='t={0:.1f} days'.format(time[0]))
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
plt.figure(3)
plt.legend()
plt.title('Cyclopentane')
plt.xlabel('Position')
plt.ylabel('Concentration')
