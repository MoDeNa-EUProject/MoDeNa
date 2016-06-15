# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:04:39 2015

@author: Pavel Ferkl
"""
from __future__ import division
import numpy
from math import radians,cos,sin,tan,pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
path='./'
plots=[1,1,0]
saveplots=[0,0,0]
data=numpy.loadtxt(path+'filmthickness.csv')
time,dr,np,vf,vs,vt=numpy.loadtxt(path+'results_1d.csv',skiprows=1,unpack=True)
points=int(np[0]) #discretization points
radius=[]
for i in range(points):
    radius.append(dr[0]*(0.5+i))
xran=len(radius)
radius=numpy.array(radius)
if plots[0]:
    nlines=5
    xpart=1
    xran=len(radius)*xpart
    fig = plt.figure(figsize=(5.0,4.0))
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    for i in range(nlines):
        j=(len(time)-1)/(nlines-1)*i
        for k in range(xran):
            radius[k]=dr[j]*(0.5+k)
        line1,=plt.plot(radius[0:xran]*1e6,data[j][0:xran]*1e6,
            label='t={0:.1e} s'.format(time[j]),lw=2)
#    plt.ylim(0,1)
#    plt.xlim(0,40)
    plt.legend(loc=2, fontsize=14)
    plt.xlabel('position', fontsize=14)
    plt.ylabel('half thickness (um)', fontsize=14)
    plt.show()
    if saveplots[0]:
        plt.savefig('profiles.pdf')

if plots[1]:
    colors=['b','g','m']
    times=[0,len(data)/20,len(data)-1]
    fig = plt.figure(figsize=(5.0,5.0))
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    Rc=radius[-1]
    for i,j in enumerate(times):
        x0=radius[0:xran]
        y0=data[times[i]][0:xran]
        a=tan(pi/2)
        c=0
        d=(x0+a*(y0-c))/(1+a**2)
        x1=2*d-x0+Rc+data[times[i]][-1]*tan(pi/6)
        y1=2*d*a-y0+2*c
        a=tan(pi/3)
        c=0
        d=(x1+a*(y1-c))/(1+a**2)
        x2=2*d-x1
        y2=2*d*a-y1+2*c
        x2=x2[::-1]
        y2=y2[::-1]
        x1=numpy.append(x1,x2)
        y1=numpy.append(y1,y2)
        if i==0:
            line1,=plt.plot(x1*1e6,y1*1e6,color=colors[i],
                label='t={0:.1e} s'.format(time[times[i]]))
        else:
            line1,=plt.plot(x1*1e6,y1*1e6,color=colors[i])
        a=tan(-pi/3)
        c=0
        d=(x1+a*(y1-c))/(1+a**2)
        x2=2*d-x1
        y2=2*d*a-y1+2*c
        if i==1:
            line2,=plt.plot(x2*1e6,y2*1e6,color=colors[i],
                label='t={0:.1e} s'.format(time[times[i]]))
        else:
            line2,=plt.plot(x2*1e6,y2*1e6,color=colors[i])
        a=tan(0)
        c=0
        d=(x1+a*(y1-c))/(1+a**2)
        x3=2*d-x1
        y3=2*d*a-y1+2*c
        if i==2:
            line3,=plt.plot(x3*1e6,y3*1e6,color=colors[i],
                label='t={0:.1e} s'.format(time[times[i]]))
        else:
            line3,=plt.plot(x3*1e6,y3*1e6,color=colors[i])
    plt.legend()
    plt.show()
    if saveplots[1]:
        plt.savefig('profiles2d.pdf')

if plots[2]:
    fig = plt.figure(figsize=(5.0,4.0))
    ax = fig.gca(projection='3d')
    radius3d, time3d = numpy.meshgrid(radius, time)
    surf = ax.plot_surface(time3d, radius3d*1e6, data*1e6, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #surf = ax.plot_surface(time3d, radius3d, data, rstride=2, cstride=2, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # ax.set_xlabel('time')
    # ax.set_ylabel('radius (um)')
    # ax.set_zlabel('thickness (um)')
    # plt.rc('xtick', labelsize=10)
    # plt.rc('ytick', labelsize=10)
    plt.show()
    if saveplots[2]:
        plt.savefig('profiles3d.png')
