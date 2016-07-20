#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wall drainage plots
@author: Pavel Ferkl
"""
from __future__ import division,unicode_literals
import numpy
from math import radians,cos,sin,tan,pi,sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.animation as animation
path='./'
plots=[1,1,0,1]
saveplots=[0,0,0,0]
fontsize=25
linewidth=4
data=numpy.loadtxt(path+'filmthickness.csv')
time,dr,np,vt,fs,hmin,hloc,hcenter,havg=\
    numpy.loadtxt(path+'results_1d.csv',skiprows=1,unpack=True)
points=int(np[0]) #discretization points
radius=[]
for i in range(points):
    radius.append(dr[0]*(0.5+i))
xran=len(radius)
radius=numpy.array(radius)
if plots[0]:
    nlines=5
    xpart=1
    xran=int(len(radius)*xpart)
    fig = plt.figure(figsize=(10.0,8.0))
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    for i in range(nlines):
        j=int((len(time)-1)/(nlines-1)*i)
        for k in range(xran):
            radius[k]=dr[j]*(0.5+k)
        line1,=plt.plot(radius[0:xran]*1e6,data[j][0:xran]*1e6,
            label='t={0:g} s'.format(time[j]),lw=linewidth)
#    plt.ylim(0,1)
#    plt.xlim(0,40)
    plt.legend(loc=2, fontsize=fontsize)
    plt.xlabel('Distance from center (μm)', fontsize=fontsize)
    plt.ylabel('Film half thickness (μm)', fontsize=fontsize)
    if saveplots[0]:
        plt.savefig('profiles.pdf')
    plt.show()

if plots[1]:
    colors=['b','g','m']
    linewidth=2
    times=[0,int(len(data)/10),len(data)-1]
    fig = plt.figure(figsize=(10.0,10.0))
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    j=numpy.argmax(dr)
    for k in range(xran):
        radius[k]=dr[j]*(0.5+k)
    Rc=radius[-1]+dr[j]/2
    plt.xlim(-Rc/sqrt(3.0)*1e6,Rc*1e6)
    plt.ylim(-Rc*1e6,Rc*1e6)
    for i,j in enumerate(times):
        xran=len(radius)
        for k in range(xran):
            radius[k]=dr[j]*(0.5+k)
        Rc=radius[-1]
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
                label='t={0:g} s'.format(time[times[i]]),lw=linewidth)
        else:
            line1,=plt.plot(x1*1e6,y1*1e6,color=colors[i],lw=linewidth)
        a=tan(-pi/3)
        c=0
        d=(x1+a*(y1-c))/(1+a**2)
        x2=2*d-x1
        y2=2*d*a-y1+2*c
        if i==1:
            line2,=plt.plot(x2*1e6,y2*1e6,color=colors[i],
                label='t={0:g} s'.format(time[times[i]]),lw=linewidth)
        else:
            line2,=plt.plot(x2*1e6,y2*1e6,color=colors[i],lw=linewidth)
        a=tan(0)
        c=0
        d=(x1+a*(y1-c))/(1+a**2)
        x3=2*d-x1
        y3=2*d*a-y1+2*c
        if i==2:
            line3,=plt.plot(x3*1e6,y3*1e6,color=colors[i],
                label='t={0:g} s'.format(time[times[i]]),lw=linewidth)
        else:
            line3,=plt.plot(x3*1e6,y3*1e6,color=colors[i],lw=linewidth)
    plt.legend(loc=1, fontsize=fontsize)
    plt.xlabel('Dimensions in (μm)', fontsize=fontsize)
    plt.ylabel('Dimensions in (μm)', fontsize=fontsize)
    if saveplots[1]:
        plt.savefig('profiles2d.pdf')
    plt.show()

if plots[2]: #not to scale when bubble grows
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

if plots[3]:
    j=numpy.argmax(dr)
    for k in range(xran):
        radius[k]=dr[j]*(0.5+k)
    fig,ax = plt.subplots(figsize=(10.0,8.0))
    ax.set_xlim(radius[0]*1e6, (radius[-1]+dr[0]/2)*1e6)
    ax.set_ylim(0, max(map(max,data))*1e6)
    line,=ax.plot(radius[0:xran]*1e6,data[0][0:xran]*1e6,
        label='t={0:.1e} s'.format(time[0]),lw=linewidth)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.legend(loc=2, fontsize=fontsize)
    plt.xlabel('position', fontsize=fontsize)
    plt.ylabel('half thickness (um)', fontsize=fontsize)
    iii=[]
    iii.append(-1)

def update(x):
    iii[0]=iii[0]+1
    if iii[0]>len(time)-1:
        iii[0]=0
    plt.legend(['t={0:.1e} s'.format(time[iii[0]])], loc=2, fontsize=fontsize)
    for k in range(xran):
        radius[k]=dr[iii[0]]*(0.5+k)
    line.set_xdata(radius[0:xran]*1e6)
    line.set_ydata(data[iii[0]][0:xran]*1e6)
    return line,

def data_gen():
    while True: yield data[0][0:xran]*1e6

if plots[3]:
    try:
        ani=animation.FuncAnimation(fig,update,data_gen,interval=100)
        plt.show()
    except:
        pass
    if saveplots[3]:
        print 'saving of animation not implemented'
