#!/usr/bin/env python
from __future__ import division
from numpy import loadtxt,pi
with open("M0.txt") as fl:
    t,m0=loadtxt(fl,unpack=True)
with open("M1.txt") as fl:
    t,m1=loadtxt(fl,unpack=True)
with open("M2.txt") as fl:
    t,m2=loadtxt(fl,unpack=True)
with open("rho_foam.txt") as fl:
    t,rho_foam=loadtxt(fl,unpack=True)
with open("wBA_g.txt") as fl:
    t,wba=loadtxt(fl,unpack=True)
with open("wCO2_g.txt") as fl:
    t,wco2=loadtxt(fl,unpack=True)
rb=(3*m1/(4*pi*m0))**(1/3)
# var=(3/(4*pi)*(m2/m0-(m1/m0)**2))**(1/3) # not working, negative values

# TODO change weight fractions to partial pressures
with open("after_foaming.txt","w") as fl:
    fl.write("{0:16s} {1:16s} {2:16s} {3:16s}\n".format(
    "#foam_density","mean_cell_diam.","pressure1","pressure2"))
    fl.write("{0:16.3e} {1:16.3e} {2:16.3e} {3:16.3e}".format(
    rho_foam[-1],rb[-1],wba[-1],wco2[-1]))
