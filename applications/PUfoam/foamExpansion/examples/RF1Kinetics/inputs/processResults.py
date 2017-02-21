#!/usr/bin/env python
from __future__ import division
from numpy import loadtxt,pi,average
with open("M0.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
m0=average(data[1:])
with open("M1.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
m1=average(data[1:])
with open("M2.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
m2=average(data[1:])
with open("rho_foam.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
rho_foam=average(data[1:])
with open("wBA_g.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
wba=average(data[1:])
with open("wCO2_g.txt") as fl:
    data=loadtxt(fl)
data=data[-1]
wco2=average(data[1:])
rb=(3*m1/(4*pi*m0))**(1/3)

with open("../../../constant/kineticsProperties") as fl:
    text=fl.read()
    for line in text.split("\n"):
        if "molecularMassCO2" in line:
            MCO2=float(line.split()[1].rstrip(";"))
        if "molecularMassBlowingAgent" in line:
            MBA=float(line.split()[1].rstrip(";"))
pCO2=(wco2/MCO2)/(wco2/MCO2 + wba/MBA)*1e5
pBA=(wba/MBA)/(wco2/MCO2 + wba/MBA)*1e5

with open("after_foaming.txt","w") as fl:
    fl.write("{0:16s} {1:16s} {2:16s} {3:16s}\n".format(
    "#foam_density","mean_cell_diam.","pressure1","pressure2"))
    fl.write("{0:16.3e} {1:16.3e} {2:16.3e} {3:16.3e}".format(
    rho_foam,rb,pBA,pCO2))
