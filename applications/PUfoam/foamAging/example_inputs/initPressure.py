# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 22:21:41 2015

@author: Pavel Ferkl
"""

p=101325
xcd=4e-4
xn=0.79
xo=1-xcd-xn
pcd=xcd*p
pn=p*xn
po=p*xo
print 'carbon dioxide pressure: {0:.3f}'.format(pcd)
print 'nitrogen pressure: {0:.3f}'.format(pn)
print 'oxygen pressure: {0:.3f}'.format(po)


