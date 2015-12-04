"""
@author: Pavel Ferkl
"""
import periodicBox
filenameIn = "PeriodicRVE.stl"
filenameOut = "PeriodicRVEDomain.stl"
# Define domain origin and size
xmin=4
ymin=4
zmin=4
dx=4
dy=4
dz=4
render=0
periodicBox.main(filenameIn,filenameOut,xmin,ymin,zmin,dx,dy,dz,render)
