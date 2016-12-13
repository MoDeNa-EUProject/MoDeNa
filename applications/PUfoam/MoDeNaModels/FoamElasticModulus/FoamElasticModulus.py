__author__ = 'IMDEA Materials-Mohammad Marvi Mashhadi: mohammad.marvi@imdea.org'
import os
import os.path
from itertools import groupby
import numpy as np
import math
from scipy.integrate import quad
import time
import commands
import glob
import random
# Main Inputs:
MU=0.39252             #Mean value of cell size distribution
SIGMA=0.1568           #Standard deviation of cell size distribution
DensityOfSolidPU=1200  #Real density of Solid Polyurethan (kg/m3) 
E='2.4e9'              #Elastic modulus of Solid Polyurethan (Pa)
####################################################################
NumOfCells=50
MeshSize=0.1
AppliedStrain=-0.01
DensityOfPUFoam=30.20  #Real density of foam (kg/m3)
PUinWall=0.70          #Amount of material in walls (in percentage %)
############################################################
#Material Elasticity properties:
G='0.9e9'
NU='0.35'
Density=DensityOfSolidPU #Artificial Density  
#############################################################
#Cross_sectional area function constants (Polynomial_power 4)
A1=2.77
A2=0.962
A3=3.4033
A4=(-0.2291)
A5=1.0255
############################################################
#####To create input file for SpherePack:
SuccessInMeshing=1
while SuccessInMeshing==1:
    SuccessInTessellation=1
    while SuccessInTessellation==1:
        NumSpheres=NumOfCells
        Rad=abs((np.random.normal(MU, SIGMA,NumSpheres)))
        Rads1=list(range(NumSpheres))
        t=0
        for i in range(NumSpheres):
            c=abs(Rad[t])
            Rads1[t]=c.astype(np.float)
            t=t+1
        Rads1=sorted(Rads1)
        v=0.00
        for i in range(NumSpheres):
            v=v+((2.00*Rads1[i])**3.00)
        centers = [[0 for i in range(3)] for j in range(NumSpheres)]
        v=v*1.40
        lc=v**(1.00/3.00)
        K=0
        while K==0:
            j=-1
            h=0
            timeout = time.time() + 10
            while NumSpheres>=j and h==0:
                if time.time()>timeout:
                    h=1
                    break
                j=j+1
                if j==NumSpheres:
                    K=1
                    break
                PickCenterX,PickCenterY,PickCenterZ=lc*random.random(),lc*random.random(),lc*random.random()
                while (lc-Rads1[j]>=PickCenterX and lc-Rads1[j]>=PickCenterY and lc-Rads1[j]>=PickCenterZ and Rads1[j]<PickCenterX and Rads1[j]<PickCenterY and Rads1[j]<PickCenterZ):
                    PickCenterX,PickCenterY,PickCenterZ=lc*random.random(),lc*random.random(),lc*random.random()
                centers[j][0],centers[j][1],centers[j][2]=PickCenterX,PickCenterY,PickCenterZ
                KeepCentreX, KeepCentreY, KeepCentreZ, KeepR=PickCenterX, PickCenterY, PickCenterZ, Rads1[j]
                if j>0:
                    for t in range(0,j):
                        if (((( ((KeepCentreX-centers[t][0])**2.00)+((KeepCentreY-centers[t][1])**2.00)+((KeepCentreZ-centers[t][2])**2.00))**0.50)-(KeepR+Rads1[t]))<0.000) and t!=j:
                            centers[j][0],centers[j][0],centers[j][0]=0,0,0
                            j=j-1
                            break
        mypath=os.getcwd()
        myfile=os.path.join(mypath,'Project01.rco')
        f=open(myfile,'w')
        for i in range(NumSpheres):
            f.write('{0:f}{1}{2:f}{3}{4:f}{5}{6:f}\n'.format(centers[i][0],'	',centers[i][1],'	',centers[i][2],'	',2.0*Rads1[i]))
        f.close()
        MAXcenters=max(centers)
        Mincenters=min(centers)
        EdgeCubeSize=[math.ceil(MAXcenters[0]-Mincenters[0]),math.ceil(MAXcenters[1]-Mincenters[1]),math.ceil(MAXcenters[2]-Mincenters[2])]
        EdgeRVESize=int(3.5*max(EdgeCubeSize)) #For NEPER: Size of edge of RVE
        ###########################################################
        os.chdir(mypath)
        while True:
            time.sleep(1)
            files = os.listdir('.')
            if ('Project01.rco' in files):
                break
        myfile12=os.path.join(mypath,'Project01.rco')
        CentersRads=np.loadtxt(myfile12,usecols = (0,1,2,3))
        os.remove('Project01.rco')
        Centers=CentersRads[:,:3]  #All centers of spheres
        Rads=CentersRads[:,3]      #All radii of spheres
        Rads=Rads/2
        MaxCenters=np.amax(Centers, axis=0)
        L2=int(0.5+MaxCenters[0])
        X=L2+Centers[:,0]
        Y=L2+Centers[:,1]
        Z=L2+Centers[:,2]
        Centers27=np.array([X,X+L2,X-L2,X+L2,X,X+L2,X,X-L2,X-L2,X-L2,X+L2,X,X,X-L2,X+L2,X+L2,X,X,X-L2,X,X,X-L2,X+L2,X+L2,X-L2,X-L2,X+L2,Y,Y+L2,Y-L2,Y+L2,Y+L2,Y,Y-L2,Y,Y-L2,Y+L2,Y-L2,Y-L2,Y+L2,Y,Y,Y,Y+L2,Y,Y,Y-L2,Y,Y+L2,Y-L2,Y+L2,Y-L2,Y+L2,Y-L2,Z,Z+L2,Z-L2,Z,Z+L2,Z+L2,Z-L2,Z-L2,Z,Z,Z,Z+L2,Z-L2,Z+L2,Z-L2,Z,Z,Z+L2,Z,Z,Z-L2,Z+L2,Z+L2,Z-L2,Z+L2,Z-L2,Z-L2])
        # Creation of input .txt files for neper
        myfile3=os.path.join(mypath,'Centers.txt')
        ff=open(myfile3,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(Centers27[i,j],Centers27[i+27,j],Centers27[i+54,j]))
        ff.close()
        myfile4=os.path.join(mypath,'Rads.txt')
        fff=open(myfile4,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                fff.write('{0:f}\n'.format(Rads[j]))
        fff.close()
        commandTessellation="neper -T -n {0:d} -domain 'cube({1:d},{2:d},{3:d})' -morpho @Centers.txt -weight @Rads.txt -regularization 1 -mloop 15 -o RVE27 -format geo -statcell vol -statedge length -statface area -statver x".format((27*NumSpheres),EdgeRVESize,EdgeRVESize,EdgeRVESize)
        os.system(commandTessellation)
        time.sleep(5)
        files = os.listdir('.')
        if ('RVE27.geo' in files):
            SuccessInTessellation=0
        os.remove('Rads.txt')
        os.remove('Centers.txt')
    ################################################################
    ######Extraction of middle Representative volume element########
    ################################################################
    myfile12=os.path.join(mypath,'RVE27.stver')
    verX=np.loadtxt(myfile12)
    NumOfNodes=len(verX)
    myfile12=os.path.join(mypath,'RVE27.stedge')
    EdgesLength=np.loadtxt(myfile12)
    NumOfEdges=len(EdgesLength)
    myfile12=os.path.join(mypath,'RVE27.stface')
    SurfacesArea=np.loadtxt(myfile12)
    NumOfSurfaces=len(SurfacesArea)
    myfile12=os.path.join(mypath,'RVE27.stcell')
    CellVolumes=np.loadtxt(myfile12)
    NumOfVolumes=len(CellVolumes)
    ####   Reading the GEO file
    myfile12=os.path.join(mypath,'RVE27.geo')
    text_file = open(myfile12, "r")
    lines = text_file.readlines()
    ####   Extraction of Nodes
    t=0
    Nodes=list(range(NumOfNodes))
    for i in range(0,NumOfNodes):
        currentline = lines[i].split("{")
        a=currentline[1].split("}")
        b=a[0].split(",")
        c=np.array(b)
        Nodes[t]=np.absolute(c.astype(np.float))
        t=t+1
    #####   Extraction of Edges
    Edges=list(range(NumOfEdges))
    r=0
    t=0
    for i in range(NumOfNodes,NumOfNodes+NumOfEdges):
        currentline = lines[i].split("{")
        a=currentline[1].split("}")
        b=a[0].split(",")
        c=np.array(b)
        Edges[t]=np.absolute(c.astype(np.float))
        t=t+1
    ####   Extraction of Faces
    Faces=list(range(NumOfSurfaces))
    r=0
    t=0
    for i in range(0,NumOfSurfaces):
        j=NumOfNodes+NumOfEdges+(2*i)
        currentline = lines[j].split("{")
        a=currentline[1].split("}")
        b=a[0].split(",")
        c=np.array(b)
        Faces[t]=np.absolute(c.astype(np.float))
        t=t+1
    ####   Extraction of volumes
    Volumes=list(range(NumOfVolumes))
    a=list(lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)])
    a[0:21]=''
    o=len(a)
    a[o-3:]=''
    lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)]="".join(a)
    currentline=lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)].split(",")
    c=np.array(currentline)
    Volumes[0]=np.absolute(c.astype(np.float))
    t=1
    for i in range(1,NumOfVolumes):
        j=NumOfNodes+NumOfEdges+(2*NumOfSurfaces)+i
        currentline = lines[j].split(" = {")
        a=currentline[2].split("}")
        b=a[0].split(",")
        c=np.array(b)
        Volumes[t]=np.absolute(c.astype(np.float))
        t=t+1
    #####################################################
    MAX0=list(range(0,NumOfCells))
    for i in range(NumOfCells):
        MAX0[i]=max(Volumes[i])
    MaxIndexOfFaces=max(MAX0)
    MAX1=list(range(0,int(MaxIndexOfFaces)))
    for i in range(int(MaxIndexOfFaces)):
        MAX1[i]=max(Faces[i])
    MaxIndexOfEdges=max(MAX1)
    MAX2=list(range(0,int(MaxIndexOfEdges)))
    for i in range(int(MaxIndexOfEdges)):
        MAX2[i]=max(Edges[i])
    MaxIndexOfNodes=max(MAX2)
    text_file.close()
    ####################################################
    # Making GEO file containing Periodic RVE
    myfile13=os.path.join(mypath,'PeriodicRVE.geo')
    text_GEO = open(myfile13, "w")
    text_GEO.write('{0}{1}{2}\n'.format('cl=',MeshSize,';'))
    for i in range(int(MaxIndexOfNodes)):
        lines[i] = lines[i].replace('};',',cl};')
        text_GEO.write('{0}'.format(lines[i]))
    for i in range(int(NumOfNodes),int(NumOfNodes+MaxIndexOfEdges)):
        text_GEO.write('{0}'.format(lines[i]))
    j=int(NumOfNodes+NumOfEdges)
    for i in range(int(NumOfNodes+NumOfEdges),int(NumOfNodes+NumOfEdges+MaxIndexOfFaces)):
        text_GEO.write('{0}{1}'.format(lines[j],lines[j+1]))
        j=2+j
    text_GEO.write('{0}{1}{2}\n'.format(' Surface Loop (1) = {',lines[int(NumOfNodes+NumOfEdges+(2*NumOfSurfaces))],'};'))
    NumOfCells=NumOfVolumes/27
    j=int(NumOfNodes+NumOfEdges+(2*NumOfSurfaces)+1)
    for i in range(int(NumOfCells-1)):
        text_GEO.write('{0}'.format(lines[j]))
        j+=1
    text_GEO.write('{0}{1}{2}{3}{4}'.format('Volume (',NumOfCells,') = {',NumOfCells,'};'))
    text_GEO.close()
    ################################ End of Geometry construction##################################
    os.remove('RVE27.geo')
    os.remove('RVE27.stcell')
    os.remove('RVE27.stedge')
    os.remove('RVE27.stface')
    os.remove('RVE27.stver')
    ######################################################################
                      # Creating Mesh Using GMSH
    ######################################################################
    commandMeshing="gmsh PeriodicRVE.geo -2 -optimize_lloyd -o PeriodicRVE -format inp -saveall"
    os.system(commandMeshing)
    time.sleep(15)
    os.remove('PeriodicRVE.geo')
    files = os.listdir('.')
    if ('PeriodicRVE.inp' in files):
        SuccessInMeshing=0
        break
######################################################################
     # Creating Corrected INP To Be processed In Next Step
######################################################################
INPgeom = open('PeriodicRVE.inp', 'r')
lines = INPgeom.readlines()
INPgeom.close()
os.remove('PeriodicRVE.inp')
text_inp = open('CorrectedINPafterGMSH.inp', "w")
text_inp.write('{0}'.format(lines[2]))
i=3
while lines[i][0]!='*':
    text_inp.write('{0}'.format(lines[i]))
    i=i+1
text_inp.write('{0}\n'.format('*ELEMENT, TYPE=B31, ELSET=auto1'))
K=0
while K==0:
    if 30>len(lines[i]):
        text_inp.write('{0}'.format(lines[i]))
    if len(lines[i])>30:
        if lines[i][15]=='C':
            K=1
    i=i+1
text_inp.write('{0}\n'.format('*ELEMENT, TYPE=S3R, ELSET=auto1'))
K=0
while K==0:
    if ((lines[i][0]!='*') and (lines[i][3]!='S')):
        text_inp.write('{0}'.format(lines[i]))
    if lines[i][3]=='S':
        K=1
    i=i+1
text_inp.write('{0}'.format('*******'))
text_inp.close()
###########################################
myfile12=os.path.join(mypath,'CorrectedINPafterGMSH.inp')
INPgeom = open(myfile12, "r")
lines = INPgeom.readlines()
INPgeom.close()
os.remove('CorrectedINPafterGMSH.inp')
#Extraction of Nodes:
i=1
r=lines[1][0]
Nodes=np.array([])
while r!='*':
    b = lines[i].split(",")
    c=np.array(b)
    Node=np.absolute(c.astype(np.float))
    Nodes=np.append(Nodes,Node)
    i=i+1
    r=lines[i][0]
i=i+1
r=lines[i][0]
Beams=np.array([])
while r!='*':
    b = lines[i].split(",")
    c=np.array(b)
    Beam=np.absolute(c.astype(np.float))
    Beams=np.append(Beams,Beam)
    i=i+1
    r=lines[i][0]
i=i+1
r=lines[i][0]
Shells=np.array([])
while r!='*':
    b = lines[i].split(",")
    c=np.array(b)
    Shell=np.absolute(c.astype(np.float))
    Shells=np.append(Shells,Shell)
    i=i+1
    r=lines[i][0]
##########################################################################
#Calculation of Total area of shells
k=1
AreaTotalShells=0
for i in range(int(len(Shells)/4)):
    x1=Nodes[1+int(4*(Shells[k]-1))]
    y1=Nodes[2+int(4*(Shells[k]-1))]
    z1=Nodes[3+int(4*(Shells[k]-1))]
    x2=Nodes[1+int(4*(Shells[k+1]-1))]
    y2=Nodes[2+int(4*(Shells[k+1]-1))]
    z2=Nodes[3+int(4*(Shells[k+1]-1))]
    x3=Nodes[1+int(4*(Shells[k+2]-1))]
    y3=Nodes[2+int(4*(Shells[k+2]-1))]
    z3=Nodes[3+int(4*(Shells[k+2]-1))]
    L12=(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))**0.5
    L13=(((x1-x3)**2)+((y1-y3)**2)+((z1-z3)**2))**0.5
    L32=(((x3-x2)**2)+((y3-y2)**2)+((z3-z2)**2))**0.5
    s=(L12+L13+L32)/2
    AreaTotalShells=((s*(s-L32)*(s-L13)*(s-L12))**0.5)+AreaTotalShells
    k=k+4
#Calculation of Total length of Beam Elements:
BeamLengths=np.array([])
k=1
#LengthTotalBeams=0
for i in range(int(len(Beams)/3)):
   x1=Nodes[1+int(4*(Beams[k]-1))]
   y1=Nodes[2+int(4*(Beams[k]-1))]
   z1=Nodes[3+int(4*(Beams[k]-1))]
   x2=Nodes[1+int(4*(Beams[k+1]-1))]
   y2=Nodes[2+int(4*(Beams[k+1]-1))]
   z2=Nodes[3+int(4*(Beams[k+1]-1))]
   BeamLengths=np.append(BeamLengths,((((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))**0.5))
   k=k+3
#####################################################################################################
###Number of beams per edge
counter=1
BeamElSlopes=np.array([])
for i in range(int((len(Beams))/3)):
    X1=Nodes[1+(4*((Beams[1+(3*(counter-1))])-1))]
    X2=Nodes[1+(4*((Beams[2+(3*(counter-1))])-1))]
    Y1=Nodes[2+(4*((Beams[1+(3*(counter-1))])-1))]
    Y2=Nodes[2+(4*((Beams[2+(3*(counter-1))])-1))]
    Z1=Nodes[3+(4*((Beams[1+(3*(counter-1))])-1))]
    Z2=Nodes[3+(4*((Beams[2+(3*(counter-1))])-1))]
    ElNumber=int(Beams[0+(3*(counter-1))])
    d=(((X1-X2)**2.0)+((Y1-Y2)**2.0)+((Z1-Z2)**2.0))**0.5
    Xslope=(X1-X2)/d
    Yslope=(Y1-Y2)/d
    Zslope=(Z1-Z2)/d
    BeamElSlopes=np.append(BeamElSlopes,int(1000*(Zslope+Xslope+Yslope)))
    counter=counter+1
ElPerBeamList=[len(list(group)) for key, group in groupby(BeamElSlopes)]
NoOfBeamPerEdge=np.asarray(ElPerBeamList)
################Extraction of Node pairs for priodic BC######################
BeamElementsNodes1=np.array([])
for i in range(len(Beams)):
    if i%3!=0:
        BeamElementsNodes1=np.append(BeamElementsNodes1,Beams[i])
BeamElementsNodes=np.unique(BeamElementsNodes1)
MaxBeamElementsNodes=int(max(BeamElementsNodes))
NodePairsXY=np.array([])
m=0
n=0
for i in range(MaxBeamElementsNodes):
    for j in range(i,MaxBeamElementsNodes):
        if Nodes[m]!=Nodes[n] and Nodes[m+1]==Nodes[n+1] and Nodes[m+2]==Nodes[n+2]:
            NodePairsXY=np.append(NodePairsXY,Nodes[m])
            NodePairsXY=np.append(NodePairsXY,Nodes[n])
        n=n+4
    m=m+4
    n=m
NodePairsXZ=np.array([])
m=0
n=0
for i in range(MaxBeamElementsNodes):
    for j in range(i,MaxBeamElementsNodes):
        if Nodes[m]!=Nodes[n] and Nodes[m+1]==Nodes[n+1] and Nodes[m+3]==Nodes[n+3]:
            NodePairsXZ=np.append(NodePairsXZ,Nodes[m])
            NodePairsXZ=np.append(NodePairsXZ,Nodes[n])
        n=n+4
    m=m+4
    n=m
NodePairsYZ=np.array([])
m=0
n=0
for i in range(MaxBeamElementsNodes):
    for j in range(i,MaxBeamElementsNodes):
        if Nodes[m]!=Nodes[n] and Nodes[m+2]==Nodes[n+2] and Nodes[m+3]==Nodes[n+3]:
            NodePairsYZ=np.append(NodePairsYZ,Nodes[m])
            NodePairsYZ=np.append(NodePairsYZ,Nodes[n])
        n=n+4
    m=m+4
    n=m
######################################################################
NodeCordinates=np.array([])
K=0
for i in range(MaxBeamElementsNodes):
    NodeCordinates=np.append(NodeCordinates,Nodes[K+1])
    NodeCordinates=np.append(NodeCordinates,Nodes[K+2])
    NodeCordinates=np.append(NodeCordinates,Nodes[K+3])
    K=K+4
UpLevelOfRVE=1+int(max(NodeCordinates))
######################################################################
#Corner Nodes:
DistCenter=np.array([])
DistX=np.array([])
DistY=np.array([])
DistZ=np.array([])
K=0
for i in range(MaxBeamElementsNodes):
    DistCenter=np.append(DistCenter,(Nodes[K+1]**2)+(Nodes[K+2]**2)+(Nodes[K+3]**2))
    DistX=np.append(DistX,((Nodes[K+1]-UpLevelOfRVE)**2)+((Nodes[K+2])**2)+((Nodes[K+3])**2))
    DistY=np.append(DistY,((Nodes[K+1])**2)+((Nodes[K+2]-UpLevelOfRVE)**2)+((Nodes[K+3])**2))
    DistZ=np.append(DistZ,((Nodes[K+1])**2)+((Nodes[K+2])**2)+((Nodes[K+3]-UpLevelOfRVE)**2))
    K=K+4
CornerNode=1+DistCenter.argmin() #Node index that should be written in inp file
CornerNodeX=1+DistX.argmin()
CornerNodeY=1+DistY.argmin()
CornerNodeZ=1+DistZ.argmin()
# Recalculation of CornerNodeX, CornerNodeY and CornerNodeZ based on CornerNode 
m=CornerNode-1
n=0
for i in range(MaxBeamElementsNodes):
    if Nodes[4*m]!=Nodes[4*n] and Nodes[(4*m)+2]==Nodes[(4*n)+2] and Nodes[(4*m)+3]==Nodes[(4*n)+3]:
        CornerNodeX=n+1
    n=n+1

n=0
for i in range(MaxBeamElementsNodes):
    if Nodes[4*m]!=Nodes[4*n] and Nodes[(4*m)+1]==Nodes[(4*n)+1] and Nodes[(4*m)+3]==Nodes[(4*n)+3]:
        CornerNodeY=n+1
    n=n+1
n=0
for i in range(MaxBeamElementsNodes):
    if Nodes[4*m]!=Nodes[4*n] and Nodes[(4*m)+1]==Nodes[(4*n)+1] and Nodes[(4*m)+2]==Nodes[(4*n)+2]:
        CornerNodeZ=n+1
    n=n+1
###Removal of repeated PBC Nodes:###
LenNodePairsXY=len(NodePairsXY)
LenNodePairsXZ=len(NodePairsXZ)
LenNodePairsYZ=len(NodePairsYZ)
K=0
L=0
for i in range(int(LenNodePairsXY)):
    for j in range(int(LenNodePairsXZ/2)):
        if NodePairsXY[K]==NodePairsXZ[L] or NodePairsXY[K]==NodePairsXZ[L+1]:
            NodePairsXZ[L]=123456
            NodePairsXZ[L+1]=123456
        L=L+2
    L=0
    K=K+1
K=0
L=0
for i in range(int(LenNodePairsXY)):
    for j in range(int(LenNodePairsYZ/2)):
        if NodePairsXY[K]==NodePairsYZ[L] or NodePairsXY[K]==NodePairsYZ[L+1]:
            NodePairsYZ[L]=123456
            NodePairsYZ[L+1]=123456
        L=L+2
    L=0
    K=K+1
K=0
L=0
for i in range(int(LenNodePairsXZ)):
    for j in range(int(LenNodePairsYZ/2)):
        if NodePairsXZ[K]==NodePairsYZ[L] or NodePairsXZ[K]==NodePairsYZ[L+1]:
            NodePairsYZ[L]=123456
            NodePairsYZ[L+1]=123456
        L=L+2
    L=0
    K=K+1
L=0
for i in range(int(LenNodePairsXY/2)):
    if CornerNode==NodePairsXY[L] or CornerNodeX==NodePairsXY[L] or CornerNodeY==NodePairsXY[L] or CornerNodeZ==NodePairsXY[L] or CornerNode==NodePairsXY[L+1] or CornerNodeX==NodePairsXY[L+1] or CornerNodeY==NodePairsXY[L+1] or CornerNodeZ==NodePairsXY[L+1]:
        NodePairsXY[L]=123456
        NodePairsXY[L+1]=123456
    L=L+2
L=0
for i in range(int(LenNodePairsXZ/2)):
    if CornerNode==NodePairsXZ[L] or CornerNodeX==NodePairsXZ[L] or CornerNodeY==NodePairsXZ[L] or CornerNodeZ==NodePairsXZ[L] or CornerNode==NodePairsXZ[L+1] or CornerNodeX==NodePairsXZ[L+1] or CornerNodeY==NodePairsXZ[L+1] or CornerNodeZ==NodePairsXZ[L+1]:
        NodePairsXZ[L]=123456
        NodePairsXZ[L+1]=123456
    L=L+2
L=0
for i in range(int(LenNodePairsYZ/2)):
    if CornerNode==NodePairsYZ[L] or CornerNodeX==NodePairsYZ[L] or CornerNodeY==NodePairsYZ[L] or CornerNodeZ==NodePairsYZ[L] or CornerNode==NodePairsYZ[L+1] or CornerNodeX==NodePairsYZ[L+1] or CornerNodeY==NodePairsYZ[L+1] or CornerNodeZ==NodePairsYZ[L+1]:
        NodePairsYZ[L]=123456
        NodePairsYZ[L+1]=123456
    L=L+2
#############RVE Actual Size############################################:
Ylength=Nodes[(4*(CornerNodeY-1))+2]-Nodes[(4*(CornerNode-1))+2]
Xlength=Nodes[(4*(CornerNodeX-1))+1]-Nodes[(4*(CornerNode-1))+1]
Zlength=Nodes[(4*(CornerNodeZ-1))+3]-Nodes[(4*(CornerNode-1))+3]
RVEarea=Xlength*Zlength
Vfoam=RVEarea*Ylength #Should be the volume of RVE
######################################################################
    #####################################################
     # Creating Corrected INP To Be processed In ABAQUS
    #####################################################
######################################################################
#How material distributed between walls and struts (Vw , Vs) :
SolidPUmass=DensityOfPUFoam*Vfoam
VsolidPUinRVE=SolidPUmass/DensityOfSolidPU
PUinStruts=1-PUinWall
Vw=VsolidPUinRVE*PUinWall       #Volume of solid PU in walls              
Vs=VsolidPUinRVE*PUinStruts     #Volume of solid PU in Struts            
Vstrut=Vs/len(NoOfBeamPerEdge)  #Volume of each strut
###################################################################
myfile13=os.path.join(mypath,'INP.inp')
text_inp = open(myfile13, "w")
#Introduction of PBC text file
text_inp.write('{0}\n'.format('****Mesh Size=0.1'))
text_inp.write('{0}\n'.format('*include, input=PBC.txt'))
# Nodes & Beams & shell Printing:
for i in range(0,len(lines)):
    text_inp.write('{0}'.format(lines[i]))
text_inp.write('\n'.format(lines[i]))
#############################################################
#Section assignment for Beam Elements
TotalVolAssignedBeams=0
BeamElCounter=1
for i in range(int(len(NoOfBeamPerEdge))):
    NumElinEdge=NoOfBeamPerEdge[i] #All of numbers of beam elements in edge
    TotalAreaStrut=0
    K=0
    DeltaX=1.0/NumElinEdge
    for j in range(int(NumElinEdge)):
        area= quad(lambda x: (A1*((x-0.5)**4))+(A2*((x-0.5)**3))+(A3*((x-0.5)**2))+(A4*((x-0.5)**1))+A5, K, K+DeltaX)
        TotalAreaStrut=TotalAreaStrut+area[0]
        K=K+DeltaX
    K=0.0
    for j in range(int(NumElinEdge)):
        text_inp.write('{0}{1}\n'.format('*Elset, elset=Beam',int(BeamElCounter-1)))
        text_inp.write('{0}\n'.format(int(Beams[3*(BeamElCounter-1)])))
        text_inp.write('{0}{1}{2}{3}{4}\n'.format('*Beam Section, elset=Beam',int(BeamElCounter-1),',material = Material-1 ',', temperature=GRADIENTS', ',section=TRAPEZOID'))
        Vbeam=quad(lambda x: (A1*((x-0.5)**4))+(A2*((x-0.5)**3))+(A3*((x-0.5)**2))+(A4*((x-0.5)**1))+A5, K, K+DeltaX)
        Vbeam1=(Vbeam[0]/TotalAreaStrut)*Vstrut
        Area=Vbeam1/BeamLengths[BeamElCounter-1]
        c=0.000001
        h=(Area*(3.0**0.5))**0.5
        a=h*(2.0/(3.00**0.5))
        d=h/3.00
        K=K+DeltaX
        text_inp.write('{0:.10f}{1}{2:.10f}{3}{4:.10f}{5}{6:.10f}\n'.format(a,', ',h,', ',c,', ',d))
        text_inp.write('{0}\n'.format('0.,0.,-1'))
        TotalVolAssignedBeams=TotalVolAssignedBeams+(Area*BeamLengths[BeamElCounter-1])   ####%%%%%%%%%%%%%%%%%%%%%%%%% Calculation  of Total volume in Struts
        BeamElCounter=BeamElCounter+1
############################################################
########Section assignment for Shell Elements###############
############################################################
#Calculation of shell thickness:
ShellThickness=Vw/AreaTotalShells
#Section assignment for shell elements:
t=(len(Shells)/4)%16
tt=(int((len(Shells)/4)-t))/16
N=0
M=1
for i in range(int(tt)):
    text_inp.write('{0}{1}\n'.format('*Elset, elset=EShell',M))
    for j in range(16):
        text_inp.write('{0}{1}'.format(int(Shells[N]),','))
        N=N+4
    M=M+1
    text_inp.write('\n')
if t>0:
    text_inp.write('{0}{1}\n'.format('*Elset, elset=EShell',M))
    for j in range(int(t)):
        text_inp.write('{0}{1}'.format(int(Shells[N]),','))
        N=N+4
text_inp.write('\n')
for i in range(1,M+1):
    text_inp.write('{0}{1}{2}\n'.format('*Shell Section, elset=EShell',i,', material=Material-1'))
    text_inp.write('{0}{1}\n'.format(ShellThickness,', 5'))
text_inp.write('{0}\n'.format('*Material, name=Material-1'))
text_inp.write('{0}\n'.format('*Density'))
text_inp.write('{0}\n'.format(Density))
text_inp.write('{0}\n'.format('*Elastic'))
text_inp.write('{0}{1}{2}\n'.format(E,', ',NU))
text_inp.write('{0}\n'.format('**'))
text_inp.write('{0}{1}\n'.format('**RVEarea',RVEarea))
text_inp.write('{0}{1}\n'.format('**Ylength',Ylength))
text_inp.write('{0}\n'.format('** STEP: Step-1'))
text_inp.write('{0}\n'.format('*Step, name=Step-1, nlgeom=YES, INC=1000'))
text_inp.write('{0}\n'.format('*Static, stabilize=0.3'))
text_inp.write('{0}\n'.format('0.001, 1., 1e-09, 0.01'))
text_inp.write('{0}\n'.format('** BOUNDARY CONDITIONS'))
text_inp.write('{0}\n'.format('*Nset, nset=MN0'))
text_inp.write('{0}\n'.format(CornerNode))
text_inp.write('{0}\n'.format('*Nset, nset=MNx'))
text_inp.write('{0}\n'.format(CornerNodeX))
text_inp.write('{0}\n'.format('*Nset, nset=MNy'))
text_inp.write('{0}\n'.format(CornerNodeY))
text_inp.write('{0}\n'.format('*Nset, nset=MNz'))
text_inp.write('{0}\n'.format(CornerNodeZ))
text_inp.write('{0}\n'.format('*Boundary'))
text_inp.write('{0}\n'.format('MN0, PINNED'))
text_inp.write('{0}\n'.format('*Boundary'))
text_inp.write('{0}\n'.format('MNx, 2, 2'))
text_inp.write('{0}\n'.format('MNx, 3, 3'))
text_inp.write('{0}\n'.format('*Boundary'))
text_inp.write('{0}\n'.format('MNy, 1, 1'))
text_inp.write('{0}{1}\n'.format('MNy, 2, 2, ',Ylength*(AppliedStrain)))
text_inp.write('{0}\n'.format('MNy, 3, 3'))
text_inp.write('{0}\n'.format('*Boundary'))
text_inp.write('{0}\n'.format('MNz, 1, 1'))
text_inp.write('{0}\n'.format('MNz, 2, 2'))
text_inp.write('{0}\n'.format('** OUTPUT REQUESTS'))
text_inp.write('{0}\n'.format('*Restart, write, frequency=0'))
text_inp.write('{0}\n'.format('** FIELD OUTPUT: F-Output-1'))
text_inp.write('{0}\n'.format('*Output, field'))
text_inp.write('{0}\n'.format('*Node Output'))
text_inp.write('{0}\n'.format('RF, U'))
text_inp.write('{0}\n'.format('** HISTORY OUTPUT: H-Output-1'))
text_inp.write('{0}\n'.format('*Output, history, variable=PRESELECT'))
text_inp.write('{0}\n'.format('*End Step'))
text_inp.close()
##################################################################################
# Writing PBC: 
myfile14=os.path.join(mypath,'PBC.txt')
t=open(myfile14, "w")
t.write('{0}\n'.format('*EQUATION'))
L=0
for j in range(int(LenNodePairsXY/2)):
    if NodePairsXY[L]!=123456 and NodePairsXY[L+1]!=123456:
        for k in range(3):
            t.write('{0}\n'.format('3'))
            t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsXY[L]),',', k+1,' ,1.00,	     ',int(NodePairsXY[L+1]),', ',k+1,' ,-1.00,	     ',int(CornerNodeZ),',', k+1,' ,-1.00'))
    L=L+2
L=0
for j in range(int(LenNodePairsXZ/2)):
    if NodePairsXZ[L]!=123456 and NodePairsXZ[L+1]!=123456:
        for k in range(3):
            t.write('{0}\n'.format('3'))
            t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsXZ[L]),',', k+1,' ,1.00,	     ',int(NodePairsXZ[L+1]),', ',k+1,' ,-1.00,	     ',int(CornerNodeY),',', k+1,' ,-1.00'))
    L=L+2
L=0
for j in range(int(LenNodePairsYZ/2)):
    if NodePairsYZ[L]!=123456 and NodePairsYZ[L+1]!=123456:
        for k in range(3):
            t.write('{0}\n'.format('3'))
            t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsYZ[L]),',', k+1,' ,1.00,	     ',int(NodePairsYZ[L+1]),', ',k+1,' ,-1.00,	     ',int(CornerNodeX),',', k+1,' ,-1.00'))
    L=L+2
t.close()
########################################################################################################################
########################################### Finite Element Analysis (ABAQUS) ###########################################
########################################################################################################################
# Getting ABAQUS directory on user's computer:
t=commands.getstatusoutput('find / -name abq6* 2> /dev/null')
currentline = t[1].split("\n")
for i in range(len(currentline)):
    y=currentline[i]
    if '/Commands' in y:
        f=i
        break
if 'abq' in currentline[f]:
    g=currentline[f].index('abq')
    D=currentline[f]
    A = D[g:]
    D = D[:g]#Directory of ABAQUS
os.chdir(D)
os.system('cp {0}/INP.inp {1}'.format(mypath,D)) #Transferring INP file to Abaqus directory
os.system('cp {0}/PBC.txt {1}'.format(mypath,D)) #Transferring INP file to Abaqus directory
while True:
    time.sleep(1)
    files = os.listdir('.')
    if ('INP.inp' in files) and ('PBC.txt' in files):
        break
os.chdir(mypath)
os.remove('INP.inp')
os.remove('PBC.txt')
os.chdir(D)
NumberOfCPUs=commands.getoutput("nproc")
CMDforABAQUS='cd {0} && ./{1} job=INP cpus={2}'.format(D,A,NumberOfCPUs)
os.system(CMDforABAQUS)
while True:
    time.sleep(1)
    files = os.listdir('.')
    if ('INP.lck' not in files) and ('INP.odb' in files):
        break
#Creation of RF2 and U2 report file:
myfile14=os.path.join(D,'AbaqusPythonInpRF2U2.py')
q=open(myfile14, "w")
q.write('{0}\n'.format('from abaqus import *'))
q.write('{0}\n'.format('from abaqusConstants import *'))
q.write('{0}\n'.format('from viewerModules import *'))
q.write('{0}\n'.format('from driverUtils import executeOnCaeStartup'))
q.write('{0}\n'.format('executeOnCaeStartup()'))
q.write('{0}\n'.format('o2 = session.openOdb('))
q.write('{0}\n'.format("    name='INP.odb')"))
q.write('{0}\n'.format("session.viewports['Viewport: 1'].setValues(displayedObject=o2)"))
q.write('{0}\n'.format("session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=("))
q.write('{0}\n'.format("    CONTOURS_ON_DEF, ))"))
q.write('{0}{1}{2}\n'.format("odb = session.odbs['",D,"INP.odb']"))
q.write('{0}\n'.format("session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF',"))
q.write('{0}\n'.format("    NODAL, ((COMPONENT, 'RF2'), )), ('U', NODAL, ((COMPONENT, 'U2'), )), ),"))
q.write('{0}{1}{2}\n'.format("    nodeLabels=(('PART-1-1', ('",CornerNodeY,"', )), ))"))
q.write('{0}{1}{2}\n'.format("x0 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: ",int(CornerNodeY),"']" ))
q.write('{0}{1}{2}\n'.format("x1 = session.xyDataObjects['U:U2 PI: PART-1-1 N: ",int(CornerNodeY),"']"))
q.write('{0}\n'.format("session.xyReportOptions.setValues()"))
q.write('{0}\n'.format("session.writeXYReport(fileName='XYdataRF2U2.rpt', xyData=(x0, x1))"))
q.close()
########################################################################################################################
########################################### Output Data Base (ODB of ABAQUS) ###########################################
########################################################################################################################
while True:
    time.sleep(1)
    files = os.listdir('.')
    if ('AbaqusPythonInpRF2U2.py' in files):
        break
os.system('cd {0} && ./abaqus cae noGUI=AbaqusPythonInpRF2U2.py'.format(D)) #Executing .py to get the XY data from ODB
while True:
    time.sleep(1)
    files = os.listdir('.')
    if ('XYdataRF2U2.rpt' in files):
        break
myfile14=os.path.join(D,'XYdataRF2U2.rpt')
text_file = open(myfile14, "r")
lines = text_file.readlines()
os.chdir(D)
for filename in glob.glob("INP*"):
    os.remove(filename)
os.remove('XYdataRF2U2.rpt')
os.remove('PBC.txt')
os.remove('AbaqusPythonInpRF2U2.py')
XYdata=np.array([])
b=1
t=4
while b==1:
    x=0
    currentline = lines[t].split(" ")
    if len(currentline)==1 and currentline[0]=='\n':
        b=2
    for i in range(len(currentline)):
        if currentline[i]!='' and x!=3 and currentline[i]!='\n':
            r=float(currentline[i])
            XYdata=np.append(XYdata,r)
            x=x+1
    t=t+1
Xdata=np.array([])
Ydata=np.array([])
t=0
for i in range(int(len(XYdata)/3)):
    Xdata=np.append(Xdata,XYdata[t+1])
    Ydata=np.append(Ydata,XYdata[t+2])
    t=t+3
Ydata=Ydata/Ylength
Xdata=Xdata/(RVEarea*1000000)
[ElasticModulusOfRVE,i] = np.polyfit(Ydata, Xdata, 1)
print('ElasticModulusOfFoam=')
print(ElasticModulusOfRVE)