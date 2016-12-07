#!/usr/bin/python
'''@cond
   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o
Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.
License
    This file is part of Modena.
    Modena is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.
    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.
    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond'''

__author__='IMDEA Materials-Mohammad Marvi Mashhadi: mohammad.marvi@imdea.org, mohamadmarvi@gmail.com'

import os,sys,datetime,re,os.path,time,commands,glob,random,math,numpy as np
from itertools import groupby
from scipy.integrate import quad
from scipy.stats import itemfreq
# ------------------------ Convinience variables ---------------------------- #
INPUT_FILE = "Inputs.in"                                # Name of the input file
OUTPUT_FILE = "Outputs.in"                             # Name of the output file
WORKING_DIR = os.getcwd()     # Directory in which the script is being executed
OS = sys.platform                               # Platform executing the script
USER = os.getenv("USER")                # Name of the user executing the script
TIMESTAMP = str(datetime.datetime.now()).split('.')[0]             # Time stamp
# --------------------------- FoamElasticModulus ------------------------------ #
def FoamElasticModulus(MU,SIGMA,strut_content,rho):
    # Main Inputs:
    AbaqusVersion='abq6101' #ABAQUS version
    NumOfCores='20' #Number of processors available in the machine(At least 1)
    ##############################################################################
    #Material properties:
    E='2.4e9'              #Elastic modulus of Solid Polyurethan (Pa)
    G='0.9e9'              #Shear modulus of Solid Polyurethan (Pa)
    NU='0.35'              #Poisson ratio of Solid Polyurethan 
    NumOfCells=40          #Number of cells in RVE (Representative Volume Element)
    MeshSize=0.1          #Size of mesh used to mesh the geometry
    AppliedStrain=-0.001   #Compressive strain 
    DensityOfSolidPU=1200  #Real density of Solid Polyurethan (kg/m3) 
    PUinWall=1-strut_content#Amount of material in walls (in percentage %)
    CubeLength=1           #Primary length of RVE (mm)
    ##############################################################################
    #Cross_sectional area function constants (Polynomial_power 4)
    A1,A2,A3,A4,A5=2.77,0.962,3.4033,-0.2291,1.0255
    ##############################################################################
    #Creating input file for spherepack:
    mypath=WORKING_DIR+'/'
    os.chdir(mypath)
    ExtraVolumeRatio=0.50
    SuccessInMeshing=1
    while SuccessInMeshing==1:
        SuccessInTessellation=1
        while SuccessInTessellation==1:
            myfile=os.path.join(mypath,'Project01.prj')
            f=open(myfile,'w')
            f.write('{0:34s}\n'.format('# RANDOM COORDINATES: 0 OR FILE: 1'))
            f.write('{0}\n{1}\n{2:d}\n{3}\n'.format('0','# NUMBER OF SPHERES',NumOfCells,'# LENGTH IN X-, Y- AND Z-DIRECTION'))
            f.write('{0:f}\n{0:f}\n{0:f}\n'.format(CubeLength))
            f.write('{0}\n{1}\n{2}\n{3}\n{4}\n'.format('# EPSILON','0.0001','# NTAU','963800000','# NOMINAL DENSITY'))
            f.write('{0}\n'.format('0.60'))
            f.write('{0}\n{1}\n{2}\n{3}\n{4}\n'.format('# MAX. AND MINIM. DIAMETER','3.0','1.0','# MAX. NUMBER OF STEPS','50000000'))
            f.write('{0}\n{1}\n{2}\n'.format('# DISTRIBUTION','3','# DISTRIBUTION PARAMETERS'))
            f.write('{0:f}\n{1:f}\n'.format(MU,SIGMA))
            f.write('{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n{8}\n'.format('0.4','-3.3','3','0.10','0.20','0.70','0.902113','3.000000','1.000000'))
            f.close()
            os.system('./spherepack')
            while True:
                time.sleep(1);files = os.listdir('.')
                if ('Project01.rco' in files):
                    break
            myfile12=os.path.join(mypath,'Project01.rco')
            CentersRads=np.loadtxt(myfile12,usecols = (0,1,2,3))
            CentersRads1=np.array([])
            ExtraVolume=ExtraVolumeRatio*CubeLength
            for i in range(int(len(CentersRads))):
                CentersRads1=np.append(CentersRads1,CentersRads[i][0]);CentersRads1=np.append(CentersRads1,CentersRads[i][1])
                CentersRads1=np.append(CentersRads1,CentersRads[i][2]);CentersRads1=np.append(CentersRads1,CentersRads[i][3])
            ExtraVolume=ExtraVolumeRatio*CubeLength
            for i in range(int(len(CentersRads))):
                if (CentersRads[i][0])<(ExtraVolume):
                    CentersRads1=np.append(CentersRads1,CubeLength+(ExtraVolume-CentersRads[i][0]))
                    CentersRads1=np.append(CentersRads1,CentersRads[i][1]);CentersRads1=np.append(CentersRads1,CentersRads[i][2])
                    CentersRads1=np.append(CentersRads1,CentersRads[i][3])
            CentersRads2=np.array([])
            for i in range(int(len(CentersRads1)/4)):
                CentersRads2=np.append(CentersRads2,CentersRads1[4*i]);CentersRads2=np.append(CentersRads2,CentersRads1[1+4*i])
                CentersRads2=np.append(CentersRads2,CentersRads1[2+4*i]);CentersRads2=np.append(CentersRads2,CentersRads1[3+4*i])
            for i in range(int(len(CentersRads1)/4)):
                if (CentersRads1[1+4*i])<(ExtraVolume):
                    CentersRads2=np.append(CentersRads2,CentersRads1[4*i])
                    CentersRads2=np.append(CentersRads2,CubeLength+(ExtraVolume-CentersRads1[1+4*i]))
                    CentersRads2=np.append(CentersRads2,CentersRads1[2+4*i]);CentersRads2=np.append(CentersRads2,CentersRads1[3+4*i])
            CentersRads3=np.array([])
            for i in range(int(len(CentersRads2)/4)):
                CentersRads3=np.append(CentersRads3,CentersRads2[4*i]);CentersRads3=np.append(CentersRads3,CentersRads2[1+4*i])
                CentersRads3=np.append(CentersRads3,CentersRads2[2+4*i]);CentersRads3=np.append(CentersRads3,CentersRads2[3+4*i])
            for i in range(int(len(CentersRads2)/4)):
                if (CentersRads2[2+4*i])<(ExtraVolume):
                    CentersRads3=np.append(CentersRads3,CentersRads2[4*i]);CentersRads3=np.append(CentersRads3,CentersRads2[1+4*i])
                    CentersRads3=np.append(CentersRads3,CubeLength+(ExtraVolume-CentersRads2[2+4*i]))
                    CentersRads3=np.append(CentersRads3,CentersRads2[3+4*i])
            myfile3=os.path.join(mypath,'Centers.txt')
            ff=open(myfile3,'w')
            for i in range(0,int(len(CentersRads3))/4):
                ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(CentersRads3[4*i],CentersRads3[1+4*i],CentersRads3[2+4*i]))
            ff.close()
            myfile4=os.path.join(mypath,'Rads.txt')
            fff=open(myfile4,'w')
            for i in range(0,int(len(CentersRads3))/4):
                fff.write('{0:f}\n'.format(CentersRads3[3+4*i]/2))
            fff.close()
            time.sleep(3)
            NumOfCells1=int(len(CentersRads3)/4)
            RVElength=CubeLength+(ExtraVolumeRatio*CubeLength)
            #commandTessellation="neper -T -n {0:d} -domain 'cube({1},{1},{1})' -morpho @{2}Centers.txt -weight @{2}Rads.txt -regularization 1 -mloop 1 -o {2}RVE27 -format geo -statcell vol -statedge length -statface area -statver x".format((NumOfCells1),RVElength,mypath)
            commandTessellation="neper -T -n {0:d} -domain 'cube({1},{1},{1})' -morphooptiini 'coo:file({2}Centers.txt),weight:file({2}Rads.txt)' -regularization 1 -mloop 1 -o {2}RVE27 -format geo -statcell vol -statedge length -statface area -statver x".format((NumOfCells1),RVElength,mypath)
            os.system(commandTessellation)
            time.sleep(5)
            files = os.listdir('.')
            if ('RVE27.geo' in files):
                SuccessInTessellation=0         
            os.remove('Project01.prj');os.remove('Project01.rst');os.remove('Project01.sco');os.remove('Project01.rco')         
            os.remove('Rads.txt');os.remove('Centers.txt')
        ##########################################################################
        myfile12=os.path.join(mypath,'RVE27.stver')
        verX=np.loadtxt(myfile12)
        NumOfNodes=len(verX)
        myfile12=os.path.join(mypath,'RVE27.stedge')
        EdgesLength=np.loadtxt(myfile12)
        NumOfEdges=len(EdgesLength)
        EdgeLengths=np.array([])
        for i in range(0,int(NumOfEdges)):
            c=np.array(EdgesLength[i])
            EdgeLengths=np.append(EdgeLengths,c)
        myfile12=os.path.join(mypath,'RVE27.stface')
        SurfacesArea=np.loadtxt(myfile12)
        NumOfSurfaces=len(SurfacesArea)
        myfile12=os.path.join(mypath,'RVE27.stcell')
        CellVolumes=np.loadtxt(myfile12)
        NumOfVolumes=len(CellVolumes)
        myfile12=os.path.join(mypath,'RVE27.geo')
        text_file = open(myfile12, "r")
        lines = text_file.readlines()
        text_file.close()
        os.remove('RVE27.geo')
        os.remove('RVE27.stcell');os.remove('RVE27.stedge');os.remove('RVE27.stface');os.remove('RVE27.stver')
        t=0
        Nodes=range(NumOfNodes)
        for i in range(0,NumOfNodes):
            currentline = lines[i].split("{")
            a=currentline[1].split("}")
            b=a[0].split(",")
            c=np.array(b)
            Nodes[t]=np.absolute(c.astype(np.float))
            t=t+1
        Edges=range(NumOfEdges)
        r=0
        t=0
        for i in range(NumOfNodes,NumOfNodes+NumOfEdges):
            currentline = lines[i].split("{")
            a=currentline[1].split("}")
            b=a[0].split(",")
            c=np.array(b)
            Edges[t]=np.absolute(c.astype(np.float))
            t=t+1
        Faces=range(NumOfSurfaces)
        r=0
        t=0
        for i in range(0,NumOfSurfaces):
            j=NumOfNodes+NumOfEdges+(2*i)
            currentline = lines[j].split("{")
            a=currentline[1].split("}")
            b=a[0].split(",")
            c=np.array(b)
            Faces[t]=(c.astype(np.float))
            t=t+1
        Volumes=range(NumOfVolumes)
        a=list(lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)])
        a[0:20]=''
        o=len(a)
        a[o-3:]=''
        lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)]="".join(a)
        currentline=lines[NumOfNodes+NumOfEdges+(2*NumOfSurfaces)].split(",")
        c=np.array(currentline)
        Volumes[0]=np.absolute(c.astype(np.float))
        t=1
        for i in range(1,2*NumOfVolumes):
            if i%2==0:
                j=NumOfNodes+NumOfEdges+(2*NumOfSurfaces)+i
                currentline = lines[j].split(" = {")
                a=currentline[1].split("}")
                b=a[0].split(",")
                c=np.array(b)
                Volumes[t]=np.absolute(c.astype(np.float))
                t=t+1
        ##########################################################################
        RemovedEdges=np.array([])
        for i in range(int(len(Edges))):
            X1=abs(Nodes[int(Edges[i][0])-1][0]);Y1=abs(Nodes[int(Edges[i][0])-1][1]);Z1=abs(Nodes[int(Edges[i][0])-1][2])
            X2=abs(Nodes[int(Edges[i][1])-1][0]);Y2=abs(Nodes[int(Edges[i][1])-1][1]);Z2=abs(Nodes[int(Edges[i][1])-1][2])
            if ((X1==0 or X1==RVElength) and X2==X1) or ((Y1==0 or Y1==RVElength) and Y2==Y1) or ((Z1==0 or Z1==RVElength) and Z2==Z1):
                RemovedEdges=np.append(RemovedEdges,i+1)
        RemovedFaces=np.array([])
        RemainedEdges=np.array([])
        for i in range(int(len(Faces))):
            FaceInsideCube=0
            for j in range(int(len(Faces[i]))):
                if abs(Faces[i][j]) in RemovedEdges:
                    FaceInsideCube=1+FaceInsideCube
            if FaceInsideCube==len(Faces[i]):
                RemovedFaces=np.append(RemovedFaces,int(i+1))
            if FaceInsideCube!=len(Faces[i]):
                for j in range(int(len(Faces[i]))):
                    RemainedEdges=np.append(RemainedEdges,abs(Faces[i][j]))
        ##########################################################################
        myfile13=os.path.join(mypath,'PeriodicRVE.geo')
        text_GEO = open(myfile13, "w")
        text_GEO.write('{0}{1}{2}\n'.format('cl=',MeshSize,';'))
        for i in range(int(len(Nodes))):
            lines[i] = lines[i].replace('};',',cl};')
            text_GEO.write('{0}'.format(lines[i]))
        for i in range(int(len(Edges))):
            if i+1 in RemainedEdges:
                text_GEO.write('{0}{1}{2}{3}{4}{5}{6}\n'.format('Line (',int(i+1),') = {',int(Edges[i][0]),',',int(Edges[i][1]),'};'))
        for i in range(int(len(Faces))):
            if i+1 not in RemovedFaces:
                text_GEO.write('{0}{1}{2}'.format('Line Loop (',int(i+1),') = {'))
                for j in range(int(len(Faces[i]))):
                    text_GEO.write('{0}'.format(int(Faces[i][j])))
                    if j+1!=len(Faces[i]):
                        text_GEO.write('{0}'.format(','))
                text_GEO.write('{0}\n'.format('};'))
                text_GEO.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}\n'.format('Plane Surface (',int(i+1),') = {',int(i+1),'}; Physical Surface (',int(i+1),') = {',int(i+1),'};'))
        text_GEO.close()
        TotalArea=0
        for i in range(0,int(len(SurfacesArea))):
            if i+1 not in RemovedFaces:
                c=np.array(SurfacesArea[i])
                TotalArea=TotalArea+(np.absolute(c.astype(np.float)))
        time.sleep(5)
        commandMeshing="gmsh PeriodicRVE.geo -2 -optimize_lloyd -o PeriodicRVE -format inp -saveall"
        os.system(commandMeshing)
        while True:
            time.sleep(4);files = os.listdir('.')
            if ('PeriodicRVE.inp' in files):
                SuccessInMeshing=0;break
        os.remove('PeriodicRVE.geo')
    # Creating Corrected INP To Be processed In Next Step
    NoOfBeamPerEdge=np.array([])
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
    i=i+1
    text_inp.write('{0}\n'.format('*ELEMENT, TYPE=B31, ELSET=auto1'))
    RemainedEdges1=np.array([])
    K=0
    N=0
    while K==0:
        if len(lines[i])>30:
            if lines[i][15]!='C':
                s=np.array(lines[i][31:int(len(lines[i])-1)])
                q=np.absolute(s.astype(np.float))
                M=0
                if q in RemovedEdges:
                    M=1
            if lines[i][15]=='C':
                K=1
        if 30>len(lines[i]) and M!=1:
            RemainedEdges1=np.append(RemainedEdges1,int(q))
            text_inp.write('{0}'.format(lines[i]))
            N=N+1
        i=i+1
    NoOfBeamPerEdge=itemfreq(RemainedEdges1)
    RemainedEdges1=np.unique(RemainedEdges1)
    text_inp.write('{0}\n'.format('*ELEMENT, TYPE=S3R, ELSET=auto1'))
    K=0
    while K==0:
        if ((lines[i][0]!='*') and (lines[i][3]!='S')):
            text_inp.write('{0}'.format(lines[i]))
        if lines[i][3]=='S':
            K=1
        i=i+1
    text_inp.write('{0}'.format('*******'))
    text_inp.close();time.sleep(10)
    myfile12=os.path.join(mypath,'CorrectedINPafterGMSH.inp')
    INPgeom = open(myfile12, "r")
    lines = INPgeom.readlines()
    INPgeom.close()
    os.remove('CorrectedINPafterGMSH.inp')
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
    NodeNums=np.array([])
    k=1
    for i in range(int(len(Shells)/4)):
        NodeNums=np.append(NodeNums,int(Shells[k]))
        NodeNums=np.append(NodeNums,int(Shells[k+1]))
        NodeNums=np.append(NodeNums,int(Shells[k+2]))
        k=k+4
    BeamLengths=np.array([])
    k=1
    for i in range(int(len(Beams)/3)):
       x1=Nodes[1+int(4*(Beams[k]-1))];y1=Nodes[2+int(4*(Beams[k]-1))];z1=Nodes[3+int(4*(Beams[k]-1))]
       x2=Nodes[1+int(4*(Beams[k+1]-1))];y2=Nodes[2+int(4*(Beams[k+1]-1))];z2=Nodes[3+int(4*(Beams[k+1]-1))]
       BeamLengths=np.append(BeamLengths,((((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))**0.5))
       k=k+3
    RemainedNodesNum=np.unique(NodeNums)
    RemainedNodesBoundary=np.array([])       
    for k in range(int(len(RemainedNodesNum))):
        if abs(Nodes[1+(4*(RemainedNodesNum[k]-1))])==RVElength or abs(Nodes[1+(4*(RemainedNodesNum[k]-1))])==0 or abs(Nodes[2+(4*(RemainedNodesNum[k]-1))])==RVElength or abs(Nodes[2+(4*(RemainedNodesNum[k]-1))])==0 or abs(Nodes[3+(4*(RemainedNodesNum[k]-1))])==RVElength or abs(Nodes[3+(4*(RemainedNodesNum[k]-1))])==0:
            RemainedNodesBoundary=np.append(RemainedNodesBoundary,int(RemainedNodesNum[k]))
            RemainedNodesBoundary=np.append(RemainedNodesBoundary,Nodes[1+(4*(RemainedNodesNum[k]-1))])
            RemainedNodesBoundary=np.append(RemainedNodesBoundary,Nodes[2+(4*(RemainedNodesNum[k]-1))])
            RemainedNodesBoundary=np.append(RemainedNodesBoundary,Nodes[3+(4*(RemainedNodesNum[k]-1))])
    A=10000
    NodePairsXY=np.array([])
    m=0
    n=0
    for i in range(int(len(RemainedNodesBoundary)/4)):
        for j in range(i,int(len(RemainedNodesBoundary)/4)):
            if RemainedNodesBoundary[m]!=RemainedNodesBoundary[n] and int(A*RemainedNodesBoundary[m+1])==int(A*RemainedNodesBoundary[n+1]) and int(A*RemainedNodesBoundary[m+2])==int(A*RemainedNodesBoundary[n+2]) and (abs(RemainedNodesBoundary[m+3]-RemainedNodesBoundary[n+3])>0.98*RVElength):
                NodePairsXY=np.append(NodePairsXY,RemainedNodesBoundary[m])
                NodePairsXY=np.append(NodePairsXY,RemainedNodesBoundary[n])
            n=n+4
        m=m+4
        n=m
    NodePairsXZ=np.array([])
    m=0
    n=0
    for i in range(int(len(RemainedNodesBoundary)/4)):
        for j in range(i,int(len(RemainedNodesBoundary)/4)):
            if RemainedNodesBoundary[m]!=RemainedNodesBoundary[n] and int(A*RemainedNodesBoundary[m+1])==int(A*RemainedNodesBoundary[n+1]) and int(A*RemainedNodesBoundary[m+3])==int(A*RemainedNodesBoundary[n+3]) and (abs(RemainedNodesBoundary[m+2]-RemainedNodesBoundary[n+2])>0.98*RVElength):# and (RemainedNodes[m+1]<1.3 or RemainedNodes[m+1]>1.7) and (RemainedNodes[m+3]<1.3 or RemainedNodes[m+3]>1.7):
                NodePairsXZ=np.append(NodePairsXZ,RemainedNodesBoundary[m])
                NodePairsXZ=np.append(NodePairsXZ,RemainedNodesBoundary[n])
            n=n+4
        m=m+4
        n=m
    NodePairsYZ=np.array([])
    m=0
    n=0
    for i in range(int(len(RemainedNodesBoundary)/4)):
        for j in range(i,int(len(RemainedNodesBoundary)/4)):
            if RemainedNodesBoundary[m]!=RemainedNodesBoundary[n] and int(A*RemainedNodesBoundary[m+2])==int(A*RemainedNodesBoundary[n+2]) and int(A*RemainedNodesBoundary[m+3])==int(A*RemainedNodesBoundary[n+3]) and (abs(RemainedNodesBoundary[m+1]-RemainedNodesBoundary[n+1])>0.98*RVElength):# and (RemainedNodes[m+3]<1.0 or RemainedNodes[m+3]>1.2) and (RemainedNodes[m+2]<1.3 or RemainedNodes[m+2]>1.7):
                NodePairsYZ=np.append(NodePairsYZ,RemainedNodesBoundary[m])
                NodePairsYZ=np.append(NodePairsYZ,RemainedNodesBoundary[n])
            n=n+4
        m=m+4
        n=m
    NodeCordinates=np.array([])
    for i in range(int(len(RemainedNodesNum))):
        K=4*(RemainedNodesNum[i]-1)
        NodeCordinates=np.append(NodeCordinates,Nodes[K+1])
        NodeCordinates=np.append(NodeCordinates,Nodes[K+2])
        NodeCordinates=np.append(NodeCordinates,Nodes[K+3])
    UpLevelOfRVE=1+int(max(NodeCordinates))
    DistCenter=np.array([]);DistX=np.array([]);DistY=np.array([]);DistZ=np.array([])
    NodeNums1=np.array([])
    for i in range(int(len(RemainedNodesNum))):
        K=4*(RemainedNodesNum[i]-1)
        NodeNums1=np.append(NodeNums1,Nodes[K])
        DistCenter=np.append(DistCenter,(Nodes[K+1]**2)+(Nodes[K+2]**2)+(Nodes[K+3]**2))
        DistX=np.append(DistX,((Nodes[K+1]-UpLevelOfRVE)**2)+((Nodes[K+2])**2)+((Nodes[K+3])**2))
        DistY=np.append(DistY,((Nodes[K+1])**2)+((Nodes[K+2]-UpLevelOfRVE)**2)+((Nodes[K+3])**2))
        DistZ=np.append(DistZ,((Nodes[K+1])**2)+((Nodes[K+2])**2)+((Nodes[K+3]-UpLevelOfRVE)**2))
    CornerNode=NodeNums1[DistCenter.argmin()] 
    CornerNodeX=NodeNums1[DistX.argmin()]
    CornerNodeY=NodeNums1[DistY.argmin()]
    CornerNodeZ=NodeNums1[DistZ.argmin()]
    LenNodePairsXY=len(NodePairsXY);LenNodePairsXZ=len(NodePairsXZ);LenNodePairsYZ=len(NodePairsYZ)
    K=0
    L=0
    for i in range(int(LenNodePairsXY)):
        for j in range(int(LenNodePairsXZ/2)):
            if NodePairsXY[K]==NodePairsXZ[L] or NodePairsXY[K]==NodePairsXZ[L+1]:
                NodePairsXZ[L]=123456;NodePairsXZ[L+1]=123456
            L=L+2
        L=0
        K=K+1
    K=0
    L=0
    for i in range(int(LenNodePairsXY)):
        for j in range(int(LenNodePairsYZ/2)):
            if NodePairsXY[K]==NodePairsYZ[L] or NodePairsXY[K]==NodePairsYZ[L+1]:
                NodePairsYZ[L]=123456;NodePairsYZ[L+1]=123456
            L=L+2
        L=0
        K=K+1
    K=0
    L=0
    for i in range(int(LenNodePairsXZ)):
        for j in range(int(LenNodePairsYZ/2)):
            if NodePairsXZ[K]==NodePairsYZ[L] or NodePairsXZ[K]==NodePairsYZ[L+1]:
                NodePairsYZ[L]=123456;NodePairsYZ[L+1]=123456
            L=L+2
        L=0
        K=K+1
    L=0
    for i in range(int(LenNodePairsXY/2)):
        if CornerNode==NodePairsXY[L] or CornerNodeX==NodePairsXY[L] or CornerNodeY==NodePairsXY[L] or CornerNodeZ==NodePairsXY[L] or CornerNode==NodePairsXY[L+1] or CornerNodeX==NodePairsXY[L+1] or CornerNodeY==NodePairsXY[L+1] or CornerNodeZ==NodePairsXY[L+1]:
            NodePairsXY[L]=123456;NodePairsXY[L+1]=123456
        L=L+2
    L=0
    for i in range(int(LenNodePairsXZ/2)):
        if CornerNode==NodePairsXZ[L] or CornerNodeX==NodePairsXZ[L] or CornerNodeY==NodePairsXZ[L] or CornerNodeZ==NodePairsXZ[L] or CornerNode==NodePairsXZ[L+1] or CornerNodeX==NodePairsXZ[L+1] or CornerNodeY==NodePairsXZ[L+1] or CornerNodeZ==NodePairsXZ[L+1]:
            NodePairsXZ[L]=123456;NodePairsXZ[L+1]=123456
        L=L+2
    L=0
    for i in range(int(LenNodePairsYZ/2)):
        if CornerNode==NodePairsYZ[L] or CornerNodeX==NodePairsYZ[L] or CornerNodeY==NodePairsYZ[L] or CornerNodeZ==NodePairsYZ[L] or CornerNode==NodePairsYZ[L+1] or CornerNodeX==NodePairsYZ[L+1] or CornerNodeY==NodePairsYZ[L+1] or CornerNodeZ==NodePairsYZ[L+1]:
            NodePairsYZ[L]=123456;NodePairsYZ[L+1]=123456
        L=L+2
    Ylength=RVElength
    RVEarea=RVElength*RVElength
    Vfoam=RVEarea*Ylength 
    #How material distributed between walls and struts (Vw , Vs) :
    SolidPUmass=rho*Vfoam           #rho is the real density of the foam (kg/m3)
    VsolidPUinRVE=SolidPUmass/DensityOfSolidPU
    PUinStruts=1-PUinWall
    Vw=VsolidPUinRVE*PUinWall       #Volume of solid PU in walls
    Vs=VsolidPUinRVE*PUinStruts     #Volume of solid PU in Struts
    TotalLength=0
    for i in range(int(len(RemainedEdges1))):
        TotalLength=TotalLength+EdgeLengths[RemainedEdges1[i]-1]
    VolumeRatioStrut=np.array([])
    for i in range(int(len(RemainedEdges1))):
        VolumeRatioStrut=np.append(VolumeRatioStrut,EdgeLengths[RemainedEdges1[i]-1]/TotalLength)
    myfile13=os.path.join(mypath,'INP.inp')
    text_inp = open(myfile13, "w")
    text_inp.write('{0}{1}\n'.format('****Mesh Size=',MeshSize))
    text_inp.write('{0}{1}\n'.format('*include, input=','PBC.txt'))
    for i in range(0,len(lines)):
        text_inp.write('{0}'.format(lines[i]))
        if lines[i][1]=='E':
            break
    k=0
    for i in range(int(len(Beams)/3)):
        text_inp.write('{0},{1},{2}\n'.format(int(Beams[k]),int(Beams[k+1]),int(Beams[k+2])))
        k=k+3
    text_inp.write('*ELEMENT, TYPE=S3R, ELSET=auto1\n')
    k=0
    for i in range(int(len(Shells)/4)):
        text_inp.write('{0},{1},{2},{3}\n'.format(int(Shells[k]),int(Shells[k+1]),int(Shells[k+2]),int(Shells[k+3])))
        k=k+4
    #Section assignment for Beam Elements
    TotalVolAssignedBeams=0
    BeamElCounter=1
    for i in range(int(len(NoOfBeamPerEdge))):
        NumElinEdge=NoOfBeamPerEdge[i][1] 
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
            Vbeam1=(Vbeam[0]/TotalAreaStrut)*(Vs*VolumeRatioStrut[i])
            Area=Vbeam1/BeamLengths[BeamElCounter-1];c=0.000001;h=(Area*(3.0**0.5))**0.5;a=h*(2.0/(3.00**0.5));d=h/3.00
            K=K+DeltaX
            text_inp.write('{0:.10f}{1}{2:.10f}{3}{4:.10f}{5}{6:.10f}\n'.format(a,', ',h,', ',c,', ',d))
            text_inp.write('{0}\n'.format('0.,0.,-1'))
            TotalVolAssignedBeams=TotalVolAssignedBeams+(Area*BeamLengths[BeamElCounter-1])
            BeamElCounter=BeamElCounter+1
    #Calculation of shell thickness:
    ShellThickness=Vw/TotalArea
    #Section assignment for shell elements:
    t=(len(Shells)/4)%16
    tt=(int((len(Shells)/4)-t))/16
    N=0
    M=1
    text_inp.write('{0}{1}\n'.format('*Elset, elset=EShell',int(0)))
    for i in range(int(tt)):
        for j in range(16):
            text_inp.write('{0}{1}'.format(int(Shells[N]),','))
            N=N+4
        M=M+1
        text_inp.write('\n')
    if t>0:
        for j in range(int(t)):
            text_inp.write('{0}{1}'.format(int(Shells[N]),','))
            N=N+4
    if t>0:
        text_inp.write('\n')
    text_inp.write('{0}{1}{2}\n'.format('*Shell Section, elset=EShell',int(0),', material=Material-1'))
    text_inp.write('{0}{1}\n'.format(ShellThickness,', 5'));text_inp.write('{0}\n'.format('*Material, name=Material-1'))
    text_inp.write('{0}\n'.format('*Density'));text_inp.write('{0}\n'.format(DensityOfSolidPU))
    text_inp.write('{0}\n'.format('*Elastic'));text_inp.write('{0}{1}{2}\n'.format(E,', ',NU))
    text_inp.write('{0}\n'.format('**'));text_inp.write('{0}{1}\n'.format('**RVEarea',RVEarea))
    text_inp.write('{0}{1}\n'.format('**Ylength',Ylength))
    text_inp.write('{0}\n'.format('** STEP: Step-1'));text_inp.write('{0}\n'.format('*Step, name=Step-1, nlgeom=YES, INC=1000'))
    text_inp.write('{0}\n'.format('*Static, stabilize=0.001'));text_inp.write('{0}\n'.format('0.001, 1., 1e-09, 0.01'))
    text_inp.write('{0}\n'.format('** BOUNDARY CONDITIONS'));text_inp.write('{0}\n'.format('*Nset, nset=Znodepairs'))
    L=0
    for j in range(int(LenNodePairsXY/2)):
        if NodePairsXY[L]!=123456 and NodePairsXY[L+1]!=123456:
            text_inp.write('{0}{1}{2}{3}\n'.format(int(NodePairsXY[L]),',',int(NodePairsXY[L+1]),','))
        L=L+2
    text_inp.write('{0}\n'.format('*Nset, nset=Ynodepairs'))
    L=0
    for j in range(int(LenNodePairsXZ/2)):
        if NodePairsXZ[L]!=123456 and NodePairsXZ[L+1]!=123456:
            text_inp.write('{0}{1}{2}{3}\n'.format(int(NodePairsXZ[L]),',',int(NodePairsXZ[L+1]),','))
        L=L+2
    text_inp.write('{0}\n'.format('*Nset, nset=Xnodepairs'))
    L=0
    for j in range(int(LenNodePairsYZ/2)):
        if NodePairsYZ[L]!=123456 and NodePairsYZ[L+1]!=123456:
            text_inp.write('{0}{1}{2}{3}\n'.format(int(NodePairsYZ[L]),',',int(NodePairsYZ[L+1]),','))
        L=L+2
    text_inp.write('{0}\n{1}\n'.format('*Nset, nset=MN0',int(CornerNode)));text_inp.write('{0}\n{1}\n'.format('*Nset, nset=MNx',int(CornerNodeX)))
    text_inp.write('{0}\n{1}\n'.format('*Nset, nset=MNy',int(CornerNodeY)));text_inp.write('{0}\n{1}\n'.format('*Nset, nset=MNz',int(CornerNodeZ)))
    text_inp.write('{0}\n{1}\n{2}\n{3}\n{4}\n'.format('*Boundary','MN0, PINNED','*Boundary','MNx, 2, 2','MNx, 3, 3'))
    text_inp.write('{0}\n{1}\n'.format('*Boundary','MNy, 1, 1'))
    text_inp.write('{0}{1}\n'.format('MNy, 2, 2, ',Ylength*(AppliedStrain)))
    text_inp.write('{0}\n{1}\n{2}\n{3}\n{4}\n'.format('MNy, 3, 3','*Boundary','MNz, 1, 1','MNz, 2, 2','** OUTPUT REQUESTS'))
    text_inp.write('{0}\n{1}\n'.format('*Restart, write, frequency=0','** FIELD OUTPUT: F-Output-1'))
    text_inp.write('{0}\n{1}\n'.format('*Output, field, variable=PRESELECT','** HISTORY OUTPUT: H-Output-1'))
    text_inp.write('{0}\n{1}\n'.format('*Output, history, variable=PRESELECT','*End Step'))
    text_inp.close()
    myfile14=os.path.join(mypath,'PBC.txt')
    t=open(myfile14, "w")
    t.write('{0}\n'.format('*EQUATION'))
    L=0
    for j in range(int(LenNodePairsXY/2)):
        if NodePairsXY[L]!=123456 and NodePairsXY[L+1]!=123456:
            for k in range(6):
                t.write('{0}\n'.format('3'))
                t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsXY[L+1]),',', k+1,' ,1.00,	     ',int(NodePairsXY[L]),', ',k+1,' ,-1.00,	     ',int(CornerNodeZ),',', k+1,' ,-1.00'))
        L=L+2
    L=0
    for j in range(int(LenNodePairsXZ/2)):
        if NodePairsXZ[L]!=123456 and NodePairsXZ[L+1]!=123456:
            for k in range(6):
                t.write('{0}\n'.format('3'))
                t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsXZ[L+1]),',', k+1,' ,1.00,	     ',int(NodePairsXZ[L]),', ',k+1,' ,-1.00,	     ',int(CornerNodeY),',', k+1,' ,-1.00'))
        L=L+2
    L=0
    for j in range(int(LenNodePairsYZ/2)):
        if NodePairsYZ[L]!=123456 and NodePairsYZ[L+1]!=123456:
            for k in range(6):
                t.write('{0}\n'.format('3'))
                t.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}\n'.format('	    ',int(NodePairsYZ[L+1]),',', k+1,' ,1.00,	     ',int(NodePairsYZ[L]),', ',k+1,' ,-1.00,	     ',int(CornerNodeX),',', k+1,' ,-1.00'))
        L=L+2
    t.close()
    ########################################################################################################################
    ########################################### Finite Element Analysis (ABAQUS) ###########################################
    ########################################################################################################################
    while True:
        time.sleep(1)
        files = os.listdir('.')
        if ('INP.inp' in files) and ('PBC.txt' in files):
            break
    CMDforABAQUS='cd {0} && {1} job=INP cpus={2}'.format(mypath,AbaqusVersion,NumOfCores)
    os.system(CMDforABAQUS)
    while True:
        time.sleep(1)
        files = os.listdir('.')
        if ('INP.lck' not in files) and ('INP.odb' in files):
            break
    #Creation of RF2 and U2 report file:
    myfile14=os.path.join(mypath,'AbaqusPythonInpRF2U2.py')
    q=open(myfile14, "w")
    q.write('{0}\n{1}\n{2}\n{3}\n'.format('from abaqus import *','from abaqusConstants import *','from viewerModules import *','from driverUtils import executeOnCaeStartup'))
    q.write('{0}\n{1}\n{2}\n'.format('executeOnCaeStartup()','o2 = session.openOdb(',"    name='INP.odb')"))
    q.write('{0}\n'.format("session.viewports['Viewport: 1'].setValues(displayedObject=o2)"))
    q.write('{0}\n'.format("session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=("))
    q.write('{0}\n'.format("    CONTOURS_ON_DEF, ))"))
    q.write('{0}{1}{2}\n'.format("odb = session.odbs['",mypath,"INP.odb']"))
    q.write('{0}\n'.format("session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF',"))
    q.write('{0}\n'.format("    NODAL, ((COMPONENT, 'RF2'), )), ('U', NODAL, ((COMPONENT, 'U2'), )), ),"))
    q.write('{0}{1}{2}\n'.format("    nodeLabels=(('PART-1-1', ('",int(CornerNodeY),"', )), ))"))
    q.write('{0}{1}{2}\n'.format("x0 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: ",int(CornerNodeY),"']" ))
    q.write('{0}{1}{2}\n'.format("x1 = session.xyDataObjects['U:U2 PI: PART-1-1 N: ",int(CornerNodeY),"']"))
    q.write('{0}\n'.format("session.xyReportOptions.setValues()"))
    q.write('{0}\n'.format("session.writeXYReport(fileName='XYdataRF2U2.rpt', xyData=(x0, x1))"))
    q.close()
    ########################################################################################################################
    ########################################### Output Data Base (ODB of ABAQUS) ###########################################
    ########################################################################################################################
    os.system('cd {0} && {1} cae noGUI=AbaqusPythonInpRF2U2.py'.format(mypath,AbaqusVersion)) #Executing .py to get the XY data from ODB
    while True:
        time.sleep(1)
        files = os.listdir('.')
        if ('XYdataRF2U2.rpt' in files):
            break
    myfile14=os.path.join(mypath,'XYdataRF2U2.rpt')
    text_file = open(myfile14, "r")
    lines = text_file.readlines()
    os.remove('INP.com');os.remove('INP.dat');os.remove('INP.inp');os.remove('INP.log');os.remove('INP.msg');os.remove('INP.prt')
    os.remove('INP.sim');os.remove('INP.sta');os.remove('INP.odb');os.remove('XYdataRF2U2.rpt');os.remove('PBC.txt');os.remove('AbaqusPythonInpRF2U2.py')
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
    return ElasticModulusOfRVE

# ------------------ Read input file && Printing E based on INPUTS------------------ #
assert os.path.isfile(INPUT_FILE), "Did not find file '%s'" %(INPUT_FILE)
print "Post processing..."
Outputfile=os.path.join(WORKING_DIR,OUTPUT_FILE)
Output=open(Outputfile,'w')
Input= np.loadtxt(INPUT_FILE,dtype=float,delimiter=",",comments='#',skiprows=0)
myfile14=os.path.join(WORKING_DIR,INPUT_FILE)
text_file_Inp = open(myfile14, "r")
if len(text_file_Inp.readlines())<2:
    text_file_Inp.close()
    MU,SIGMA,strut_content,rho=Input[0],Input[1],Input[2],Input[3]
    EModulus=FoamElasticModulus(MU,SIGMA,strut_content,rho)
    Output.write('{0:f}\n'.format(EModulus))
else:
    for No in range(int(len(Input))):
        MU,SIGMA,strut_content,rho=Input[No][0],Input[No][1],Input[No][2],Input[No][3]
        EModulus=FoamElasticModulus(MU,SIGMA,strut_content,rho)
        Output.write('{0:f}\n'.format(EModulus))
Output.close()