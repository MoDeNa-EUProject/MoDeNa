"""
@file
Packing and tesselation. Prepares representative volume element (RVE) of foam.

@author    Mohammad Marvi-Mashhadi (IMDEA Materials)
@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   foam_constr
"""
import os
import os.path
import numpy as np
import re
import time
import random
import math
mypath=os.getcwd()
## Main function.
#
#  Packing algorithm for creation of seeds for tessellation.
#  Tessellation by NEPER.
#  Saves `.geo` and `.gnu` files with tessellation for later voxelization.
#  Prepares input file for Hypermesh.
def main(MU,SIGMA,NumOfCells,filenameOut,packing,alternativePackingAlgorithm,
    tesselation,visualizeTesselation,geometry,statistics,hypermesh,deleteFiles,
    dx,dy,dz):
    #####To create input file for SpherePack
    if packing:
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
                PickCenterX,PickCenterY,PickCenterZ=\
                    lc*random.random(),lc*random.random(),lc*random.random()
                while (lc-Rads1[j]>=PickCenterX and lc-Rads1[j]>=PickCenterY
                    and lc-Rads1[j]>=PickCenterZ and Rads1[j]<PickCenterX
                    and Rads1[j]<PickCenterY and Rads1[j]<PickCenterZ):
                    PickCenterX,PickCenterY,PickCenterZ=\
                        lc*random.random(),lc*random.random(),lc*random.random()
                centers[j][0],centers[j][1],centers[j][2]=\
                    PickCenterX,PickCenterY,PickCenterZ
                KeepCentreX, KeepCentreY, KeepCentreZ, KeepR=\
                    PickCenterX, PickCenterY, PickCenterZ, Rads1[j]
                if j>0:
                    for t in range(0,j):
                        if (((( ((KeepCentreX-centers[t][0])**2.00)+
                            ((KeepCentreY-centers[t][1])**2.00)+
                            ((KeepCentreZ-centers[t][2])**2.00))**0.50)-
                            (KeepR+Rads1[t]))<0.000) and t!=j:
                            centers[j][0],centers[j][0],centers[j][0]=0,0,0
                            j=j-1
                            break
        mypath=os.getcwd()
        myfile=os.path.join(mypath,'Project01.rco')
        f=open(myfile,'w')
        for i in range(NumSpheres):
            f.write('{0:f}{1}{2:f}{3}{4:f}{5}{6:f}\n'.format(
                centers[i][0],'	',centers[i][1],'	',centers[i][2],'	',
                2.0*Rads1[i]))
        f.close()
        MAXcenters=max(centers)
        Mincenters=min(centers)
        EdgeCubeSize=[math.ceil(MAXcenters[0]-Mincenters[0]),
            math.ceil(MAXcenters[1]-Mincenters[1]),
            math.ceil(MAXcenters[2]-Mincenters[2])]
        EdgeRVESize=int(3.0*max(EdgeCubeSize)) #For NEPER: Size of edge of RVE
        if alternativePackingAlgorithm:
            os.system('./spherepack')
    if tesselation:
        myfile12=os.path.join(mypath,'Project01.rco')
        CentersRads=np.loadtxt(myfile12,usecols = (0,1,2,3))
        Centers=CentersRads[:,:3]  #All centers of spheres
        Rads=CentersRads[:,3]   #All radii of spheres
        Rads=Rads/2
        MaxCenters=np.amax(Centers, axis=0)
        L2=int(0.5+MaxCenters[0])
        print (L2)
        X=L2+Centers[:,0]
        Y=L2+Centers[:,1]
        Z=L2+Centers[:,2]
        Centers27=np.array([
            X,X+L2,X-L2,X+L2,X,X+L2,X,X-L2,X-L2,X-L2,X+L2,X,X,X-L2,X+L2,
            X+L2,X,X,X-L2,X,X,X-L2,X+L2,X+L2,X-L2,X-L2,X+L2,Y,Y+L2,Y-L2,
            Y+L2,Y+L2,Y,Y-L2,Y,Y-L2,Y+L2,Y-L2,Y-L2,Y+L2,Y,Y,Y,Y+L2,Y,Y,
            Y-L2,Y,Y+L2,Y-L2,Y+L2,Y-L2,Y+L2,Y-L2,Z,Z+L2,Z-L2,Z,Z+L2,Z+L2,
            Z-L2,Z-L2,Z,Z,Z,Z+L2,Z-L2,Z+L2,Z-L2,Z,Z,Z+L2,Z,Z,Z-L2,Z+L2,
            Z+L2,Z-L2,Z+L2,Z-L2,Z-L2])
        # Creation of input .txt files for neper
        myfile3=os.path.join(mypath,'Centers.txt')
        ff=open(myfile3,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(Centers27[i,j],
                    Centers27[i+27,j],Centers27[i+54,j]))
        ff.close()
        myfile4=os.path.join(mypath,'Rads.txt')
        fff=open(myfile4,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                fff.write('{0:f}\n'.format(Rads[j]))
        fff.close()
        commandTessellation="neper -T -n {0:d} -domain \
            'cube({1:d},{2:d},{3:d})' -morpho @Centers.txt -weight @Rads.txt \
            -regularization 1 -mloop 15 \
            -o RVE27 -format tess,geo \
            -statcell vol -statedge length -statface area \
            -statver x".format((27*NumSpheres),EdgeRVESize,EdgeRVESize,
            EdgeRVESize)
        os.system(commandTessellation)
        if visualizeTesselation: # needs POV-Ray
            commandVisualization="neper -V RVE27.tess -datacellcol ori \
                -datacelltrs 0.5 -showseed all -dataseedrad @Rads.txt \
                -dataseedtrs 1.0 -print RVE27"
            os.system(commandVisualization)
    ################################################################
    ######Extraction of middle Representative volume element########
    ################################################################
    if geometry:
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
        #print NumOfNodes+NumOfEdges+(2*NumOfSurfaces)
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
        myfile13=os.path.join(mypath,filenameOut+".geo")
        text_GEO = open(myfile13, "w")
        for i in range(int(MaxIndexOfNodes)):
            text_GEO.write('{0}'.format(lines[i]))
        for i in range(int(NumOfNodes),int(NumOfNodes+MaxIndexOfEdges)):
            text_GEO.write('{0}'.format(lines[i]))
        j=int(NumOfNodes+NumOfEdges)
        for i in range(int(NumOfNodes+NumOfEdges),
            int(NumOfNodes+NumOfEdges+MaxIndexOfFaces)):
            text_GEO.write('{0}{1}'.format(lines[j],lines[j+1]))
            j=2+j
        text_GEO.write('{0}{1}{2}\n'.format(' Surface Loop (1) = {',
            lines[int(NumOfNodes+NumOfEdges+(2*NumOfSurfaces))],'};'))
        NumOfCells=NumOfVolumes/27
        j=int(NumOfNodes+NumOfEdges+(2*NumOfSurfaces)+1)
        for i in range(int(NumOfCells-1)):
            text_GEO.write('{0}'.format(lines[j]))
            j+=1
        text_GEO.write('{0}{1}{2}{3}{4}'.format(
            'Volume (',NumOfCells,') = {',NumOfCells,'};'))
        text_GEO.close()
        ####################################################
        # Making gnuplot file containing Periodic RVE
        myfile13=os.path.join(mypath,filenameOut+".gnu")
        text_gnu = open(myfile13, "w")
        for i in range(int(NumOfNodes),int(NumOfNodes+MaxIndexOfEdges)):
            a=re.findall("[-+]?\d+[\.]?\d*", lines[i])
            b=int(a[1])-1
            c=int(a[2])-1
            text_gnu.write('{0}\t{1}\t{2}\n'.format(
                (Nodes[b][0]-L2)/L2,(Nodes[b][1]-L2)/L2,(Nodes[b][2]-L2)/L2))
            text_gnu.write('{0}\t{1}\t{2}\n\n\n'.format(
                (Nodes[c][0]-L2)/L2,(Nodes[c][1]-L2)/L2,(Nodes[c][2]-L2)/L2))
        text_gnu.close()
    # End of Geometry construction
    ################################################################
    ######Statistical info for Periodic RVE#########################
    ################################################################
    if statistics:
        myfile14=os.path.join(mypath,'EdgeLengths_PeriodicRVE.txt')
        EdgesLengthsOfMiddleRVE= open(myfile14, "w")
        for i in range(int(MaxIndexOfEdges)):
            EdgesLengthsOfMiddleRVE.write('{0}\n'.format(EdgesLength[i]))
        EdgesLengthsOfMiddleRVE.close()
        myfile15=os.path.join(mypath,'FaceAreas_PeriodicRVE.txt')
        AreaOfFaceOfMiddleRVE= open(myfile15, "w")
        for i in range(int(MaxIndexOfFaces)):
            AreaOfFaceOfMiddleRVE.write('{0}\n'.format(SurfacesArea[i]))
        AreaOfFaceOfMiddleRVE.close()
        myfile16=os.path.join(mypath,'CellVolumes_PeriodicRVE.txt')
        VolumeOfCellOfMiddleRVE= open(myfile16, "w")
        for i in range(int(NumOfCells)):
            VolumeOfCellOfMiddleRVE.write('{0}\n'.format(CellVolumes[i]))
        VolumeOfCellOfMiddleRVE.close()
    ################################################################
    ######Creation of input for HyperMesh###########################
    ################################################################
    if hypermesh:
        myfile17=os.path.join(mypath,'HyperMeshInput.cmf')
        HyperMeshInput= open(myfile17, "w")
        HyperMeshInput.write('{0}\n'.format('*cleanuptoleranceset(0.001)'))
        HyperMeshInput.write('{0}\n'.format('*toleranceset(0.001)'))
        for i in range(int(MaxIndexOfNodes)):
            HyperMeshInput.write('{0}{1}{2}{3}{4}{5}{6}\n'.format(
                '*createnode(',Nodes[i][0],',',Nodes[i][1],',',Nodes[i][2],
                ',0,0,0)'))
        for i in range(int(MaxIndexOfEdges)):
            HyperMeshInput.write('{0}\n{1}{2}{3}{4}\n'.format(
                '*linecreatefromnodes(1,0,150,5,179)','*createlist(nodes,1) ',
                int(Edges[i][0]),' ',int(Edges[i][1])))
        for i in range(int(MaxIndexOfFaces)):
            HyperMeshInput.write('{0}\n'.format('*surfacemode(4)'))
            HyperMeshInput.write('{0}'.format('*createmark(lines,1) '))
            for j in range(len(Faces[i])):
                HyperMeshInput.write('{0}\t'.format(int(Faces[i][j])))
            HyperMeshInput.write('\n{0}\n{1}\n'.format(
                '*createplane(1,1,0,0,0,0,0)','*splinesurface(lines,1,1,1,1)'))
        HyperMeshInput.write('{0}\n'.format('*createmark(lines,1) "all"'))
        HyperMeshInput.write('{0}\n'.format('*deletemark(lines,1)'))
        HyperMeshInput.write('{0}\n'.format('*settopologyedgedisplay(2,0)'))
        HyperMeshInput.write('{0}\n'.format('*plot()'))
        HyperMeshInput.write('{0}\n'.format('*settopologyedgedisplay(3,0)'))
        HyperMeshInput.write('{0}\n'.format('*plot()'))
        HyperMeshInput.write('{0}\n'.format('*settopologyedgedisplay(1,0)'))
        HyperMeshInput.write('{0}\n'.format('*plot()'))
        Tol=0.001
        for i in range(NumOfCells):
            HyperMeshInput.write('{0}\n'.format(
                '*createmark(surfaces,1) "displayed"'))
            HyperMeshInput.write('{0}{1}{2}{3}{4}\n'.format(
                '*selfstitchcombine(1,82,',Tol,',',Tol,')'))
            Tol=0.0001+Tol
        HyperMeshInput.close()
    ################################################################
    ######Removing unrequired files#################################
    ################################################################
    if deleteFiles:
        os.chdir(mypath)
        os.remove('Project01.prj')
        os.remove('Project01.rco')
        os.remove('Project01.rst')
        os.remove('Project01.sco')
        os.remove('Rads.txt')
        os.remove('Centers.txt')
        os.remove('RVE27.geo')
        os.remove('RVE27.stcell')
        os.remove('RVE27.stedge')
        os.remove('RVE27.stface')
        os.remove('RVE27.stver')
        os.remove('RVE27.tess')
