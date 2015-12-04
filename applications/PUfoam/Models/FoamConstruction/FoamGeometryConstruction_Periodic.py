__author__ = 'Mohammad Marvi-Mashhadi (IMDEA Materials)'
import os
import os.path
import numpy as np
mypath=os.getcwd()
def main(MU,SIGMA,NumOfCells,filenameOut,packing,tesselation,geometry,
    statistics,hypermesh,deleteFiles):
    #####To create input file for SpherePack
    if packing:
        myfile=os.path.join(mypath,'Project01.prj')
        f=open(myfile,'w')
        f.write('{0:34s}\n'.format('# RANDOM COORDINATES: 0 OR FILE: 1'))
        f.write('{0}\n'.format('0'))
        f.write('{0}\n'.format('# NUMBER OF SPHERES'))
        f.write('{0:d}\n'.format(NumOfCells))
        f.write('{0}\n'.format('# LENGTH IN X-, Y- AND Z-DIRECTION'))
        f.write('{0:f}\n'.format(4.0))
        f.write('{0:f}\n'.format(4.0))
        f.write('{0:f}\n'.format(4.0))
        f.write('{0}\n'.format('# EPSILON'))
        f.write('{0}\n'.format('0.0001'))
        f.write('{0}\n'.format('# NTAU'))
        f.write('{0}\n'.format('963800000'))
        f.write('{0}\n'.format('# NOMINAL DENSITY'))
        f.write('{0}\n'.format('0.25'))
        f.write('{0}\n'.format('# MAX. AND MINIM. DIAMETER'))
        f.write('{0}\n'.format('3.0'))
        f.write('{0}\n'.format('1.0'))
        f.write('{0}\n'.format('# MAX. NUMBER OF STEPS'))
        f.write('{0}\n'.format('50000000'))
        f.write('{0}\n'.format('# DISTRIBUTION'))
        f.write('{0}\n'.format('3'))
        f.write('{0}\n'.format('# DISTRIBUTION PARAMETERS'))
        f.write('{0:f}\n'.format(MU))
        f.write('{0:f}\n'.format(SIGMA))
        f.write('{0}\n'.format('0.4'))
        f.write('{0}\n'.format('-3.3'))
        f.write('{0}\n'.format('3'))
        f.write('{0}\n'.format('0.10'))
        f.write('{0}\n'.format('0.20'))
        f.write('{0}\n'.format('0.70'))
        f.write('{0}\n'.format('0.902113'))
        f.write('{0}\n'.format('3.000000'))
        f.write('{0}\n'.format('1.000000'))
        f.close()
        ###########################################################
        os.chdir(mypath)
        command1='wine SpherePackFB.exe'
        os.system(command1)
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
        Centers27=np.array([X,X+L2,X-L2,X+L2,X,X+L2,X,X-L2,X-L2,X-L2,X+L2,X,X,X-L2,X+L2,X+L2,X,X,X-L2,X,X,X-L2,X+L2,X+L2,X-L2,X-L2,X+L2,Y,Y+L2,Y-L2,Y+L2,Y+L2,Y,Y-L2,Y,Y-L2,Y+L2,Y-L2,Y-L2,Y+L2,Y,Y,Y,Y+L2,Y,Y,Y-L2,Y,Y+L2,Y-L2,Y+L2,Y-L2,Y+L2,Y-L2,Z,Z+L2,Z-L2,Z,Z+L2,Z+L2,Z-L2,Z-L2,Z,Z,Z,Z+L2,Z-L2,Z+L2,Z-L2,Z,Z,Z+L2,Z,Z,Z-L2,Z+L2,Z+L2,Z-L2,Z+L2,Z-L2,Z-L2])
        # Creation of input .txt files for neper
        myfile3=os.path.join(mypath,'Centers.txt')
        ff=open(myfile3,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(Centers27[i,j],Centers27[i+27,j],Centers27[i+54,j]))
        ff.close()
        #mypath4=r'/home/mohammad/'
        myfile4=os.path.join(mypath,'Rads.txt')
        fff=open(myfile4,'w')
        for i in range(0,27):
            for j in range(0,NumOfCells):
                fff.write('{0:f}\n'.format(Rads[j]))
        fff.close()
        commandTessellation="neper -T -n {0:d} -domain 'cube(12,12,12)' -morpho @Centers.txt -weight @Rads.txt -o RVE27 -format geo -statcell vol -statedge length -statface area -statver x".format((27*NumOfCells))
        os.system(commandTessellation)
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
        ####################################################
        # Creating periodic domain - not working yet
        # minDomain=[4,4,4]
        # maxDomain=[8,8,8]
        # def inDomain(point):
        #     inDomain=1
        #     for i,pos in enumerate(point):
        #         if (pos<minDomain[i] or pos>maxDomain[i]):
        #             inDomain=0
        #             return inDomain
        #     return inDomain
        # NodesInDomain=[]
        # NodeIndex=[]
        # NodeIndexNew=range(len(Nodes))
        # for i,point in enumerate(Nodes):
        #     if (inDomain(point)):
        #         NodesInDomain.append(point)
        #         NodeIndex.append(i)
        #         NodeIndexNew[i]=len(NodeIndex)-1
        # EdgesInDomain=[]
        # for i,index in enumerate(Edges):
        #     if (index[0] in NodeIndex or index[1] in NodeIndex):
        #         EdgesInDomain.append([NodeIndexNew[int(index[0])],NodeIndexNew[int(index[1])]])
        # geofile=open('PeriodicRVEdomain.geo','w')
        # for i,pos in enumerate(NodesInDomain):
        #     geofile.write('Point ({0}) = {{{1},{2},{3}}};\n'.format(i,pos[0],pos[1],pos[2]))
        # for i,index in enumerate(EdgesInDomain):
        #     geofile.write('Line ({0}) = {{{1},{2}}};\n'.format(i,index[0],index[1]))
        # geofile.close()
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
            HyperMeshInput.write('{0}{1}{2}{3}{4}{5}{6}\n'.format('*createnode(',Nodes[i][0],',',Nodes[i][1],',',Nodes[i][2],',0,0,0)'))
        for i in range(int(MaxIndexOfEdges)):
            HyperMeshInput.write('{0}\n{1}{2}{3}{4}\n'.format('*linecreatefromnodes(1,0,150,5,179)','*createlist(nodes,1) ',int(Edges[i][0]),' ',int(Edges[i][1])))
        for i in range(int(MaxIndexOfFaces)):
            HyperMeshInput.write('{0}\n'.format('*surfacemode(4)'))
            HyperMeshInput.write('{0}'.format('*createmark(lines,1) '))
            for j in range(len(Faces[i])):
                HyperMeshInput.write('{0}\t'.format(int(Faces[i][j])))
            HyperMeshInput.write('\n{0}\n{1}\n'.format('*createplane(1,1,0,0,0,0,0)','*splinesurface(lines,1,1,1,1)'))
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
            HyperMeshInput.write('{0}\n'.format('*createmark(surfaces,1) "displayed"'))
            HyperMeshInput.write('{0}{1}{2}{3}{4}\n'.format('*selfstitchcombine(1,82,',Tol,',',Tol,')'))
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
