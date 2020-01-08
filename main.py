"""
/***************************************************************************
Soil parameters selection
        begin                : 2018-09-24
        copyright            : (C) 2018 by Giacomo Titti and Giulia Bossi, CNR-IRPI, Padova
        email                : giacomo.titti@irpi.cnr.it
 ***************************************************************************/

/***************************************************************************
    Soil parameters selection
    Copyright (C) 2018 by Giacomo Titti and Giulia Bossi, CNR-IRPI, Padova

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/
"""
#coding=utf-8
import sys
import numpy as np
import gdal
import os
import math
import re
import csv
import subprocess

class TerrainFinder:
    def __init__(self):
        self.pathRoot=''
        self.pathHome=''
        self.pathGrid=''
        self.pathCentroids=''
        self.pathDEM=''
        self.pathLow=self.pathHome+'/cubesUp.txt'
        self.pathUp=self.pathHome+'/cubesLow.txt'
        self.pathcsvfile1step='/solution1step.txt'
        self.pathcsvfile2step='/solution2step.txt'
        self.variables={}
        self.variables['v1step']=1
        self.variables['v2step']=5
        test.pathGPS0=''
        test.pathGPS1=''
        test.pathInclinometer=''
        self.stringsG=[]
        self.stringsI=[]

    def grid_construction(self):
        self.grid = np.array(list(csv.reader(open(self.pathGrid, "rb"), delimiter=" "))).astype("float")
        self.centroids = np.array(list(csv.reader(open(self.pathCentroids, "rb"), delimiter=" "))).astype("float")
        landslide=self.centroids
        landslide=np.delete(landslide, 1, axis=1)
        cubesLow=np.zeros((len(landslide),2,2),dtype=int)
        cubesUp=np.zeros((len(landslide),2,2),dtype=int)
        count=0
        for row in landslide:
            low=np.zeros((2,2),dtype=int)
            up=np.zeros((2,2),dtype=int)
            eucledian=np.array([])
            cube=np.array([])
            cubeUp=np.array([])
            cubeLow=np.array([])
            cubeUpNorth=np.array([])
            cubeUpSouth=np.array([])
            cubeLowNorth=np.array([])
            cubeLowSouth=np.array([])
            cubeUpNorthWest=np.array([])
            cubeUpNorthEast=np.array([])
            cubeUpSouthWest=np.array([])
            cubeUpSouthEast=np.array([])
            cubeLowNorthWest=np.array([])
            cubeLowNorthEast=np.array([])
            cubeLowSouthWest=np.array([])
            cubeLowSouthEast=np.array([])
            for rows in self.grid:
                dist=np.float(math.sqrt(np.sum((row[1:]-rows[1:])**2)))
                eucledian=np.append(eucledian,dist)
            a=np.vstack((eucledian,self.grid.T)).T
            a=a[np.argsort(a[:,0])]
            cube=a[:8,1:]
            cube=cube[np.argsort(cube[:,3])]
            ####
            cubeLow=cube[:4,:]
            cubeUp=cube[4:,:]
            cubeLow=cubeLow[np.argsort(cubeLow[:,2])]
            cubeUp=cubeUp[np.argsort(cubeUp[:,2])]
            ####
            cubeLowSouth=cubeLow[:2,:]
            cubeLowNorth=cubeLow[2:,:]
            cubeUpSouth=cubeUp[:2,:]
            cubeUpNorth=cubeUp[2:,:]
            cubeLowSouth=cubeLowSouth[np.argsort(cubeLowSouth[:,1])]
            cubeLowNorth=cubeLowNorth[np.argsort(cubeLowNorth[:,1])]
            cubeUpSouth=cubeUpSouth[np.argsort(cubeUpSouth[:,1])]
            cubeUpNorth=cubeUpNorth[np.argsort(cubeUpNorth[:,1])]
            ####
            cubeLowSouthWest=cubeLowSouth[0,:]
            cubeLowSouthEast=cubeLowSouth[1,:]
            cubeLowNorthWest=cubeLowNorth[0,:]
            cubeLowNorthEast=cubeLowNorth[1,:]
            cubeUpSouthWest=cubeUpSouth[0,:]
            cubeUpSouthEast=cubeUpSouth[1,:]
            cubeUpNorthWest=cubeUpNorth[0,:]
            cubeUpNorthEast=cubeUpNorth[1,:]
            ####
            low=np.array([[cubeLowNorthWest[0],cubeLowNorthEast[0]],[cubeLowSouthWest[0],cubeLowSouthEast[0]]]).astype(int)
            up=np.array([[cubeUpNorthWest[0],cubeUpNorthEast[0]],[cubeUpSouthWest[0],cubeUpSouthEast[0]]]).astype(int)
            cubesLow[count,:,:]=low
            cubesUp[count,:,:]=up
            count=count+1
        with open(self.pathLow,'wb') as f:
            for a in cubesLow:
               np.savetxt(f, a, fmt='%10d')
               f.write('\n')
        with open(self.pathUp,'wb') as f:
            for a in cubesUp:
               np.savetxt(f, a, fmt='%10d')
               f.write('\n')

    def input_data(self):
        self.PM1=len(self.stringsG)
        self.PM2=len(self.stringsI)
        self.dem=self.pathDEM
        self.grid = np.array(list(csv.reader(open(self.pathGrid, "rb"), delimiter=" "))).astype("float")
        self.centroids = np.array(list(csv.reader(open(self.pathCentroids, "rb"), delimiter=" "))).astype("float")
        f=np.loadtxt(self.pathUp).astype(int)
        self.cubesL=np.zeros((len(f)/2,2,2),dtype=int)
        countf=0
        for i in range(0,len(f),2):
            self.cubesL[countf,:,:]=f[i:i+2,:]
            countf=countf+1
        ff=np.loadtxt(self.pathLow).astype(int)
        self.cubesU=np.zeros((len(ff)/2,2,2),dtype=int)
        countff=0
        for i in range(0,len(ff),2):
            self.cubesU[countff,:,:]=ff[i:i+2,:]
            countff=countff+1

    def CentroidPixel(self):
        src_ds=gdal.Open(self.dem)
        gt=src_ds.GetGeoTransform()
        rb=src_ds.GetRasterBand(1)
        arr = np.float16(rb.ReadAsArray())
        [rows, columns]=np.shape(arr)#dem dimensions
        gdal.UseExceptions()#it doesn't print to screen everytime point is outside grid
        mx=None
        my=None
        self.q=None
        mx=self.first[0]
        my=self.first[1]
        self.q=self.first[2]
        self.LowerLeft=np.array([gt[0],gt[3]-(rows*gt[1])])
        self.LowerLeft[0]=self.LowerLeft[0]
        self.LowerLeft[1]=self.LowerLeft[1]
        self.cx=None
        self.cy=None
        self.cx=mx-self.LowerLeft[0]
        self.cy=my-self.LowerLeft[1]
        src_ds=None

    def NearestCentroidFlug3d(self):
        self.ref=None
        self.ref=10000
        id=None
        id=1
        count=0
        for row in self.centroids:
            if row[1]==1 or row[1]==2:
                diff=None
                vect=None
                vect=math.sqrt((self.cx-row[2])**2+(self.cy-row[3])**2+(self.q-row[4])**2)
                if vect < self.ref:
                    self.ref=None
                    self.ref=vect
                    self.id_ref=None
                    self.id_ref=row[0]
                    row_ref=np.array([])
                    row_ref=row
                    diff_ref=None
                    diff_ref=diff
                    self.xyz=np.array([])
                    self.xyz=np.array([row[2], row[3], row[4]])
                rows=np.array([row[0],row[2], row[3], row[4]])
                count=count+1
            id=id+1
        if row_ref[1]==0:
            print('wrong point of reference!!!!')
        elif row_ref[1]==1:
            print('landslide n°1')
        else:
            print('landslide n°2')

    def PolarCoordinatesflac(self):
        #print(lineflac)
        if self.lineflac[1]==0 or self.lineflac[2]==0:
            b=0
        else:
            b=abs(math.atan(float(self.lineflac[2])/float(self.lineflac[1])))
        #print(lineflac,b)
        if self.lineflac[1]>0 and self.lineflac[2]>0:
            self.gamma=math.pi/2.-b
        elif self.lineflac[1]>0 and self.lineflac[2]<0:
            self.gamma=math.pi-(math.pi/2.-b)
        elif self.lineflac[1]<0 and self.lineflac[2]<0:
            self.gamma=(3./2.)*math.pi-b
        elif self.lineflac[1]<0 and self.lineflac[2]>0:
            self.gamma=2.*math.pi-((math.pi/2.)-b)
        elif self.lineflac[1]>0 and self.lineflac[2]==0:
            self.gamma=math.pi/2.
        elif self.lineflac[1]==0 and self.lineflac[2]<0:
            self.gamma=math.pi
        elif self.lineflac[1]<0 and self.lineflac[2]==0:
            self.gamma=(3./2.)*math.pi
        elif self.lineflac[1]==0 and self.lineflac[2]>0:
            self.gamma=0
        elif self.lineflac[1]==0 and self.lineflac[2]==0:
            self.gamma=None
        xy=math.sqrt((math.pow(self.lineflac[1], 2))+(math.pow(self.lineflac[2], 2)))
        if xy==0 or self.lineflac[3]==0:
            z=0
        else:
            z=abs(math.atan(float(self.lineflac[3])/float(xy)))
        ####
        if xy>0 and self.lineflac[3]>0:
            self.omega=math.pi/2.-z
        elif xy>0 and self.lineflac[3]<0:
            self.omega=math.pi-(math.pi/2.-z)
        elif xy<0 and self.lineflac[3]<0:
            self.omega=(3./2.)*math.pi-z
        elif xy<0 and self.lineflac[3]>0:
            self.omega=2.*math.pi-(math.pi/2.-z)
        elif xy>0 and self.lineflac[3]==0:
            self.omega=math.pi/2.
        elif xy==0 and self.lineflac[3]<0:
            self.omega=math.pi
        elif xy<0 and self.lineflac[3]==0:
            self.omega=(3./2.)*math.pi
        elif xy==0 and self.lineflac[3]>0:
            self.omega=0
        elif xy==0 and self.lineflac[3]==0:
            self.omega=None

    def PolarCoordinates(self):
        if self.line[0]==0 or self.line[1]==0:
            b=0
        else:
            b=abs(math.atan(float(self.line[1])/float(self.line[0])))
        if self.line[0]>0 and self.line[1]>0:
            self.alfa=math.pi/2.-b
        elif self.line[0]>0 and self.line[1]<0:
            self.alfa=math.pi-(math.pi/2.-b)
        elif self.line[0]<0 and self.line[1]<0:
            self.alfa=(3./2.*np.pi)-b
        elif self.line[0]<0 and self.line[1]>0:
            self.alfa=2.*math.pi-(math.pi/2.-b)
        elif self.line[0]>0 and self.line[1]==0:
            self.alfa=math.pi/2.
        elif self.line[0]==0 and self.line[1]<0:
            self.alfa=math.pi
        elif self.line[0]<0 and self.line[1]==0:
            self.alfa=(3./2.)*math.pi
        elif self.line[0]==0 and self.line[1]>0:
            self.alfa=0
        elif self.line[0]==0 and self.line[1]==0:
            self.alfa=None
        ####
        if len(self.line)==3:
            xy=math.sqrt((math.pow(self.line[0], 2))+(math.pow(self.line[1], 2)))
            if xy==0 or self.line[2]==0:
                z=0
            else:
                z=abs(math.atan(float(self.line[2])/float(xy)))
            if xy>0 and self.line[2]>0:
                self.beta=math.pi/2.-z
            elif xy>0 and self.line[2]<0:
                self.beta=math.pi-(math.pi/2.-z)
            elif xy<0 and self.line[2]<0:
                self.beta=(3./2.)*math.pi-z
            elif xy<0 and self.line[2]>0:
                self.beta=2.*math.pi-(math.pi/2.-z)
            elif xy>0 and self.line[2]==0:
                self.beta=math.pi/2.
            elif xy==0 and self.line[2]<0:
                self.beta=math.pi
            elif xy<0 and self.line[2]==0:
                self.beta=(3./2.)*math.pi
            elif xy==0 and self.line[2]>0:
                self.beta=0
            elif xy==0 and self.line[2]==0:
                self.beta=None
        else:
            self.beta=0
    def ReadPerm(self):
        path=os.path.join(self.pathData)
        perm=np.zeros(3,dtype=float)
        perm.fill(100000)
        ListPerm=np.array([(0.,0.,0.,0.,0.,'txt')],dtype='float,float,float,float,float,object')#id,ratio,pl,ver,path
        List2step=np.array([(0.,0.,0.,0.,0.,0.,'txt')],dtype='float,float,float,float,float,float,object')
        AngleErrors=np.zeros((6,2),dtype=float)
        RatioErrors=np.zeros((3,3),dtype=float)
        cc=0
        for filename in os.listdir(path):
            c=0
            if filename.endswith(".log"):
                filenameL=np.array([filename])
                filenamelist=filenameL
                if np.any(self.name==filenamelist):
                    print('next perm')
                else:
                    f = open((os.path.join(path, filename)),"r")
                    lines = f.readlines()
                    CheckAngles=np.zeros((len(self.id), 2), dtype=float)
                    CheckDisp=np.zeros(len(self.id), dtype=float)
                    dispGPS=np.array([])
                    for ii in range(len(self.id)):
                        idL=np.reshape(self.cubesL[self.id[ii]-1,:,:], (1,4)).astype(int)
                        idU=np.reshape(self.cubesU[self.id[ii]-1,:,:], (1,4)).astype(int)
                        gammaV=np.array([])
                        omegaV=np.array([])
                        displacement=np.array([])
                        displaceI=np.array([])
                        for iiii in range(4):
                            #lower quadrant options
                            x=lines[idL[0,iiii]+10-1]
                            List=re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", x)
                            self.lineflac=np.array([])
                            self.lineflac=(np.array([float(i) for i in List]))
                            displace=float(math.sqrt((self.lineflac[1])**2+(self.lineflac[2])**2+(self.lineflac[3])**2))
                            displacement=np.append(displacement,displace)
                            self.PolarCoordinatesflac()
                            gammaV=np.append(gammaV,self.gamma)
                            omegaV=np.append(omegaV,self.omega)
                            #upper quadrant options
                            y=lines[idU[0,iiii]+10-1]
                            List=re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", y)
                            self.lineflac=np.array([])
                            self.lineflac=(np.array([float(i) for i in List]))
                            displace=float(math.sqrt((self.lineflac[1])**2+(self.lineflac[2])**2+(self.lineflac[3])**2))
                            displacement=np.append(displacement,displace)
                            self.PolarCoordinatesflac()
                            gammaV=np.append(gammaV,self.gamma)
                            omegaV=np.append(omegaV,self.omega)
                            ####
                        sumPlane=float((np.sum(gammaV*displacement))/(np.sum(displacement)))#weighted average
                        sumVert=float((np.sum(omegaV*displacement))/(np.sum(displacement)))#weighted average
                        dispT=float((np.sum(displacement))/(len(displacement)))#weighted average
                        CheckDisp[ii]=dispT
                        CheckAngles[ii, 0]=sumPlane
                        CheckAngles[ii, 1]=sumVert
                        if ii<=len(self.idGPS)-1:
                            modeledDisp=np.sum(displacement)/len(displacement)
                            dispGPS=np.append(dispGPS,modeledDisp)
                    div=np.zeros((len(self.idGPS),len(self.idGPS)),dtype=float)
                    divMeasured=np.zeros((len(self.idGPS),len(self.idGPS)),dtype=float)
                    for s in range(len(self.idGPS)):
                        for ss in range(len(self.idGPS)):
                            div[ss,s]=float(dispGPS[s]/dispGPS[ss])
                            divMeasured[ss,s]=float(self.dispGPSmeasured[s]/self.dispGPSmeasured[ss])
                    fitDiv=(abs(divMeasured-div))
                    ##########removing None value where the points are not moving
                    fit=np.array([])
                    fit1=np.array([])
                    fitt=np.zeros(4)
                    idx=None
                    idx=np.argwhere(np.isnan(CheckAngles))
                    fit=abs(self.mAngle-CheckAngles)
                    if len(idx)>0:
                        print('some points are not valid!!!!!!!!',stop)
                    fitPlaneGps=fit[:self.PM1,0]#planar errors of GPS
                    fitPlaneInc=fit[self.PM1:,0]#planar errors of Invlinometers
                    fitVertGps=fit[:self.PM1,1]#vertical errors of GPS
                    [row,col]=np.shape(fitDiv)
                    dim=row*col
                    ##########################################################
                    #based on average errors
                    #fitPlaneGps=np.array([0])#if we are not considering xy plane of GPS
                    #fitPlaneInc=np.array([0])#if we are not considering xy plane of Inclinometers
                    fitVertGps=np.array([0])#if we are not considering vertical plane GPS
                    #fitDiv=np.array([0])###if we are not considering ratio
                    ##################################################average of the parameters
                    fitt[0]=np.sum((fitPlaneGps/float(np.pi*2))/float(len(fitPlaneGps)))
                    fitt[1]=np.sum((fitPlaneInc/float(np.pi*2))/float(len(fitPlaneInc)))
                    fitt[2]=np.sum((fitVertGps/float(np.pi*2))/float(len(fitVertGps)))
                    fitt[3]=float(np.sum(fitDiv)/dim)
                    ###########################################################
                    #collection of permutations
                    rowPerm=np.array([(cc,fitt[0],fitt[1],fitt[2],fitt[3],filename)], dtype='float,float,float,float,float,object')
                    row2step=np.array([(cc,float((fitt[0]+fitt[1])/2.),fitt[0],fitt[1],fitt[2],fitt[3],filename)], dtype='float,float,float,float,float,float,object')
                    ListPerm=np.append(ListPerm,rowPerm)
                    List2step=np.append(List2step,row2step)
                    AngleErrors=np.dstack((AngleErrors,fit))
                    RatioErrors=np.dstack((RatioErrors,fitDiv))
                    cc=cc+1
        #########################
        ListPerm=ListPerm[1:]#the first is zeros()
        List2step=List2step[1:]#the first is zeros()
        AngleErrors=AngleErrors[:,:,1:]#the first is zeros()
        RatioErrors=RatioErrors[:,:,1:]#the first is zeros()
        ######################
        self.Perm10=[]
        self.Perm1=[]
        ee=np.zeros((1,2),dtype=int)
        eee=np.array([])
        for n in range(len(List2step)):
            aa=List2step[n][0]#count
            ff=List2step[n][1]#h
            dd=List2step[n][5]#r
            eee=np.array([aa,int(ff*1000000)])
            ee=np.vstack((ee,eee))
        ee=ee[1:,:]
        ee=sorted(ee, key=lambda a_entry: a_entry[1])
        ee=np.asarray(ee)
        col=ee[:,1]
        p=int(np.percentile(col, 25))
        idx=np.where(col>=p)
        idx2=int(ee[idx[0][0],0])
        ListPermSorted=sorted(List2step, key=lambda b_entry: b_entry[self.variables['v1step']])
        [item for item in ListPermSorted if item[0] == idx2]
        idx3=[x for x, y in enumerate(ListPermSorted) if y[0] == idx2]
        self.Perm10=ListPermSorted[:idx3[0]]
        ListPerm10=sorted(self.Perm10, key=lambda c_entry: c_entry[self.variables['v2step']])
        self.Perm1=ListPerm10[:self.v2step]#best
        #########################

        #################core functions###########################################
    def GPS(self):
        print('GPS')
        #################################GPS
        self.idGPS=np.zeros((self.PM1), dtype=np.int)
        self.dispGPSmeasured=np.zeros((self.PM1), dtype=float)
        self.mAngleGPS=np.zeros((self.PM1, 2))
        coord=np.array([])
        positionXYZ=np.array([])
        for i in range(self.PM1):
            #open txt file of GPS coordinates
            #fr=['fr1', 'fr2', 'fr2']
            GPS0=os.path.join(self.pathGPS0+'/xxx.txt')
            GPS2=GPS0.replace('xxx', self.stringsG[i])
            #GPS2=GPS1.replace('nnn', fr[i])
            self.first=None
            self.first=np.loadtxt(GPS2)
            positionXYZ=np.append(positionXYZ, self.first)
            GPSfirst0=os.path.join(self.pathGPS0+'/xxx.txt')
            GPSfirst2=GPSfirst0.replace('xxx', self.stringsG[i])
            #GPSfirst2=GPSfirst1.replace('nnn', fr[i])
            self.GPSfirst=np.loadtxt(GPSfirst2)
            GPSsecond0=os.path.join(self.pathGPS1+'/xxx.txt')
            GPSsecond2=GPSsecond0.replace('xxx', self.stringsG[i])
            #GPSsecond2=GPSsecond1.replace('nnn', fr[i])
            self.GPSsecond=np.loadtxt(GPSsecond2)
            self.CentroidPixel()
            self.NearestCentroidFlug3d()
            coord=np.append(coord, self.xyz)
            print('Distance from centroid of reference: %s' %self.ref, 'id centroid of reference: %s' %self.id_ref.astype(int))
            self.line=np.array([])
            self.line=self.GPSsecond-self.GPSfirst
            self.PolarCoordinates()#GPS polar coordinates
            self.mAngleGPS[i, :]=np.array([self.alfa, self.beta])
            self.idGPS[i]=self.id_ref
            self.dispGPSmeasured[i]=float(math.sqrt((self.line[0])**2+(self.line[1])**2+(self.line[2])**2))
        position0=np.reshape(positionXYZ, (self.PM1, 3))
        points=np.reshape(coord, (self.PM1, 3))

    def inclinometers(self):
        print('inclinometers')
        #################################inclinometers
        self.idInc=np.zeros((self.PM2), dtype=np.int)
        self.mAngleInc=np.zeros((self.PM2, 2))
        coord=np.array([])
        positionXYZ=np.array([])
        for i in range(self.PM2):
            self.first=None
            #open txt file of inclinometers coordinates
            #fr=['fr1', 'fr2', 'fr2']
            first0=os.path.join(self.pathInclinometer+'/xxx.txt')
            first2=first0.replace('xxx', self.stringsI[i])
            #first2=first1.replace('nnn', fr[i])
            self.first=np.loadtxt(first2)
            positionXYZ=np.append(positionXYZ, self.first)
            #depth0=os.path.join(self.pathRoot+'/depth_slip_surf_per_inclinometer')
            #depth2=depth0.replace('xxx', self.stringsI[i])
            #depth2=depth1.replace('nnn', fr[i])
            depth=self.first[3]
            self.first[2]=self.first[2]-depth
            self.CentroidPixel()
            self.q=self.q#reference system altitude
            self.NearestCentroidFlug3d()
            coord=np.append(coord, self.xyz)
            print('Distance from centroid of reference: %s' %self.ref, 'id centroid of reference: %s' %self.id_ref.astype(int))
            #disp0=os.path.join(self.pathRoot+'/input_inclinometers')
            #disp2=disp0.replace('xxx', self.stringsI[i])
            #disp2=disp1.replace('nnn', fr[i])
            disp=np.loadtxt(disp2)
            self.line=np.array([self.first[4],self.first[5]])
            #self.line=self.first[3]
            self.PolarCoordinates()
            self.mAngleInc[i, :]=np.array([self.alfa, self.beta])
            self.idInc[i]=self.id_ref
        position0=np.reshape(positionXYZ, (self.PM2, 3))
        points=np.reshape(coord, (self.PM2, 3))######only plane xy

    def GPS_inclinometers(self):
        self.input_data()
        PMs=self.PM1+self.PM2
        self.GPS()
        self.inclinometers()
        self.id=np.append(self.idGPS,self.idInc)
        self.mAngle=np.vstack((self.mAngleGPS,self.mAngleInc))
        #self.name='perm'
        self.name=np.array([])
        self.ReadPerm()
        ########################top first quartile, first step
        self.csvfile1 =os.path.join(self.pathHome+self.pathcsvfile2step)
        with open(self.csvfile1, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in self.Perm1:
                writer.writerow([val])
        subprocess.check_call([self.pathHome+'/myscript.sh',os.path.join(self.pathHome+self.pathcsvfile2step)])
        ###########################top 10, second step
        self.csvfile10 =os.path.join(self.pathHome+self.pathcsvfile1step)
        with open(self.csvfile10, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in self.Perm10:
                writer.writerow([val])
        subprocess.check_call([self.pathHome+'/myscript.sh',os.path.join(self.pathHome+self.pathcsvfile1step)])
