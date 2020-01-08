"""
/***************************************************************************
Soil parameters selection
        begin                : 2018-09-24
        copyright            : (C) 2018 by Giacomo Titti and Giulia Bossi,
                               CNR-IRPI, Padova
        email                : giacomo.titti@irpi.cnr.it
 ***************************************************************************/

/***************************************************************************
    Soil parameters selection
    Copyright (C) 2018 by Giacomo Titti and Giulia Bossi CNR-IRPI, Padova

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
from main import TerrainFinder

test = TerrainFinder()
##############################
test.pathRoot='/to/folder/test'#head input folder
test.pathHome='/to/work/folder'#work folder
test.pathGrid='/to/test/grid.csv'# grid file from FLAC3D
test.pathCentroids='/to/test/centroid.csv'# centroids of grid file from FLAD3D
test.pathDEM='/to/test/DEM.tif'# DEM used for numerical modeling
test.pathGPS0='/to/tets/GPS0'# GPS benchmarks location (x,y,z) at time 0, file xxx.txt
test.pathGPS1='/to/test/GPS1'# GPS benchmarks location (x,y,z) at time 1, file xxx.txt
test.pathInclinometer='/to/test/Inclinometer' #Inclinmeters location depth and displacement (x,y,z,depth,dispx,dispy), file xxx.txt
test.pathData='/to/test/perm'#log files of modeled displacements from FLAC3D
test.v2step=10#numerosity of top solutions
test.stringsG=['001'.'002']# ['001', '002', '003',...]
test.stringsI=['001','002']# ['001', '002', '003',...]

##############################
test.grid_construction()
test.GPS_inclinometers()
print('end')
