#!/usr/bin/env python
from __future__ import print_function
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""Display sky map geometry as a 3d plot
"""
import math
import numpy
import sys
import pickle

from mpl_toolkits.mplot3d import Axes3D # used by fig.gca
import matplotlib.pyplot as plt

import lsst.afw.geom as afwGeom

def reportSkyMapInfo(skyMap):
    paramDict = skyMap.config.toDict()
    paramNameList = sorted(paramDict.iterkeys())
    print("Sky Map parameters:")
    for paramName in paramNameList:
        param = paramDict[paramName]
        print("skyMap.config.%s = %s" % (paramName, param))

    print("\nSkyMap has %d tracts:" % (len(skyMap)))
    for tractInfo in skyMap:
        wcs = tractInfo.getWcs()
        posBox = afwGeom.Box2D(tractInfo.getBBox())
        pixelPosList = (
            posBox.getMin(),
            afwGeom.Point2D(posBox.getMaxX(), posBox.getMinY()),
            posBox.getMax(),
            afwGeom.Point2D(posBox.getMinX(), posBox.getMaxY()),
        )
        skyPosList = [wcs.pixelToSky(pos).getPosition(afwGeom.degrees) for pos in pixelPosList]
        posStrList = ["(%0.3f, %0.3f)" % tuple(skyPos) for skyPos in skyPosList]
        print("tract %s has corners %s (RA, Dec deg) and %s x %s patches" % \
            (tractInfo.getId(), ", ".join(posStrList), \
            tractInfo.getNumPatches()[0], tractInfo.getNumPatches()[1]))

def plotSkyMap(skyMap):
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_axis_off()

    # make sure a complete 1x1x1 cube is shown -- what I really want is to constrain the aspect ratio
    # but that is not yet supported for 3D plots
    for direction in (-1, 1):
        for point in numpy.diag(direction * numpy.array([1,1,1])):
            ax.plot([point[0]], [point[1]], [point[2]], 'w')    

    for tractInfo in skyMap:
        # display outer edge; scale to be approximately in the same plane as the inner region
        wcs = tractInfo.getWcs()
        posBox = afwGeom.Box2D(tractInfo.getBBox())
        xRange = posBox.getMinX(), posBox.getMaxX()
        yRange = posBox.getMinY(), posBox.getMaxY()
        
        numX = min(50, max(1, ((xRange[1] - xRange[0]) // 100)))
        numY = min(50, max(1, ((yRange[1] - yRange[0]) // 100)))
        
        outerPixPosList = \
              [(x, yRange[0]) for x in numpy.linspace(xRange[0], xRange[1], num=numX, endpoint=False)] \
            + [(xRange[1], y) for y in numpy.linspace(yRange[0], yRange[1], num=numY, endpoint=False)] \
            + [(x, yRange[1]) for x in numpy.linspace(xRange[1], xRange[0], num=numX, endpoint=False)] \
            + [(xRange[0], y) for y in numpy.linspace(yRange[1], yRange[0], num=numY, endpoint=False)]
        outerPixPosList.append(outerPixPosList[0])
        
        outerPoints = [numpy.array(wcs.pixelToSky(p[0], p[1]).getVector()) for p in outerPixPosList]
        outX, outY, outZ = zip(*outerPoints)
        ax.plot(outX, outY, outZ)
    
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("To use: plotSkyMap path-to-skymap-pickle")
        sys.exit(1)
    
    pathToSkyMap = sys.argv[1]
    
    with file(pathToSkyMap, "r") as f:
        skyMap = pickle.load(f)
        reportSkyMapInfo(skyMap)
        plotSkyMap(skyMap)
