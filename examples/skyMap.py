#!/usr/bin/env python
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
import math
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.skymap

for overlapDeg in (0.0001, 0.33, 1.0, 3.5):
    print "overlap = %.1f degrees" % (overlapDeg)
    skyMap = lsst.skymap.SkyMap(overlap=overlapDeg * math.pi / 180.0)
    totNumPix = 0
    for i in range(12):
        skyTileInfo = skyMap.getSkyTileInfo(i)
        ctrCoord = skyTileInfo.getCtrCoord()
        ctrPos = ctrCoord.getPosition(afwCoord.DEGREES)
        dimensions = skyTileInfo.getDimensions()
        numPix = dimensions[0] * dimensions[1]
        totNumPix += numPix
        wcs = skyTileInfo.getWcs()
        leftCoord = wcs.pixelToSky(0.0, dimensions[1]/2.0)
        rightCoord = wcs.pixelToSky(dimensions[0]-1, dimensions[1]/2.0)
        topCoord = wcs.pixelToSky(dimensions[0]/2.0, dimensions[1] - 1)
        bottomCoord = wcs.pixelToSky(dimensions[0]/2.0, 0)
        xSpan = leftCoord.angularSeparation(rightCoord, afwCoord.DEGREES)
        ySpan = bottomCoord.angularSeparation(topCoord, afwCoord.DEGREES)
        print "sky tile %2d center RA/Dec = %6.1f, %6.1f; dimensions = %.2e x %.2e pixels (%0.1e total); span = %.1f x %.1f deg" % \
            (i, ctrPos[0], ctrPos[1], dimensions[0], dimensions[1], numPix, xSpan, ySpan)
    
    nomPixelArea = skyMap.getPixelScale()**2 # nominal area of a pixel in deg^2
    numPixToTileSphere = 4 * math.pi * (180.0 / math.pi)**2 / nomPixelArea
    print "total # pixels = %.1e, pixels to tile sphere = %.1e, extra storage (tot pix/pix to tile) = %.1f" % \
        (totNumPix, numPixToTileSphere, totNumPix / float(numPixToTileSphere))
