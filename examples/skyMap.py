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

_RadPerDeg = math.pi / 180.0

for overlapDeg in (0.0, 0.33, 1.0, 3.5):
    print "overlap = %s degrees" % (overlapDeg)
    skyMap = lsst.skymap.SkyMap(overlap=afwGeom.Angle(overlapDeg, afwGeom.degrees))
    totNumPix = 0
    print "Tract  Ctr RA  Ctr Dec    Rows        Cols       # Pix   Width  Height"
    print " ID     (deg)   (deg)     (pix)       (pix)              (deg)  (deg)"
    for i in range(12):
        skyTractInfo = skyMap.getSkyTractInfo(i)
        bbox = skyTractInfo.getBBox()
        dimensions = bbox.getDimensions()
        numPix = dimensions[0] * dimensions[1]
        totNumPix += numPix
        wcs = skyTractInfo.getWcs()
        posBBox = afwGeom.Box2D(bbox)
        ctrPixPos = posBBox.getCenter()
        ctrCoord = wcs.pixelToSky(ctrPixPos)
        ctrSkyPosDeg = ctrCoord.getPosition(afwGeom.degrees)
        leftCoord   = wcs.pixelToSky(afwGeom.Point2D(posBBox.getMinX(), ctrPixPos[1]))
        rightCoord  = wcs.pixelToSky(afwGeom.Point2D(posBBox.getMaxX(), ctrPixPos[1]))
        topCoord    = wcs.pixelToSky(afwGeom.Point2D(ctrPixPos[0], posBBox.getMinY()))
        bottomCoord = wcs.pixelToSky(afwGeom.Point2D(ctrPixPos[0], posBBox.getMaxY()))
        xSpan = leftCoord.angularSeparation(rightCoord).asDegrees()
        ySpan = bottomCoord.angularSeparation(topCoord).asDegrees()
        print "%3d   %7.1f %7.1f %10.2e  %10.2e %10.1e %6.1f %6.1f" % \
            (i, ctrSkyPosDeg[0], ctrSkyPosDeg[1], dimensions[0], dimensions[1], numPix, xSpan, ySpan)
    
    nomPixelArea = skyMap.getPixelScale().asRadians()**2 # nominal area of a pixel in rad^2
    numPixToTileSphere = 4 * math.pi / nomPixelArea
    print "total # pixels = %.1e\npixels to tile sphere = %.1e\nextra storage (tot pix/pix to tile) = %.1f\n" % \
        (totNumPix, numPixToTileSphere, totNumPix / float(numPixToTileSphere))
