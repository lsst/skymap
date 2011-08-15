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
import lsst.afw.image as afwImage

_RadPerDeg = math.pi / 180.0

class SkyTileInfo(object):
    def __init__(self, ind, ctrCoord, vertexCoordList, overlap, wcsFactory):
        """
        
        Inputs:
        - ind: index of sky tile
        - ctrCoord: Coord of center of sky tile
        - vertexCoordList: list of Coords of vertices that define edge of optimal area
        - overlap: minimum overlap between adjacent sky tiles (rad)
        - wcsFactory: a _WcsFactory object
        """
        
#         print "SkyTileInfo(ind=%s, ctrCoord=%s, overlap=%0.1f)" % \
#             (ind, ctrCoord.getPosition(afwCoord.DEGREES), overlap)
        self._ind = ind
        self._ctrCoord = ctrCoord
        self._vertexCoordList = vertexCoordList
        self._overlap = float(overlap)
        wcsFactory
        
        DebugMinInd = 0

        # start by computing everything relative to the center of the lower left pixel;
        # compute a WCS and use it to determine the required size of the tile
        basicWcs = wcsFactory.makeWcs(ctrInd=(0, 0), ctrCoord=self._ctrCoord)
        minPos = [0, 0]
        maxPos = [0, 0]
        
        # for now ignore overlap and just compute position at each vertex
        ctrSphPos = self._ctrCoord.getPosition(afwCoord.DEGREES)
        if ind < DebugMinInd:
            print "center position =", self._ctrCoord.getPosition(afwCoord.DEGREES)
        for vertexCoord in self._vertexCoordList:
            vertexDeg = vertexCoord.getPosition(afwCoord.DEGREES)
            if self._overlap == 0:
                numToTest = 1
            else:
                numToTest = 4
            for i in range(numToTest):
                offAngle = 90.0 * i
                offCoord = vertexCoord.clone()
                offCoord.offset(offAngle * _RadPerDeg, self._overlap / 2.0)
                pixPos = basicWcs.skyToPixel(offCoord)
                if ind < DebugMinInd:
                    print "vertexPos=%s, offPos=%s, xy=%s" % (vertexDeg, offCoord.getPosition(afwCoord.DEGREES), pixPos)
                for j, pixPosVal in enumerate(pixPos):
                    minPos[j] = min(minPos[j], pixPosVal)
                    maxPos[j] = max(maxPos[j], pixPosVal)
        minInd = [afwImage.positionToIndex(val) for val in minPos]
        maxInd = [afwImage.positionToIndex(val) for val in maxPos]
        if ind < DebugMinInd:
            print "minInd=%s, maxInd=%s" % (minInd, maxInd)
        
        self._dimensions = tuple(1 + maxInd[i] - minInd[i] for i in range(2))
        
        # now offset things so the lower left pixel is 0,0 and compute a new WCS using that
        self._ctrPixInd = tuple(-ind for ind in minInd)
        
        self._wcs = wcsFactory.makeWcs(ctrInd=self._ctrPixInd, ctrCoord=self._ctrCoord)

    def getCtrCoord(self):
        return self._ctrCoord

    def getCtrIndex(self):
        return self._ctrPixInd
    
    def getDimensions(self):
        return self._dimensions

    def getIndex(self):
        """Return sky tile index"""
        return self._ind
    
    def getWcs(self):
        return self._wcs
