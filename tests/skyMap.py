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

from __future__ import with_statement
"""Test SkyMap class
"""
import os
import sys
import math
import unittest

import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy
import lsst.pex.policy as pexPolicy
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.skymap as skymap

_RadPerDeg = math.pi / 180.0

# dodecahedron properties
_NumFaces = 12
_Phi = (1.0 + math.sqrt(5.0)) / 2.0
_DihedralAngle = 2.0 * math.atan(_Phi) / _RadPerDeg
_NeighborAngularSeparation = 180.0 - _DihedralAngle

class SkyMapTestCase(unittest.TestCase):
    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        sm = skymap.SkyMap()
        self.assertEquals(sm.getNumSkyTiles(), _NumFaces)
        self.assertEquals(sm.getOverlap(), 3.5 * _RadPerDeg)
        self.assertEquals(sm.getProjection(), "STG")
        
        for overlap in (0.0, 0.01, 0.1):
            sm = skymap.SkyMap(overlap = overlap)
            self.assertEquals(sm.getOverlap(), overlap)
            for tileId in range(sm.getNumSkyTiles()):
                tileInfo = sm.getSkyTileInfo(tileId)
                self.assertEquals(tileInfo.getOverlap(), overlap)
        
        for pixelScale in (1e-9, 1e-8, 1e-7):
            sm = skymap.SkyMap(pixelScale = pixelScale)
            self.assertEquals(sm.getPixelScale(), pixelScale)
        
        for projection in ("STG", "TAN", "MOL"):
            sm = skymap.SkyMap(projection = projection)
            self.assertEquals(sm.getProjection(), projection)
    
    def testTileSeparation(self):
        """Confirm that each sky tile has the proper distance to other tiles
        """
        sm = skymap.SkyMap()
        numSkyTiles = sm.getNumSkyTiles()
        tileInfoList = []
        for tileId in range(numSkyTiles):
            tileInfo = sm.getSkyTileInfo(tileId)
            self.assertEquals(tileInfo.getId(), tileId)
            tileInfoList.append(tileInfo)
        
        for tileInfo in tileInfoList:
            ctrCoord = tileInfo.getCtrCoord()
            distList = []
            for tileInfo1 in tileInfoList:
                otherCtrCoord = tileInfo1.getCtrCoord()
                distList.append(ctrCoord.angularSeparation(otherCtrCoord, afwCoord.DEGREES))
            distList.sort()
            self.assertAlmostEquals(distList[0], 0.0)
            for dist in distList[1:6]:
                self.assertAlmostEquals(dist, _NeighborAngularSeparation)
            self.assertAlmostEquals(distList[11], 180.0)
    
    def testGetSkyTileId(self):
        """Test the getSkyTileId method
        """
        sm = skymap.SkyMap()
        numSkyTiles = sm.getNumSkyTiles()
        tileInfoList = []
        for tileId in range(numSkyTiles):
            tileInfo = sm.getSkyTileInfo(tileId)
            tileInfoList.append(tileInfo)
        
        for tileInfo0 in tileInfoList:
            tileId0 = tileInfo0.getId()
            ctrCoord0 = tileInfo0.getCtrCoord()
            vector0 = numpy.array(ctrCoord0.getVector())
            
            # make a list of all 5 nearest neighbors
            nbrTileList = []
            for otherTileInfo in tileInfoList:
                otherCtrCoord = otherTileInfo.getCtrCoord()
                dist = ctrCoord0.angularSeparation(otherCtrCoord, afwCoord.DEGREES)
                if abs(dist - _NeighborAngularSeparation) < 0.1:
                    nbrTileList.append(otherTileInfo)
            self.assertEquals(len(nbrTileList), 5)
            
            for tileInfo1 in nbrTileList:
                tileId1 = tileInfo1.getId()
                ctrCoord1 = tileInfo1.getCtrCoord()
                vector1 = numpy.array(ctrCoord1.getVector())
                for tileInfo2 in nbrTileList[tileInfo1.getId():]:
                    dist = ctrCoord1.angularSeparation(tileInfo2.getCtrCoord(), afwCoord.DEGREES)
                    if abs(dist - _NeighborAngularSeparation) > 0.1:
                        continue
                    tileId2 = tileInfo2.getId()
                    ctrCoord2 = tileInfo2.getCtrCoord()
                    vector2 = numpy.array(ctrCoord2.getVector())
                
                    # sky tiles 0, 1 and 2 form a triangle of nearest neighbors
                    # explore the boundary between tile0 and tile1
                    # and also the boundary between tile0 and tile2
                    for deltaFrac in (-0.001, 0.001):
                        isNearest0 = deltaFrac > 0.0
                        
                        for exploreBoundary1 in (True, False):
                            # if exploreBoundary1, explore boundary between tile 0 and tile 1,
                            # else explore the boundary between tile 0 and tile 2
                        
                            if isNearest0:
                                expectedTileId = tileId0
                            elif exploreBoundary1:
                                expectedTileId = tileId1
                            else:
                                expectedTileId = tileId2
                            
                            for farFrac in (0.0, 0.05, 0.3, (1.0/3.0) - 0.01):
                                # farFrac is the fraction of the tile center vector point whose boundary
                                # is not being explored; it must be less than 1/3;
                                # remFrac is the remaining fraction, which is divided between tile 0
                                # and the tile whose boundary is being explored
                                remFrac = 1.0 - farFrac
                                frac0 = (remFrac / 2.0) + deltaFrac
                                boundaryFrac = (remFrac / 2.0) - deltaFrac
                                
                                if exploreBoundary1:
                                    frac2 = farFrac
                                    frac1 = boundaryFrac
                                else:
                                    frac1 = farFrac
                                    frac2 = boundaryFrac
    
                                testVector = (vector0 * frac0) + (vector1 * frac1) + (vector2 * frac2)
                                vecLen = math.sqrt(numpy.sum(testVector**2))
                                testVector /= vecLen
                                lsstVec = afwGeom.Point3D(testVector)
                                testCoord = afwCoord.IcrsCoord(lsstVec)
                                nearestTileId = sm.getSkyTileId(testCoord)
    
                                if expectedTileId != nearestTileId:
                                    nearestTileInfo = sm.getSkyTileInfo(nearestTileId)
                                    nearestCtrCoord = nearestTileInfo.getCtrCoord()
                                    nearestVector = nearestCtrCoord.getVector()
    
                                    print "tileId0=%s; tileId1=%s; tileId2=%s; nearestTileId=%s" % \
                                        (tileId0, tileId1, tileId2, nearestTileId)
                                    print "vector0=%s; vector1=%s; vector2=%s; nearestVector=%s" % \
                                         (vector0, vector1, vector2, nearestVector)
                                    print "frac0=%s; frac1=%s; frac2=%s" % (frac0, frac1, frac2)
                                    print "testVector=", testVector
    
                                    print "dist0=%s; dist1=%s; dist2=%s; nearDist=%s" % (
                                        testCoord.angularSeparation(ctrCoord0, afwCoord.DEGREES),
                                        testCoord.angularSeparation(ctrCoord1, afwCoord.DEGREES),
                                        testCoord.angularSeparation(ctrCoord2, afwCoord.DEGREES),
                                        testCoord.angularSeparation(nearestCtrCoord, afwCoord.DEGREES),
                                    )
                                    self.fail("Expected nearest tileId=%s; got tileId=%s" % \
                                        (expectedTileId, nearestTileId))
                    


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(SkyMapTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
