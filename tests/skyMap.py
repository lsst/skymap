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
_NumPatches = 12
_Phi = (1.0 + math.sqrt(5.0)) / 2.0
_DihedralAngle = 2.0 * math.atan(_Phi) / _RadPerDeg
_NeighborAngularSeparation = 180.0 - _DihedralAngle

class SkyMapTestCase(unittest.TestCase):
    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        sm = skymap.SkyMap()
        self.assertEqual(sm.getNumSkyPatches(), _NumPatches)
        self.assertEqual(sm.getOverlap(), 3.5 * _RadPerDeg)
        self.assertEqual(sm.getProjection(), "STG")
        
        for overlap in (0.0, 0.01, 0.1): # degrees
            sm = skymap.SkyMap(overlap = afwGeom.Angle(overlap, afwGeom.degrees))
            self.assertEqual(sm.getOverlap().asDegrees(), overlap)
            for patchId in range(sm.getNumSkyPatches()):
                patchInfo = sm.getSkyPatchInfo(patchId)
                self.assertAlmostEqual(patchInfo.getOverlap().asDegrees(), overlap)
        
        for pixelScale in (0.01, 0.1, 1.0): # arcseconds/pixel
            sm = skymap.SkyMap(pixelScale = afwGeom.Angle(pixelScale, afwGeom.arcseconds))
            self.assertAlmostEqual(sm.getPixelScale().asArcseconds(), pixelScale)
        
        for projection in ("STG", "TAN", "MOL"):
            sm = skymap.SkyMap(projection = projection)
            self.assertEqual(sm.getProjection(), projection)
    
    def testPatchSeparation(self):
        """Confirm that each sky patch has the proper distance to other patches
        """
        sm = skymap.SkyMap()
        numSkyPatches = sm.getNumSkyPatches()
        patchInfoList = []
        for patchId in range(numSkyPatches):
            patchInfo = sm.getSkyPatchInfo(patchId)
            self.assertEqual(patchInfo.getId(), patchId)
            patchInfoList.append(patchInfo)
        
        for patchInfo in patchInfoList:
            ctrCoord = patchInfo.getCtrCoord()
            distList = []
            for patchInfo1 in patchInfoList:
                otherCtrCoord = patchInfo1.getCtrCoord()
                distList.append(ctrCoord.angularSeparation(otherCtrCoord).asDegrees())
            distList.sort()
            self.assertAlmostEquals(distList[0], 0.0)
            for dist in distList[1:6]:
                self.assertAlmostEquals(dist, _NeighborAngularSeparation)
            self.assertAlmostEquals(distList[11], 180.0)
    
    def testGetSkyPatchId(self):
        """Test the getSkyPatchId method
        """
        sm = skymap.SkyMap()
        numSkyPatches = sm.getNumSkyPatches()
        patchInfoList = []
        for patchId in range(numSkyPatches):
            patchInfo = sm.getSkyPatchInfo(patchId)
            patchInfoList.append(patchInfo)
        
        for patchInfo0 in patchInfoList:
            patchId0 = patchInfo0.getId()
            ctrCoord0 = patchInfo0.getCtrCoord()
            vector0 = numpy.array(ctrCoord0.getVector())
            
            # make a list of all 5 nearest neighbors
            nbrPatchList = []
            for otherPatchInfo in patchInfoList:
                otherCtrCoord = otherPatchInfo.getCtrCoord()
                dist = ctrCoord0.angularSeparation(otherCtrCoord).asDegrees()
                if abs(dist - _NeighborAngularSeparation) < 0.1:
                    nbrPatchList.append(otherPatchInfo)
            self.assertEqual(len(nbrPatchList), 5)
            
            for patchInfo1 in nbrPatchList:
                patchId1 = patchInfo1.getId()
                ctrCoord1 = patchInfo1.getCtrCoord()
                vector1 = numpy.array(ctrCoord1.getVector())
                for patchInfo2 in nbrPatchList[patchInfo1.getId():]:
                    dist = ctrCoord1.angularSeparation(patchInfo2.getCtrCoord()).asDegrees()
                    if abs(dist - _NeighborAngularSeparation) > 0.1:
                        continue
                    patchId2 = patchInfo2.getId()
                    ctrCoord2 = patchInfo2.getCtrCoord()
                    vector2 = numpy.array(ctrCoord2.getVector())
                
                    # sky patches 0, 1 and 2 form a triangle of nearest neighbors
                    # explore the boundary between patch 0 and patch 1
                    # and also the boundary between patch 0 and patch 2
                    for deltaFrac in (-0.001, 0.001):
                        isNearest0 = deltaFrac > 0.0
                        
                        for exploreBoundary1 in (True, False):
                            # if exploreBoundary1, explore boundary between patch 0 and patch 1,
                            # else explore the boundary between patch 0 and patch 2
                        
                            if isNearest0:
                                expectedPatchId = patchId0
                            elif exploreBoundary1:
                                expectedPatchId = patchId1
                            else:
                                expectedPatchId = patchId2
                            
                            for farFrac in (0.0, 0.05, 0.3, (1.0/3.0) - 0.01):
                                # farFrac is the fraction of the patch center vector point whose boundary
                                # is not being explored; it must be less than 1/3;
                                # remFrac is the remaining fraction, which is divided between patch 0
                                # and the patch whose boundary is being explored
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
                                nearestPatchId = sm.getSkyPatchId(testCoord)
    
                                if expectedPatchId != nearestPatchId:
                                    nearestPatchInfo = sm.getSkyPatchInfo(nearestPatchId)
                                    nearestCtrCoord = nearestPatchInfo.getCtrCoord()
                                    nearestVector = nearestCtrCoord.getVector()
    
                                    print "patchId0=%s; patchId1=%s; patchId2=%s; nearestPatchId=%s" % \
                                        (patchId0, patchId1, patchId2, nearestPatchId)
                                    print "vector0=%s; vector1=%s; vector2=%s; nearestVector=%s" % \
                                         (vector0, vector1, vector2, nearestVector)
                                    print "frac0=%s; frac1=%s; frac2=%s" % (frac0, frac1, frac2)
                                    print "testVector=", testVector
    
                                    print "dist0=%s; dist1=%s; dist2=%s; nearDist=%s" % (
                                        testCoord.angularSeparation(ctrCoord0).asDegrees(),
                                        testCoord.angularSeparation(ctrCoord1).asDegrees(),
                                        testCoord.angularSeparation(ctrCoord2).asDegrees(),
                                        testCoord.angularSeparation(nearestCtrCoord).asDegrees(),
                                    )
                                    self.fail("Expected nearest patchId=%s; got patchId=%s" % \
                                        (expectedPatchId, nearestPatchId))
                    


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
