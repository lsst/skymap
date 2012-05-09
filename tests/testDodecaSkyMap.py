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
"""Test DodecaSkyMap class
"""
import itertools
import os
import sys
import math
import pickle
import unittest

import numpy

import lsst.utils.tests as utilsTests
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from lsst.skymap import DodecaSkyMap, skyMapRegistry

# dodecahedron properties
_NumTracts = 12
_Phi = (1.0 + math.sqrt(5.0)) / 2.0
_DihedralAngle = afwGeom.Angle(2.0 * math.atan(_Phi), afwGeom.radians)
_NeighborAngularSeparation = 180.0 - _DihedralAngle.asDegrees()

class DodecaSkyMapTestCase(unittest.TestCase):
    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        for tractOverlap in (0.0, 0.01, 0.1): # degrees
            config = DodecaSkyMap.ConfigClass()
            config.tractOverlap = tractOverlap
            skyMap = DodecaSkyMap(config)
            for tractInfo in skyMap:
                self.assertAlmostEqual(tractInfo.getTractOverlap().asDegrees(), tractOverlap)
            self.assertEqual(len(skyMap), _NumTracts)

        for patchBorder in (0, 101):
            config = DodecaSkyMap.ConfigClass()
            config.patchBorder = patchBorder
            skyMap = DodecaSkyMap(config)
            for tractInfo in skyMap:
                self.assertEqual(tractInfo.getPatchBorder(), patchBorder)
            self.assertEqual(len(skyMap), _NumTracts)
 
        skyMapClass = skyMapRegistry["dodeca"]
        for xInnerDim in (1005, 5062):
            for yInnerDim in (2032, 5431):
                config = skyMapClass.ConfigClass()
                config.patchInnerDimensions = (xInnerDim, yInnerDim)
                skyMap = DodecaSkyMap(config)
                for tractInfo in skyMap:
                    self.assertEqual(tuple(tractInfo.getPatchInnerDimensions()), (xInnerDim, yInnerDim))
                self.assertEqual(len(skyMap), _NumTracts)
       
    def testPickle(self):
        """Test that pickling and unpickling restores the original exactly
        """
        skyMap = DodecaSkyMap()
        pickleStr = pickle.dumps(skyMap)
        unpickledSkyMap = pickle.loads(pickleStr)
        self.assertEqual(skyMap.getVersion(), unpickledSkyMap.getVersion())
        self.assertEqual(len(skyMap), len(unpickledSkyMap))
        for configName in (
            "patchInnerDimensions",
            "patchBorder",
            "projection",
            "pixelScale",
            "tractOverlap",
            "withTractsOnPoles",
        ):
            self.assertEqual(getattr(skyMap.config, configName), getattr(unpickledSkyMap.config, configName))
        for tractInfo, unpickledTractInfo in itertools.izip(skyMap, unpickledSkyMap):
            for getterName in (
                "getBBox",
                "getCtrCoord",
                "getId",
                "getNumPatches",
                "getPatchBorder",
                "getPatchInnerDimensions",
                "getTractOverlap",
                "getVertexList",
                "getWcs",
            ):
                self.assertEqual(getattr(tractInfo, getterName)(), getattr(unpickledTractInfo, getterName)())
            
            # test WCS at a few locations
            wcs = tractInfo.getWcs()
            unpickledWcs = unpickledTractInfo.getWcs()
            for x in (-1000.0, 0.0, 1000.0):
                for y in (-532.5, 0.5, 532.5):
                    pixelPos = afwGeom.Point2D(x, y)
                    skyPos = wcs.pixelToSky(pixelPos)
                    unpickledSkyPos = unpickledWcs.pixelToSky(pixelPos)
                    self.assertEqual(skyPos, unpickledSkyPos)
            
            # compare a few patches
            numPatches = tractInfo.getNumPatches()
            patchBorder = skyMap.config.patchBorder
            for xInd in (0, 1, numPatches[0]/2, numPatches[0]-2, numPatches[0]-1):
                for yInd in (0, 1, numPatches[1]/2, numPatches[1]-2, numPatches[1]-1):
                    patchInfo = tractInfo.getPatchInfo((xInd, yInd))
                    unpickledPatchInfo = unpickledTractInfo.getPatchInfo((xInd, yInd))
                    self.assertEqual(patchInfo, unpickledPatchInfo)
                    
                    # check inner and outer bbox (nothing to do with pickle,
                    # but a convenient place for the test)
                    innerBBox = patchInfo.getInnerBBox()
                    outerBBox = patchInfo.getOuterBBox()
                    
                    if xInd == 0:
                        self.assertEqual(innerBBox.getMinX(), outerBBox.getMinX())
                    else:
                        self.assertEqual(innerBBox.getMinX() - patchBorder, outerBBox.getMinX())
                    if yInd == 0:
                        self.assertEqual(innerBBox.getMinY(), outerBBox.getMinY())
                    else:
                        self.assertEqual(innerBBox.getMinY() - patchBorder, outerBBox.getMinY())
                        
                    if xInd == numPatches[0] - 1:
                        self.assertEqual(innerBBox.getMaxX(), outerBBox.getMaxX())
                    else:
                        self.assertEqual(innerBBox.getMaxX() + patchBorder, outerBBox.getMaxX())
                    if yInd == numPatches[1] - 1:
                        self.assertEqual(innerBBox.getMaxY(), outerBBox.getMaxY())
                    else:
                        self.assertEqual(innerBBox.getMaxY() + patchBorder, outerBBox.getMaxY())
                    
    
    def testTractSeparation(self):
        """Confirm that each sky tract has the proper distance to other tracts
        """
        skyMap = DodecaSkyMap()
        for tractId, tractInfo in enumerate(skyMap):
            self.assertEqual(tractInfo.getId(), tractId)
        
            ctrCoord = tractInfo.getCtrCoord()
            distList = []
            for tractInfo1 in skyMap:
                otherCtrCoord = tractInfo1.getCtrCoord()
                distList.append(ctrCoord.angularSeparation(otherCtrCoord).asDegrees())
            distList.sort()
            self.assertAlmostEquals(distList[0], 0.0)
            for dist in distList[1:6]:
                self.assertAlmostEquals(dist, _NeighborAngularSeparation)
            self.assertAlmostEquals(distList[11], 180.0)
    
    def testFindTract(self):
        """Test the findTract method
        """
        skyMap = DodecaSkyMap()
        for tractInfo0 in skyMap:
            tractId0 = tractInfo0.getId()
            ctrCoord0 = tractInfo0.getCtrCoord()
            vector0 = numpy.array(ctrCoord0.getVector())
            
            # make a list of all 5 nearest neighbors
            nbrTractList = []
            for otherTractInfo in skyMap:
                otherCtrCoord = otherTractInfo.getCtrCoord()
                dist = ctrCoord0.angularSeparation(otherCtrCoord).asDegrees()
                if abs(dist - _NeighborAngularSeparation) < 0.1:
                    nbrTractList.append(otherTractInfo)
            self.assertEqual(len(nbrTractList), 5)
            
            for tractInfo1 in nbrTractList:
                tractId1 = tractInfo1.getId()
                ctrCoord1 = tractInfo1.getCtrCoord()
                vector1 = numpy.array(ctrCoord1.getVector())
                for tractInfo2 in nbrTractList[tractInfo1.getId():]:
                    dist = ctrCoord1.angularSeparation(tractInfo2.getCtrCoord()).asDegrees()
                    if abs(dist - _NeighborAngularSeparation) > 0.1:
                        continue
                    tractId2 = tractInfo2.getId()
                    ctrCoord2 = tractInfo2.getCtrCoord()
                    vector2 = numpy.array(ctrCoord2.getVector())
                
                    # sky tracts 0, 1 and 2 form a triangle of nearest neighbors
                    # explore the boundary between tract 0 and tract 1
                    # and also the boundary between tract 0 and tract 2
                    for deltaFrac in (-0.001, 0.001):
                        isNearest0 = deltaFrac > 0.0
                        
                        for exploreBoundary1 in (True, False):
                            # if exploreBoundary1, explore boundary between tract 0 and tract 1,
                            # else explore the boundary between tract 0 and tract 2
                        
                            if isNearest0:
                                expectedTractId = tractId0
                            elif exploreBoundary1:
                                expectedTractId = tractId1
                            else:
                                expectedTractId = tractId2
                            
                            for farFrac in (0.0, 0.05, 0.3, (1.0/3.0) - 0.01):
                                # farFrac is the fraction of the tract center vector point whose boundary
                                # is not being explored; it must be less than 1/3;
                                # remFrac is the remaining fraction, which is divided between tract 0
                                # and the tract whose boundary is being explored
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
                                nearestTractInfo = skyMap.findTract(testCoord)
                                nearestTractId = nearestTractInfo.getId()
    
                                if expectedTractId != nearestTractId:
                                    nearestCtrCoord = nearestTractInfo.getCtrCoord()
                                    nearestVector = nearestCtrCoord.getVector()
    
                                    print "tractId0=%s; tractId1=%s; tractId2=%s; nearestTractId=%s" % \
                                        (tractId0, tractId1, tractId2, nearestTractId)
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
                                    self.fail("Expected nearest tractId=%s; got tractId=%s" % \
                                        (expectedTractId, nearestTractId))
                                
                                patchInfo = nearestTractInfo.findPatch(testCoord)
                                pixelInd = afwGeom.Point2I(
                                    nearestTractInfo.getWcs().skyToPixel(testCoord.toIcrs()))
                                self.assertTrue(patchInfo.getInnerBBox().contains(pixelInd))


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(DodecaSkyMapTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
