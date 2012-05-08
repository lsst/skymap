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
"""Test EquatSkyMap class
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
from lsst.skymap import EquatSkyMap

class EquatSkyMapTestCase(unittest.TestCase):
    def getNeighborTracts(self, skyMap, tractId):
        """Return previous and next tractInfo
        """
        if not 0 <= tractId < len(skyMap):
            raise RuntimeError("Invalid tractId=%s" % (tractId,))
        if tractId > 0:
            prevTract = skyMap[tractId - 1]
        else:
            prevTract = skyMap[-1]
        if tractId + 1 < len(skyMap):
            nextTract = skyMap[tractId + 1]
        else:
            nextTract = skyMap[0]
        return (prevTract, nextTract)

    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        for tractOverlap in (0.0, 0.01, 0.1): # degrees
            skyMap = EquatSkyMap(tractOverlap = afwGeom.Angle(tractOverlap, afwGeom.degrees))
            self.assertEqual(skyMap.getTractOverlap().asDegrees(), tractOverlap)
            for tractInfo in skyMap:
                self.assertAlmostEqual(tractInfo.getTractOverlap().asDegrees(), tractOverlap)
        
        for pixelScale in (0.01, 0.1, 1.0): # arcseconds/pixel
            skyMap = EquatSkyMap(pixelScale = afwGeom.Angle(pixelScale, afwGeom.arcseconds))
            self.assertAlmostEqual(skyMap.getPixelScale().asArcseconds(), pixelScale)
        
        for projection in ("CEA", "TAN", "MOL"):
            # use numTracts > 2 to avoid "TAN" failing in wcslib
            skyMap = EquatSkyMap(numTracts = 4, projection = projection)
            self.assertEqual(skyMap.getProjection(), projection)

        for numTracts in (1, 2, 4, 25):
            skyMap = EquatSkyMap(numTracts = numTracts)
            self.assertEqual(len(skyMap), numTracts)
        
        for minDec in (-45.0, -3.7, 2.5):
            for maxDec in (-2.5, 3.7, 69.0):
                if maxDec <= minDec:
                    continue
                decRange = (afwGeom.Angle(minDec, afwGeom.degrees), afwGeom.Angle(maxDec, afwGeom.degrees))
                skyMap = EquatSkyMap(decRange = decRange)
                self.assertEqual(skyMap.getDecRange(), decRange)

    def testPickle(self):
        """Test that pickling and unpickling restores the original exactly
        """
        skyMap = EquatSkyMap()
        pickleStr = pickle.dumps(skyMap)
        unpickledSkyMap = pickle.loads(pickleStr)
        self.assertEqual(len(skyMap), len(unpickledSkyMap))
        for getterName in (
            "getPatchInnerDimensions",
            "getPatchBorder",
            "getProjection",
            "getPixelScale",
            "getTractOverlap",
            "getVersion",
        ):
            self.assertEqual(getattr(skyMap, getterName)(), getattr(unpickledSkyMap, getterName)())
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
            patchBorder = tractInfo.getPatchBorder()
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
        for numTracts in (2, 4, 25):
            for minDec in (-45, -2.5, 32):
                for deltaDec in (2, 17):
                    maxDec = minDec + deltaDec
                    decRange = (afwGeom.Angle(minDec, afwGeom.degrees), afwGeom.Angle(maxDec, afwGeom.degrees))
                    
                    skyMap = EquatSkyMap(numTracts = numTracts, decRange = decRange)
                    predDeltaRa = 360.0 / numTracts
            
                    for tractId, tractInfo in enumerate(skyMap):
                        self.assertEqual(tractInfo.getId(), tractId)
                        
                        prevTract, nextTract = self.getNeighborTracts(skyMap, tractId)

                        deltaRa = tractInfo.getCtrCoord().getRa() - prevTract.getCtrCoord().getRa()
                        deltaRa = deltaRa.asDegrees()
                        raError = abs(deltaRa - predDeltaRa) % 180.0
                        self.assertAlmostEquals(raError, 0.0)
                        
                        deltaRa = nextTract.getCtrCoord().getRa() - tractInfo.getCtrCoord().getRa()
                        deltaRa = deltaRa.asDegrees()
                        raError = abs(deltaRa - predDeltaRa) % 180.0
                        self.assertAlmostEquals(raError, 0.0)
    
    def testFindTract(self):
        """Test the findTract method
        """
        for numTracts in (2, 4):
            skyMap = EquatSkyMap(numTracts = numTracts)
            decRange = skyMap.getDecRange()
            decList = (
                (decRange[0] * 0.999) + (decRange[1] * 0.901),
                (decRange[0] * 0.500) + (decRange[1] * 0.500),
                (decRange[0] * 0.091) + (decRange[1] * 0.999),
            )
            for tractInfo0 in skyMap:
                tractId0 = tractInfo0.getId()
                ctrCoord0 = tractInfo0.getCtrCoord()
                vector0 = numpy.array(ctrCoord0.getVector())
                
                for tractInfo1 in self.getNeighborTracts(skyMap, tractId0):
                
                    tractId1 = tractInfo1.getId()
                    ctrCoord1 = tractInfo1.getCtrCoord()
                    vector1 = numpy.array(ctrCoord1.getVector())
                
                    for deltaFrac in (-0.001, 0.001):
                        # this fuss is because Point3D does not support * float
                        v0 = [v * (0.5 + deltaFrac) for v in ctrCoord0.getVector()]
                        v1 = [v * (0.5 - deltaFrac) for v in ctrCoord1.getVector()]
                        testVec = afwGeom.Point3D(*(v0[i] + v1[i] for i in range(3)))
                        testRa = afwCoord.IcrsCoord(testVec).getRa()
    
                        if deltaFrac > 0.0:
                            expectedTractId = tractId0
                        else:
                            expectedTractId = tractId1
                        
                        for testDec in decList:
                            testCoord = afwCoord.IcrsCoord(testRa, testDec)
                        
                            nearestTractInfo = skyMap.findTract(testCoord)
                            nearestTractId = nearestTractInfo.getId()
    
                            self.assertEqual(nearestTractId, expectedTractId)
                            
                            patchInfo = nearestTractInfo.findPatch(testCoord)
                            pixelInd = afwGeom.Point2I(
                                nearestTractInfo.getWcs().skyToPixel(testCoord.toIcrs()))
                            self.assertTrue(patchInfo.getInnerBBox().contains(pixelInd))
                
                # find a point outside the tract and make sure it fails
                tractInfo = tractInfo0
                wcs = tractInfo.getWcs()
                bbox = afwGeom.Box2D(tractInfo.getBBox())
                outerPixPosList = [
                    bbox.getMin() - afwGeom.Extent2D(1, 1),
                    afwGeom.Point2D(bbox.getMaxX(), bbox.getMinY()) - afwGeom.Extent2D(1, 1),
                    bbox.getMax() + afwGeom.Extent2D(1, 1),
                    afwGeom.Point2D(bbox.getMinX(), bbox.getMaxY()) + afwGeom.Extent2D(1, 1),
                ]
                for outerPixPos in outerPixPosList:
                    testCoord = wcs.pixelToSky(outerPixPos)
                    self.assertRaises(RuntimeError, tractInfo.findPatch, testCoord)
            

def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(EquatSkyMapTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
