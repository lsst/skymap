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
from lsst.skymap import EquatSkyMap, skyMapRegistry


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

    def testDefaults(self):
        """Test important default values
        """
        config = EquatSkyMap.ConfigClass()
        skyMap = EquatSkyMap(config)
        self.assertEqual(skyMap.config.projection, "CEA")

    def testSymmetry(self):
        """Verify that the projection is symmetrical about the equator
        """
        for minDec in (-5.0, -1.0, 0.5):
            maxDec = minDec + 2.0
            config = EquatSkyMap.ConfigClass()
            config.decRange = minDec, maxDec
            skyMap = EquatSkyMap(config)
            for tractInfo in skyMap[0:1]:
                numPatches = tractInfo.getNumPatches()
                midXIndex = numPatches[0] / 2
                minPixelPosList = []
                maxPixelPosList = []
                maxYInd = numPatches[1] - 1
                for xInd in (0, midXIndex, numPatches[0] - 1):
                    minDecPatchInfo = tractInfo.getPatchInfo((xInd, 0))
                    minDecPosBox = afwGeom.Box2D(minDecPatchInfo.getOuterBBox())
                    minPixelPosList += [
                        minDecPosBox.getMin(),
                        afwGeom.Point2D(minDecPosBox.getMaxX(), minDecPosBox.getMinY()),
                    ]

                    maxDecPatchInfo = tractInfo.getPatchInfo((xInd, maxYInd))
                    maxDecPosBox = afwGeom.Box2D(maxDecPatchInfo.getOuterBBox())
                    maxPixelPosList += [
                        maxDecPosBox.getMax(),
                        afwGeom.Point2D(maxDecPosBox.getMinX(), maxDecPosBox.getMaxY()),
                    ]
                wcs = tractInfo.getWcs()
                minDecList = [wcs.pixelToSky(pos).getPosition(afwGeom.degrees)[1] for pos in minPixelPosList]
                maxDecList = [wcs.pixelToSky(pos).getPosition(afwGeom.degrees)[1] for pos in maxPixelPosList]
                self.assertTrue(numpy.allclose(minDecList, minDecList[0]))
                self.assertTrue(numpy.allclose(maxDecList, maxDecList[0]))
                self.assertTrue(minDecList[0] <= minDec)
                self.assertTrue(maxDecList[0] >= maxDec)

    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        for numTracts in (1, 2, 4, 25):
            config = EquatSkyMap.ConfigClass()
            config.numTracts = numTracts
            skyMap = EquatSkyMap(config)
            self.assertEqual(len(skyMap), numTracts)

        for tractOverlap in (0.0, 0.01, 0.1): # degrees
            config = EquatSkyMap.ConfigClass()
            config.tractOverlap = tractOverlap
            skyMap = EquatSkyMap(config)
            for tractInfo in skyMap:
                self.assertAlmostEqual(tractInfo.getTractOverlap().asDegrees(), tractOverlap)
            self.assertEqual(len(skyMap), skyMap.config.numTracts)

        for patchBorder in (0, 101):
            config = EquatSkyMap.ConfigClass()
            config.patchBorder = patchBorder
            skyMap = EquatSkyMap(config)
            for tractInfo in skyMap:
                self.assertEqual(tractInfo.getPatchBorder(), patchBorder)
            self.assertEqual(len(skyMap), skyMap.config.numTracts)

        skyMapClass = skyMapRegistry["equat"]
        for xInnerDim in (1005, 5062):
            for yInnerDim in (2032, 5431):
                config = skyMapClass.ConfigClass()
                config.patchInnerDimensions = (xInnerDim, yInnerDim)
                skyMap = EquatSkyMap(config)
                for tractInfo in skyMap:
                    self.assertEqual(tuple(tractInfo.getPatchInnerDimensions()), (xInnerDim, yInnerDim))
                self.assertEqual(len(skyMap), skyMap.config.numTracts)

    def testPickle(self):
        """Test that pickling and unpickling restores the original exactly
        """
        skyMap = EquatSkyMap()
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
            "numTracts",
            "decRange",
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
        for numTracts in (2, 4, 25):
            for minDec in (-45, -2.5, 32):
                for deltaDec in (2, 17):
                    maxDec = minDec + deltaDec
                    config = EquatSkyMap.ConfigClass()
                    config.numTracts = numTracts
                    config.decRange = (minDec, maxDec)
                    skyMap = EquatSkyMap(config)
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
            config = EquatSkyMap.ConfigClass()
            config.numTracts = numTracts
            skyMap = EquatSkyMap(config)
            decRange = skyMap.config.decRange
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

                        for testDecDeg in decList:
                            testDec = afwGeom.Angle(testDecDeg, afwGeom.degrees)
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
                    self.assertRaises(LookupError, tractInfo.findPatch, testCoord)


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
