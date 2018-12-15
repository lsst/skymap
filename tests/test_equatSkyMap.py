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
import unittest

import numpy

import lsst.afw.geom as afwGeom
import lsst.utils.tests

from lsst.skymap import EquatSkyMap
from helper import skyMapTestCase


class EquatSkyMapTestCase(skyMapTestCase.SkyMapTestCase):
    def setUp(self):
        self.setAttributes(
            SkyMapClass=EquatSkyMap,
            name="equat",
            numTracts=4,
            neighborAngularSeparation=90*afwGeom.degrees,
            numNeighbors=2,
        )

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

    def test_gettingPatchInfo(self):
        skyMap = EquatSkyMap(EquatSkyMap.ConfigClass())
        tract = skyMap[0]
        # Look up patchInfo with a coordinate pair
        patchIndex = (3, 5)
        pairPatchInfo = tract.getPatchInfo(patchIndex)
        # Fetch the sequential index
        sequentialPatchIndex = tract.getSequentialPatchIndex(pairPatchInfo)
        # Turn the sequential index back into a pair
        returnedPatchIndex = tract.getPatchIndexPair(sequentialPatchIndex)
        # Test that the different indexes match
        self.assertEqual(patchIndex, returnedPatchIndex)
        # verify patch info can be retrieved with both indexes
        sequentialPatchInfo = tract.getPatchInfo(sequentialPatchIndex)
        self.assertEqual(pairPatchInfo, sequentialPatchInfo)

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
                midXIndex = numPatches[0]//2
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

    def testMoreBasicAttributes(self):
        """Confirm that constructor attributes are available.
        """
        defaultSkyMap = self.getSkyMap()
        for numTracts in (1, 2, 4, 25):
            config = EquatSkyMap.ConfigClass()
            config.numTracts = numTracts
            skyMap = EquatSkyMap(config)
            self.assertEqual(len(skyMap), numTracts)
            if numTracts == defaultSkyMap.config.numTracts:
                self.assertEqual(skyMap, defaultSkyMap)
            else:
                self.assertNotEqual(skyMap, defaultSkyMap)

        for decRange in ([-1.3, 0.5], [6.1, 6.8]):
            config = EquatSkyMap.ConfigClass()
            config.decRange = decRange
            skyMap = EquatSkyMap(config)
            self.assertNotEqual(skyMap, defaultSkyMap)

    def testFindTract(self):
        """Test the SkyMap.findTract method
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

                for tractInfo1 in self.getNeighborTracts(skyMap, tractId0):

                    tractId1 = tractInfo1.getId()
                    ctrCoord1 = tractInfo1.getCtrCoord()

                    for deltaFrac in (-0.001, 0.001):
                        v0 = ctrCoord0.getVector() * (0.5 + deltaFrac)
                        v1 = ctrCoord1.getVector() * (0.5 - deltaFrac)
                        testVec = v0 + v1
                        testRa = afwGeom.SpherePoint(testVec).getRa()

                        if deltaFrac > 0.0:
                            expectedTractId = tractId0
                        else:
                            expectedTractId = tractId1

                        for testDecDeg in decList:
                            testDec = afwGeom.Angle(testDecDeg, afwGeom.degrees)
                            testCoord = afwGeom.SpherePoint(testRa, testDec)

                            nearestTractInfo = skyMap.findTract(testCoord)
                            nearestTractId = nearestTractInfo.getId()

                            self.assertEqual(nearestTractId, expectedTractId)

                            patchInfo = nearestTractInfo.findPatch(testCoord)
                            pixelInd = afwGeom.Point2I(nearestTractInfo.getWcs().skyToPixel(testCoord))
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
