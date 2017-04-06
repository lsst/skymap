#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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
from __future__ import absolute_import, division, print_function
from builtins import zip
import pickle

import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.utils.tests

from lsst.skymap import skyMapRegistry


class SkyMapTestCase(lsst.utils.tests.TestCase):
    """An abstract base class for testing a SkyMap.

    To use, subclass and set the following class variables:
    _SkyMapClass: the particular SkyMap subclass to test
    _SkyMapConfig: a SkyMapConfig instance to use in the test, or None to use _SkyMapClass.ConfigClass()
    _SkyMapName: the name of the particular SkyMap class in the registry
    _NumTracts: the number of tracts to expect (for the default configuration)
    _NeighborAngularSeparation: Expected angular separation between tracts
    _numNeighbors: Number of neighbors that should be within the expected angular separation
    """
    _SkyMapConfig = None

    def getSkyMap(self, config=None):
        """Provide an instance of the skymap"""
        if config is None:
            config = self.getConfig()
        return self._SkyMapClass(config=config)

    def getConfig(self):
        """Provide an instance of the configuration class"""
        if self._SkyMapConfig is None:
            return self._SkyMapClass.ConfigClass()
        # Want to return a copy of self._SkyMapConfig, so it can be modified.
        # However, there is no Config.copy() method, so this is more complicated than desirable.
        return pickle.loads(pickle.dumps(self._SkyMapConfig))

    def testRegistry(self):
        """Confirm that the skymap can be retrieved from the registry"""
        self.assertEqual(skyMapRegistry[self._SkyMapName], self._SkyMapClass)

    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        for tractOverlap in (0.0, 0.01, 0.1):  # degrees
            config = self.getConfig()
            config.tractOverlap = tractOverlap
            skyMap = self.getSkyMap(config)
            for tractInfo in skyMap:
                self.assertAlmostEqual(tractInfo.getTractOverlap().asDegrees(), tractOverlap)
            self.assertEqual(len(skyMap), self._NumTracts)

        for patchBorder in (0, 101):
            config = self.getConfig()
            config.patchBorder = patchBorder
            skyMap = self.getSkyMap(config)
            for tractInfo in skyMap:
                self.assertEqual(tractInfo.getPatchBorder(), patchBorder)
            self.assertEqual(len(skyMap), self._NumTracts)

        for xInnerDim in (1005, 5062):
            for yInnerDim in (2032, 5431):
                config = self.getConfig()
                config.patchInnerDimensions = (xInnerDim, yInnerDim)
                skyMap = self.getSkyMap(config)
                for tractInfo in skyMap:
                    self.assertEqual(tuple(tractInfo.getPatchInnerDimensions()), (xInnerDim, yInnerDim))
                self.assertEqual(len(skyMap), self._NumTracts)

    def assertUnpickledTractInfo(self, unpickled, original, patchBorder):
        """Assert that an unpickled TractInfo is functionally identical to the original

        @param unpickled      The unpickled TractInfo
        @param original       The original TractInfo
        @param patchBorder    Border around each patch, from SkyMap.config.patchBorder
        """
        for getterName in ("getBBox",
                           "getCtrCoord",
                           "getId",
                           "getNumPatches",
                           "getPatchBorder",
                           "getPatchInnerDimensions",
                           "getTractOverlap",
                           "getVertexList",
                           "getWcs",
                           ):
            self.assertEqual(getattr(original, getterName)(), getattr(unpickled, getterName)())

        # test WCS at a few locations
        wcs = original.getWcs()
        unpickledWcs = unpickled.getWcs()
        for x in (-1000.0, 0.0, 1000.0):
            for y in (-532.5, 0.5, 532.5):
                pixelPos = afwGeom.Point2D(x, y)
                skyPos = wcs.pixelToSky(pixelPos)
                unpickledSkyPos = unpickledWcs.pixelToSky(pixelPos)
                self.assertEqual(skyPos, unpickledSkyPos)

        # compare a few patches
        numPatches = original.getNumPatches()
        for xInd in (0, 1, numPatches[0]//2, numPatches[0]-2, numPatches[0]-1):
            for yInd in (0, 1, numPatches[1]//2, numPatches[1]-2, numPatches[1]-1):
                patchInfo = original.getPatchInfo((xInd, yInd))
                unpickledPatchInfo = unpickled.getPatchInfo((xInd, yInd))
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

    def testPickle(self):
        """Test that pickling and unpickling restores the original exactly
        """
        skyMap = self.getSkyMap()
        pickleStr = pickle.dumps(skyMap)
        unpickledSkyMap = pickle.loads(pickleStr)
        self.assertEqual(len(skyMap), len(unpickledSkyMap))
        self.assertEqual(unpickledSkyMap.config, skyMap.config)
        for tractInfo, unpickledTractInfo in zip(skyMap, unpickledSkyMap):
            self.assertUnpickledTractInfo(unpickledTractInfo, tractInfo, skyMap.config.patchBorder)

    def testTractSeparation(self):
        """Confirm that each sky tract has the proper distance to other tracts
        """
        skyMap = self.getSkyMap()
        for tractId, tractInfo in enumerate(skyMap):
            self.assertEqual(tractInfo.getId(), tractId)

            ctrCoord = tractInfo.getCtrCoord()
            distList = []
            for tractInfo1 in skyMap:
                otherCtrCoord = tractInfo1.getCtrCoord()
                distList.append(ctrCoord.angularSeparation(otherCtrCoord))
            distList.sort()
            self.assertEqual(distList[0], 0.0)
            for dist in distList[1:self._numNeighbors]:
                self.assertAnglesAlmostEqual(dist, self._NeighborAngularSeparation)

    def testFindPatchList(self):
        """Test findTract.findPatchList
        """
        skyMap = self.getSkyMap()
        for tractId in (0, 5):
            tractInfo = skyMap[tractId]
            wcs = tractInfo.getWcs()
            numPatches = tractInfo.getNumPatches()
            border = tractInfo.getPatchBorder()
            for patchInd in ((0, 0),
                             (0, 1),
                             (5, 0),
                             (5, 6),
                             (numPatches[0] - 2, numPatches[1] - 1),
                             (numPatches[0] - 1, numPatches[1] - 2),
                             (numPatches[0] - 1, numPatches[1] - 1),
                             ):
                patchInfo = tractInfo.getPatchInfo(patchInd)
                patchIndex = patchInfo.getIndex()
                bbox = patchInfo.getInnerBBox()
                bbox.grow(-(border+1))
                coordList = getCornerCoords(wcs=wcs, bbox=bbox)
                patchInfoList = tractInfo.findPatchList(coordList)
                self.assertEqual(len(patchInfoList), 1)
                self.assertEqual(patchInfoList[0].getIndex(), patchIndex)

                # grow to include neighbors and test again
                bbox.grow(2)
                predFoundIndexSet = set()
                for dx in (-1, 0, 1):
                    nbrX = patchIndex[0] + dx
                    if not 0 <= nbrX < numPatches[0]:
                        continue
                    for dy in (-1, 0, 1):
                        nbrY = patchIndex[1] + dy
                        if not 0 <= nbrY < numPatches[1]:
                            continue
                        nbrInd = (nbrX, nbrY)
                        predFoundIndexSet.add(nbrInd)
                coordList = getCornerCoords(wcs=wcs, bbox=bbox)
                patchInfoList = tractInfo.findPatchList(coordList)
                self.assertEqual(len(patchInfoList), len(predFoundIndexSet))
                foundIndexSet = set(patchInfo.getIndex() for patchInfo in patchInfoList)
                self.assertEqual(foundIndexSet, predFoundIndexSet)

    def testFindTractPatchList(self):
        """Test findTractPatchList

        Note: this test uses single points for speed and to avoid really large regions.
        Note that findPatchList is being tested elsewhere.
        """
        skyMap = self.getSkyMap()
        for tractId in (1, 3, 7):
            tractInfo = skyMap[tractId]
            self.assertTractPatchListOk(
                skyMap=skyMap,
                coordList=[tractInfo.getCtrCoord()],
                knownTractId=tractId,
            )
            self.assertClosestTractPatchList(skyMap, [tractInfo.getCtrCoord()], tractId)

            vertices = tractInfo.getVertexList()
            if len(vertices) > 0:
                self.assertTractPatchListOk(
                    skyMap=skyMap,
                    coordList=[tractInfo.getVertexList()[0]],
                    knownTractId=tractId,
                )

            if len(vertices) > 2:
                self.assertTractPatchListOk(
                    skyMap=skyMap,
                    coordList=[tractInfo.getVertexList()[2]],
                    knownTractId=tractId,
                )

    def testTractContains(self):
        """Test that TractInfo.contains works"""
        skyMap = self.getSkyMap()
        for tract in skyMap:
            coord = tract.getCtrCoord()
            self.assertTrue(tract.contains(coord))
            opposite = afwCoord.IcrsCoord(coord.getLongitude() + 12*afwGeom.hours, -1*coord.getLatitude())
            self.assertFalse(tract.contains(opposite))

    def assertTractPatchListOk(self, skyMap, coordList, knownTractId):
        """Assert that findTractPatchList produces the correct results

        @param[in] skyMap: sky map to test
        @param[in] coordList: coordList of region to search for
        @param[in] knownTractId: this tractId must appear in the found list
        """
        tractPatchList = skyMap.findTractPatchList(coordList)
        tractPatchDict = dict((tp[0].getId(), tp[1]) for tp in tractPatchList)
        self.assertTrue(knownTractId in tractPatchDict)
        for tractInfo in skyMap:
            tractId = tractInfo.getId()
            patchList = tractInfo.findPatchList(coordList)
            if patchList:
                self.assertTrue(tractId in tractPatchDict)
                self.assertTrue(len(patchList) == len(tractPatchDict[tractId]))
            else:
                self.assertTrue(tractId not in tractPatchDict)

    def assertClosestTractPatchList(self, skyMap, coordList, knownTractId):
        if not hasattr(skyMap, "findClosestTractPatchList"):
            self.skipTest("This skymap doesn't implement findClosestTractPatchList")
        tractPatchList = skyMap.findClosestTractPatchList(coordList)
        self.assertEqual(len(coordList), len(tractPatchList))  # One tract+patchList per coordinate
        for coord, (tract, patchList) in zip(coordList, tractPatchList):
            self.assertEqual(tract.getId(), knownTractId)
            self.assertEqual(patchList, tract.findPatchList([coord]))


def getCornerCoords(wcs, bbox):
    """Return the coords of the four corners of a bounding box
    """
    bbox = afwGeom.Box2D(bbox)  # mak
    cornerPosList = (
        bbox.getMin(),
        afwGeom.Point2D(bbox.getMaxX(), bbox.getMinY()),
        bbox.getMax(),
        afwGeom.Point2D(bbox.getMinX(), bbox.getMaxY()),
    )
    return [wcs.pixelToSky(cp).toIcrs() for cp in cornerPosList]
