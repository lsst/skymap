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
import itertools
import pickle

import numpy as np

import lsst.geom as geom
import lsst.utils.tests

from lsst.skymap import skyMapRegistry


def checkDm14809(testcase, skymap):
    """Test that DM-14809 has been fixed

    The observed behaviour was:

        skyMap.findTract(skyMap[9712].getCtrCoord()).getId() != 9712

    and

        skyMap[1].getCtrCoord() == skyMap[11].getCtrCoord()

    In order to be thorough, we generalise these over the entire skymap.
    """
    # Check that the tract found for central coordinate of a tract is that tract
    expect = [tract.getId() for tract in skymap]
    got = [skymap.findTract(tract.getCtrCoord()).getId() for tract in skymap]
    testcase.assertListEqual(got, expect)

    # Check that the tract central coordinates are unique
    # Round to integer arcminutes so differences are relatively immune to small numerical inaccuracies
    centers = set([(int(coord.getRa().asArcminutes()), int(coord.getDec().asArcminutes())) for
                   coord in (tract.getCtrCoord() for tract in skymap)])
    testcase.assertEqual(len(centers), len(skymap))


class SkyMapTestCase(lsst.utils.tests.TestCase):
    """An abstract base class for testing a SkyMap.

    To use, subclass and call `setAttributes` from `setUp`
    """
    def setAttributes(self, *,
                      SkyMapClass,
                      name,
                      numTracts,
                      config=None,
                      neighborAngularSeparation=None,
                      numNeighbors=None):
        """Initialize the test (call from setUp in the subclass)

        Parameters
        ----------
        SkyMapClass : subclass of `lsst.skymap.BaseSkyMap`
            Class of sky map to test
        name : `str`
            Name of sky map in sky map registry
        numTracts : `int`
            Number of tracts in the default configuration
        config : subclass of `lsst.skymap.SkyMapConfig`, optional
            Default configuration used by `getSkyMap`;
            if None use SkyMapClass.ConfigClass()
        neighborAngularSeparation : `lsst.geom.Angle`, optional
            Expected angular separation between tracts;
            if None then angular separation is not tested unless your
            subclass of SkyMapTestCaseoverrides `testTractSeparation`.
        numNeighbors : `int` or `None`
            Number of neighbors that should be within
            ``neighborAngularSeparation``;
            Required if ``neighborAngularSeparation`` is not None;
            ignored otherwise.
        """
        self.SkyMapClass = SkyMapClass
        self.config = config
        self.name = name
        self.numTracts = numTracts
        self.neighborAngularSeparation = neighborAngularSeparation
        self.numNeighbors = numNeighbors
        np.random.seed(47)

    def getSkyMap(self, config=None):
        """Provide an instance of the skymap"""
        if config is None:
            config = self.getConfig()
        return self.SkyMapClass(config=config)

    def getConfig(self):
        """Provide an instance of the configuration class"""
        if self.config is None:
            return self.SkyMapClass.ConfigClass()
        # Want to return a copy of self.config, so it can be modified.
        # However, there is no Config.copy() method, so this is more complicated than desirable.
        return pickle.loads(pickle.dumps(self.config))

    def testRegistry(self):
        """Confirm that the skymap can be retrieved from the registry"""
        self.assertEqual(skyMapRegistry[self.name], self.SkyMapClass)

    def testBasicAttributes(self):
        """Confirm that constructor attributes are available
        """
        defaultSkyMap = self.getSkyMap()
        for tractOverlap in (0.0, 0.01, 0.1):  # degrees
            config = self.getConfig()
            config.tractOverlap = tractOverlap
            skyMap = self.getSkyMap(config)
            for tractInfo in skyMap:
                self.assertAlmostEqual(tractInfo.getTractOverlap().asDegrees(), tractOverlap)
            self.assertEqual(len(skyMap), self.numTracts)
            self.assertNotEqual(skyMap, defaultSkyMap)

        for patchBorder in (0, 101):
            config = self.getConfig()
            config.patchBorder = patchBorder
            skyMap = self.getSkyMap(config)
            for tractInfo in skyMap:
                self.assertEqual(tractInfo.getPatchBorder(), patchBorder)
            self.assertEqual(len(skyMap), self.numTracts)
            self.assertNotEqual(skyMap, defaultSkyMap)

        for xInnerDim in (1005, 5062):
            for yInnerDim in (2032, 5431):
                config = self.getConfig()
                config.patchInnerDimensions = (xInnerDim, yInnerDim)
                skyMap = self.getSkyMap(config)
                for tractInfo in skyMap:
                    self.assertEqual(tuple(tractInfo.getPatchInnerDimensions()), (xInnerDim, yInnerDim))
                self.assertEqual(len(skyMap), self.numTracts)
                self.assertNotEqual(skyMap, defaultSkyMap)

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
                pixelPos = geom.Point2D(x, y)
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
        self.assertEqual(skyMap, unpickledSkyMap)
        for tractInfo, unpickledTractInfo in zip(skyMap, unpickledSkyMap):
            self.assertUnpickledTractInfo(unpickledTractInfo, tractInfo, skyMap.config.patchBorder)

    def testTractSeparation(self):
        """Confirm that each sky tract has the proper distance to other tracts
        """
        if self.neighborAngularSeparation is None:
            self.skipTest("Not testing angular separation for %s: neighborAngularSeparation is None" %
                          (self.SkyMapClass.__name__,))
        skyMap = self.getSkyMap()
        for tractId, tractInfo in enumerate(skyMap):
            self.assertEqual(tractInfo.getId(), tractId)

            ctrCoord = tractInfo.getCtrCoord()
            distList = []
            for tractInfo1 in skyMap:
                otherCtrCoord = tractInfo1.getCtrCoord()
                distList.append(ctrCoord.separation(otherCtrCoord))
            distList.sort()
            self.assertEqual(distList[0], 0.0)
            for dist in distList[1:self.numNeighbors]:
                self.assertAnglesAlmostEqual(dist, self.neighborAngularSeparation)

    def testFindPatchList(self):
        """Test TractInfo.findPatchList
        """
        skyMap = self.getSkyMap()
        # pick two arbitrary tracts
        for tractId in np.random.choice(len(skyMap), 2):
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
        # pick 3 arbitrary tracts
        for tractId in np.random.choice(len(skyMap), 3):
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
            opposite = geom.SpherePoint(coord.getLongitude() + 12*geom.hours, -1*coord.getLatitude())
            self.assertFalse(tract.contains(opposite))

    def testTractInfoGetPolygon(self):
        skyMap = self.getSkyMap()
        for tractInfo in skyMap:
            centerCoord = tractInfo.getCtrCoord()
            self.assertPolygonOk(polygon=tractInfo.getInnerSkyPolygon(),
                                 vertexList=tractInfo.getVertexList(),
                                 centerCoord=centerCoord)
            self.assertBBoxPolygonOk(polygon=tractInfo.getOuterSkyPolygon(),
                                     bbox=tractInfo.getBBox(), wcs=tractInfo.getWcs())

    def testPatchInfoGetPolygon(self):
        skyMap = self.getSkyMap()
        numPatches = skyMap[0].getNumPatches()

        def getIndices(numItems):
            """Return up to 3 indices for testing"""
            if numItems > 2:
                return (0, 1, numItems-1)
            elif numItems > 1:
                return (0, 1)
            return (0,)

        for tractInfo in skyMap:
            wcs = tractInfo.getWcs()
            for patchInd in itertools.product(getIndices(numPatches[0]), getIndices(numPatches[1])):
                with self.subTest(patchInd=patchInd):
                    patchInfo = tractInfo.getPatchInfo(patchInd)
                    self.assertBBoxPolygonOk(polygon=patchInfo.getInnerSkyPolygon(tractWcs=wcs),
                                             bbox=patchInfo.getInnerBBox(), wcs=wcs)
                    self.assertBBoxPolygonOk(polygon=patchInfo.getOuterSkyPolygon(tractWcs=wcs),
                                             bbox=patchInfo.getOuterBBox(), wcs=wcs)

    def testDm14809(self):
        """Generic version of test that DM-14809 has been fixed"""
        checkDm14809(self, self.getSkyMap())

    def testNumbering(self):
        """Check the numbering of tracts matches the indexing"""
        skymap = self.getSkyMap()
        expect = list(range(len(skymap)))
        got = [tt.getId() for tt in skymap]
        self.assertEqual(got, expect)

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

    def assertBBoxPolygonOk(self, polygon, bbox, wcs):
        """Assert that an on-sky polygon from a pixel bbox
        covers the expected region.

        Parameters
        ----------
        polygon : `lsst.sphgeom.ConvexPolygon`
            On-sky polygon
        vertexList : `iterable` of `lsst.geom.SpherePoint`
            Vertices of polygon
        centerCoord : `lsst.geom.SpherePoint`
            A coord approximately in the center of the region
        """
        bboxd = geom.Box2D(bbox)
        centerPixel = bboxd.getCenter()
        centerCoord = wcs.pixelToSky(centerPixel)
        skyCorners = getCornerCoords(wcs=wcs, bbox=bbox)
        self.assertPolygonOk(polygon=polygon, vertexList=skyCorners, centerCoord=centerCoord)

    def assertPolygonOk(self, polygon, vertexList, centerCoord):
        """Assert that an on-sky polygon from covers the expected region.

        Parameters
        ----------
        polygon : `lsst.sphgeom.ConvexPolygon`
            On-sky polygon
        vertexList : `iterable` of `lsst.geom.SpherePoint`
            Vertices of polygon
        centerCoord : `lsst.geom.SpherePoint`
            A coord approximately in the center of the region
        """
        shiftAngle = 0.01*geom.arcseconds
        self.assertTrue(polygon.contains(centerCoord.getVector()))
        for vertex in vertexList:
            bearingToCenter = vertex.bearingTo(centerCoord)
            cornerShiftedIn = vertex.offset(bearing=bearingToCenter, amount=shiftAngle)
            cornerShiftedOut = vertex.offset(bearing=bearingToCenter, amount=-shiftAngle)
            self.assertTrue(polygon.contains(cornerShiftedIn.getVector()))
            self.assertFalse(polygon.contains(cornerShiftedOut.getVector()))


def getCornerCoords(wcs, bbox):
    """Return the coords of the four corners of a bounding box
    """
    cornerPosList = geom.Box2D(bbox).getCorners()
    return wcs.pixelToSky(cornerPosList)
