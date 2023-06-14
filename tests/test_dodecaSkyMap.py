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
import math
import unittest

import numpy

import lsst.sphgeom
import lsst.geom as geom
import lsst.utils.tests

from lsst.skymap import DodecaSkyMap
from helper import skyMapTestCase


# dodecahedron properties
_Phi = (1.0 + math.sqrt(5.0)) / 2.0
_DihedralAngle = geom.Angle(2.0 * math.atan(_Phi), geom.radians)


class DodecaSkyMapTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        self.setAttributes(
            SkyMapClass=DodecaSkyMap,
            name="dodeca",
            numNeighbors=6,
            numTracts=12,
            neighborAngularSeparation=180*geom.degrees - _DihedralAngle,
        )

    def testSha1Compare(self):
        """Test that DodecaSkyMap's extra state is included in its hash."""
        defaultSkyMap = self.getSkyMap()
        config = self.getConfig()
        config.withTractsOnPoles = True
        skyMap = self.getSkyMap(config=config)
        self.assertNotEqual(skyMap, defaultSkyMap)

    def testFindTract(self):
        """Test findTract and tractInfo.findPatch
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
                dist = ctrCoord0.separation(otherCtrCoord)
                if abs(dist - self.neighborAngularSeparation) < geom.Angle(0.1, geom.degrees):
                    nbrTractList.append(otherTractInfo)
            self.assertEqual(len(nbrTractList), 5)

            for tractInfo1 in nbrTractList:
                tractId1 = tractInfo1.getId()
                ctrCoord1 = tractInfo1.getCtrCoord()
                vector1 = numpy.array(ctrCoord1.getVector())
                for tractInfo2 in nbrTractList[tractInfo1.getId():]:
                    dist = ctrCoord1.separation(tractInfo2.getCtrCoord())
                    if abs(dist - self.neighborAngularSeparation) > geom.Angle(0.1, geom.degrees):
                        continue
                    tractId2 = tractInfo2.getId()
                    ctrCoord2 = tractInfo2.getCtrCoord()
                    vector2 = numpy.array(ctrCoord2.getVector())

                    # Sky tracts 0, 1 and 2 form a triangle of nearest
                    # neighbors explore the boundary between tract 0 and tract
                    # 1 and also the boundary between tract 0 and tract 2.
                    for deltaFrac in (-0.001, 0.001):
                        isNearest0 = deltaFrac > 0.0

                        for exploreBoundary1 in (True, False):
                            # If exploreBoundary1, explore boundary between
                            # tract 0 and tract 1,
                            # else explore the boundary between tract 0 and
                            # tract 2.

                            if isNearest0:
                                expectedTractId = tractId0
                            elif exploreBoundary1:
                                expectedTractId = tractId1
                            else:
                                expectedTractId = tractId2

                            for farFrac in (0.0, 0.05, 0.3, (1.0/3.0) - 0.01):
                                # farFrac is the fraction of the tract center
                                # vector point whose boundary is not being
                                # explored; it must be less than 1/3;
                                # remFrac is the remaining fraction, which is
                                # divided between tract 0 and the tract whose
                                # boundary is being explored
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
                                testCoord = geom.SpherePoint(lsst.sphgeom.Vector3d(*testVector))
                                nearestTractInfo = skyMap.findTract(testCoord)
                                nearestTractId = nearestTractInfo.getId()

                                if expectedTractId != nearestTractId:
                                    nearestCtrCoord = nearestTractInfo.getCtrCoord()
                                    nearestVector = nearestCtrCoord.getVector()

                                    print("tractId0=%s; tractId1=%s; tractId2=%s; nearestTractId=%s" %
                                          (tractId0, tractId1, tractId2, nearestTractId))
                                    print("vector0=%s; vector1=%s; vector2=%s; nearestVector=%s" %
                                          (vector0, vector1, vector2, nearestVector))
                                    print("frac0=%s; frac1=%s; frac2=%s" % (frac0, frac1, frac2))
                                    print("testVector=", testVector)

                                    print("dist0=%s; dist1=%s; dist2=%s; nearDist=%s" % (
                                        testCoord.separation(ctrCoord0).asDegrees(),
                                        testCoord.separation(ctrCoord1).asDegrees(),
                                        testCoord.separation(ctrCoord2).asDegrees(),
                                        testCoord.separation(nearestCtrCoord).asDegrees(),
                                    ))
                                    self.fail("Expected nearest tractId=%s; got tractId=%s" %
                                              (expectedTractId, nearestTractId))

                                patchInfo = nearestTractInfo.findPatch(testCoord)
                                pixelInd = geom.Point2I(
                                    nearestTractInfo.getWcs().skyToPixel(testCoord))
                                self.assertTrue(patchInfo.getInnerBBox().contains(pixelInd))


def getCornerCoords(wcs, bbox):
    """Return the coords of the four corners of a bounding box
    """
    bbox = geom.Box2D(bbox)  # mak
    cornerPosList = (
        bbox.getMin(),
        geom.Point2D(bbox.getMaxX(), bbox.getMinY()),
        bbox.getMax(),
        geom.Point2D(bbox.getMinX(), bbox.getMaxY()),
    )
    return wcs.pixelToSky(cornerPosList)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
