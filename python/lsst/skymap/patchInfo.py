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

__all__ = ["PatchInfo", "makeSkyPolygonFromBBox"]

from lsst.sphgeom import ConvexPolygon
from lsst.geom import Box2D


def makeSkyPolygonFromBBox(bbox, wcs):
    """Make an on-sky polygon from a bbox and a SkyWcs

    Parameters
    ----------
    bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
        Bounding box of region, in pixel coordinates
    wcs : `lsst.afw.geom.SkyWcs`
        Celestial WCS

    Returns
    -------
    polygon : `lsst.sphgeom.ConvexPolygon`
        On-sky region
    """
    pixelPoints = Box2D(bbox).getCorners()
    skyPoints = wcs.pixelToSky(pixelPoints)
    return ConvexPolygon.convexHull([sp.getVector() for sp in skyPoints])


class PatchInfo:
    """Information about a patch within a tract of a sky map.

    See `TractInfo` for more information.

    Parameters
    ----------
    index : `tuple` of `int`
        x,y index of patch (a pair of ints)
    innerBBox : `lsst.geom.Box2I`
        inner bounding box
    outerBBox : `lsst.geom.Box2I`
        inner bounding box
    """

    def __init__(self, index, innerBBox, outerBBox):
        self._index = index
        self._innerBBox = innerBBox
        self._outerBBox = outerBBox
        if not outerBBox.contains(innerBBox):
            raise RuntimeError("outerBBox=%s does not contain innerBBox=%s" % (outerBBox, innerBBox))

    def getIndex(self):
        """Return patch index: a tuple of (x, y)

        Returns
        -------
        result : `tuple` of `int`
            Patch index (x, y).
        """
        return self._index

    def getInnerBBox(self):
        """Get inner bounding box.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The inner bounding Box.
        """
        return self._innerBBox

    def getOuterBBox(self):
        """Get outer bounding box.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The outer bounding Box.
        """
        return self._outerBBox

    def getInnerSkyPolygon(self, tractWcs):
        """Get the inner on-sky region.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The inner sky region.
        """
        return makeSkyPolygonFromBBox(bbox=self.getInnerBBox(), wcs=tractWcs)

    def getOuterSkyPolygon(self, tractWcs):
        """Get the outer on-sky region.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The outer sky region.
        """
        return makeSkyPolygonFromBBox(bbox=self.getOuterBBox(), wcs=tractWcs)

    def __eq__(self, rhs):
        return (self.getIndex() == rhs.getIndex()) \
            and (self.getInnerBBox() == rhs.getInnerBBox()) \
            and (self.getOuterBBox() == rhs.getOuterBBox())

    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    def __str__(self):
        return "PatchInfo(index=%s)" % (self.getIndex(),)

    def __repr__(self):
        return "PatchInfo(index=%s, innerBBox=%s, outerBBox=%s)" % \
            (self.getIndex(), self.getInnerBBox(), self.getOuterBBox())
