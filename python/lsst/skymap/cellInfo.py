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

__all__ = ["CellInfo"]

from .detail import makeSkyPolygonFromBBox


class CellInfo:
    """Information about a cell within a patch of a tract of a sky map.

    See `PatchInfo` and `TractInfo` for more information.

    Parameters
    ----------
    index : `lsst.skymap.Index2D`
        x,y index of a cell (a pair of ints)
    innerBBox : `lsst.geom.Box2I`
        Inner bounding box.
    outerBBox : `lsst.geom.Box2I`
        Outer bounding box.
    sequentialIndex : `int`
        Cell sequential index.
    tractWcs : `lsst.afw.geom.SkyWcs`
        Tract WCS object.
    """
    def __init__(self, index, innerBBox, outerBBox, sequentialIndex, tractWcs):
        self._index = index
        self._sequentialIndex = sequentialIndex
        self._innerBBox = innerBBox
        self._outerBBox = outerBBox
        self._wcs = tractWcs
        if not outerBBox.contains(innerBBox):
            raise RuntimeError("outerBBox=%s does not contain innerBBox=%s" % (outerBBox, innerBBox))

    def getIndex(self):
        """Return cell index: a tuple of (x, y)

        Returns
        -------
        result : `lsst.skymap.Index2D`
            Patch index (x, y).
        """
        return self._index

    index = property(getIndex)

    def getSequentialIndex(self):
        """Return cell sequential index.

        Returns
        -------
        result : `int`
            Sequential cell index.
        """
        return self._sequentialIndex

    sequential_index = property(getSequentialIndex)

    def getWcs(self):
        """Return the associated tract wcs

        Returns
        -------
        wcs : `lsst.afw.geom.SkyWcs`
            Tract WCS.
        """
        return self._wcs

    wcs = property(getWcs)

    def getInnerBBox(self):
        """Get inner bounding box.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The inner bounding Box.
        """
        return self._innerBBox

    inner_bbox = property(getInnerBBox)

    def getOuterBBox(self):
        """Get outer bounding box.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The outer bounding Box.
        """
        return self._outerBBox

    outer_bbox = property(getOuterBBox)

    def getInnerSkyPolygon(self, tractWcs=None):
        """Get the inner on-sky region.

        Parameters
        ----------
        tractWcs : `lsst.afw.image.SkyWcs`, optional
            WCS for the associated tract.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The inner sky region.
        """
        _tractWcs = tractWcs if tractWcs is not None else self._wcs
        return makeSkyPolygonFromBBox(bbox=self.getInnerBBox(), wcs=_tractWcs)

    @property
    def inner_sky_polygon(self):
        return self.getInnerSkyPolygon()

    def getOuterSkyPolygon(self, tractWcs=None):
        """Get the outer on-sky region.

        Parameters
        ----------
        tractWcs : `lsst.afw.image.SkyWcs`, optional
            WCS for the associated tract.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The outer sky region.
        """
        _tractWcs = tractWcs if tractWcs is not None else self._wcs
        return makeSkyPolygonFromBBox(bbox=self.getOuterBBox(), wcs=_tractWcs)

    @property
    def outer_sky_polygon(self):
        return self.getOuterSkyPolygon()

    def __eq__(self, rhs):
        return (self.getIndex() == rhs.getIndex()) \
            and (self.getInnerBBox() == rhs.getInnerBBox()) \
            and (self.getOuterBBox() == rhs.getOuterBBox())

    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    def __str__(self):
        return "CellInfo(index=%s)" % (self.getIndex(),)

    def __repr__(self):
        return "CellInfo(index=%s, innerBBox=%s, outerBBox=%s)" % \
            (self.getIndex(), self.getInnerBBox(), self.getOuterBBox())
