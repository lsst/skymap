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

import numbers

from lsst.geom import Extent2I, Point2I, Box2I
from .detail import makeSkyPolygonFromBBox
from .cellInfo import CellInfo


class PatchInfo:
    """Information about a patch within a tract of a sky map.

    If cellInnerDimensions and cellBorder are set then the patch
    will be gridded with cells.

    See `TractInfo` for more information.

    Parameters
    ----------
    index : `tuple` of `int`
        x,y index of patch (a pair of ints)
    innerBBox : `lsst.geom.Box2I`
        inner bounding box
    outerBBox : `lsst.geom.Box2I`
        inner bounding box
    cellInnerDimensions : `lsst.geom.Extent2I`, optional
        Inner dimensions of each cell (x,y pixels).
    cellBorder : `int`, optional
        Cell border size (pixels).
    numCellsPerPatchInner : `int`, optional
        Number of cells per inner patch region.
    """

    def __init__(self, index, innerBBox, outerBBox,
                 cellInnerDimensions=Extent2I(0, 0), cellBorder=0,
                 numCellsPerPatchInner=0):
        self._index = index
        self._innerBBox = innerBBox
        self._outerBBox = outerBBox
        if not outerBBox.contains(innerBBox):
            raise RuntimeError("outerBBox=%s does not contain innerBBox=%s" % (outerBBox, innerBBox))
        self._cellInnerDimensions = cellInnerDimensions
        self._cellBorder = cellBorder
        if numCellsPerPatchInner == 0:
            self._numCells = Extent2I(0, 0)
        else:
            # There is one extra boundary cell on each side
            self._numCells = Extent2I(numCellsPerPatchInner + 2,
                                      numCellsPerPatchInner + 2)

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

        Parameters
        ----------
        tractWcs : `lsst.afw.image.SkyWcs`
            WCS for the associated tract.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The inner sky region.
        """
        return makeSkyPolygonFromBBox(bbox=self.getInnerBBox(), wcs=tractWcs)

    def getOuterSkyPolygon(self, tractWcs):
        """Get the outer on-sky region.

        Parameters
        ----------
        tractWcs : `lsst.afw.image.SkyWcs`
            WCS for the associated tract.

        Returns
        -------
        result : `lsst.sphgeom.ConvexPolygon`
            The outer sky region.
        """
        return makeSkyPolygonFromBBox(bbox=self.getOuterBBox(), wcs=tractWcs)

    def getNumCells(self):
        """Get the number of cells in x, y.

        May return (0, 0) if no cells are defined.

        Returns
        -------
        result : `lsst.geom.Extent2I`
            The number of cells in x, y.
        """
        return self._numCells

    def getCellBorder(self):
        return self._cellBorder

    def getCellInfo(self, index):
        """Return information for the specified cell.

        Parameters
        ----------
        index : `tuple` of `int`
            Index of cell, as a pair of ints;
            or a sequential index as returned by getSequentialCellIndex;
            negative values are not supported.

        Returns
        -------
        result : `lsst.skymap.CellInfo`
            The cell info for that index.

        Raises
        ------
        IndexError
            If index is out of range.
        """
        if isinstance(index, numbers.Number):
            index = self.getCellIndexPair(index)
        if (not 0 <= index[0] < self._numCells[0]) \
                or (not 0 <= index[1] < self._numCells[1]):
            raise IndexError("Cell index %s is not in range [0-%d, 0-%d]" %
                             (index, self._numCells[0] - 1, self._numCells[1] - 1))
        # We offset the index by 1 because the cells start outside the inner
        # dimensions.
        innerMin = Point2I(*[(index[i] - 1)*self._cellInnerDimensions[i]
                             for i in range(2)])


        innerBBox = Box2I(innerMin, self._cellInnerDimensions)
        outerBBox = Box2I(innerBBox)
        outerBBox.grow(self._cellBorder)

        return CellInfo(
            index=index,
            innerBBox=innerBBox,
            outerBBox=outerBBox
        )

    def getCellInnerDimensions(self):
        """Get dimensions of inner region of the cells (all are the same)
        """
        return self._cellInnerDimensions

    def getSequentialCellIndex(self, cellInfo):
        """Return a single integer that uniquely identifies the given
        cell within this patch.
        """
        x, y = cellInfo.getIndex()
        nx, ny = self.getNumCells()
        return nx*y + x

    def getCellIndexPair(self, sequentialIndex):
        nx, ny = self.getNumCells()
        x = sequentialIndex % nx
        y = (sequentialIndex - x) // nx
        return (x, y)

    def __iter__(self):
        xNum, yNum = self.getNumCells()
        for y in range(yNum):
            for x in range(xNum):
                yield self.getCellInfo((x, y))

    def __len__(self):
        xNum, yNum = self.getNumCells()
        return xNum*yNum

    def __getitem__(self, index):
        return self.getCellInfo(index)

    def __eq__(self, rhs):
        return (self.getIndex() == rhs.getIndex()) \
            and (self.getInnerBBox() == rhs.getInnerBBox()) \
            and (self.getOuterBBox() == rhs.getOuterBBox()) \
            and (self.getNumCells() == rhs.getNumCells()) \
            and (self.getCellBorder() == rhs.getCellBorder())

    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    def __str__(self):
        return "PatchInfo(index=%s)" % (self.getIndex(),)

    def __repr__(self):
        if self.getNumCells()[0] > 0:
            return ("PatchInfo(index=%s, innerBBox=%s, outerBBox=%s, cellInnerDimensions=%s, "
                    "cellBorder=%s, numCellsPerPatchInner=%s)") % \
                (self.getIndex(), self.getInnerBBox(), self.getOuterBBox(),
                 self.getCellInnerDimensions(), self.getCellBorder(),
                 self.getNumCells()[0])
        else:
            return "PatchInfo(index=%s, innerBBox=%s, outerBBox=%s)" % \
                (self.getIndex(), self.getInnerBBox(), self.getOuterBBox())
