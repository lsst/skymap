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

__all__ = ["PatchInfo"]

import numbers
from collections.abc import Iterable

from lsst.geom import Extent2I, Point2I, Box2I
from .detail import makeSkyPolygonFromBBox, Index2D
from .cellInfo import CellInfo


class PatchInfo:
    """Information about a patch within a tract of a sky map.

    If cellInnerDimensions and cellBorder are set then the patch
    will be gridded with cells.

    See `TractInfo` for more information.

    Parameters
    ----------
    index : `lsst.skymap.Index2D`
        x,y index of patch (a pair of ints)
    innerBBox : `lsst.geom.Box2I`
        inner bounding box
    outerBBox : `lsst.geom.Box2I`
        inner bounding box
    sequentialIndex : `int`
        Patch sequential index
    tractWcs : `lsst.afw.geom.SkyWcs`
        Tract WCS object.
    cellInnerDimensions : `Iterable` [`int`, `int`] or `lsst.geom.Extent2I`, optional
        Inner dimensions of each cell (x,y pixels).
    cellBorder : `int`, optional
        Cell border size (pixels).
    numCellsPerPatchInner : `int`, optional
        Number of cells per inner patch region.
    numCellsInPatchBorder : `int`, optional
        Number of cells in the patch border.
    """

    def __init__(self, index, innerBBox, outerBBox, sequentialIndex,
                 tractWcs,
                 cellInnerDimensions=(0, 0), cellBorder=0,
                 numCellsPerPatchInner=0, numCellsInPatchBorder=0):
        self._index = index
        self._sequentialIndex = sequentialIndex
        self._innerBBox = innerBBox
        self._outerBBox = outerBBox
        self._wcs = tractWcs
        if not outerBBox.contains(innerBBox):
            raise RuntimeError("outerBBox=%s does not contain innerBBox=%s" % (outerBBox, innerBBox))
        if not isinstance(cellInnerDimensions, (Iterable, Extent2I)):
            raise ValueError("Input cellInnerDimensions is not an iterable.")
        if len(cellInnerDimensions) != 2:
            raise ValueError("Input cellInnerDimensions does not have two values.")
        self._cellInnerDimensions = Extent2I(*cellInnerDimensions)
        self._cellBorder = cellBorder
        self._numCellsInPatchBorder = numCellsInPatchBorder
        if numCellsPerPatchInner == 0:
            self._numCells = Index2D(x=0, y=0)
        else:
            # There are numCellsInPatchBorder extra boundary cell on each side
            self._numCells = Index2D(x=numCellsPerPatchInner + 2*numCellsInPatchBorder,
                                     y=numCellsPerPatchInner + 2*numCellsInPatchBorder)

    def getIndex(self):
        """Return patch index: a tuple of (x, y)

        Returns
        -------
        result : `lsst.skymap.Index2D`
            Patch index (x, y).
        """
        return self._index

    index = property(getIndex)

    def getSequentialIndex(self):
        """Return patch sequential index.

        Returns
        -------
        result : `int`
            Sequential patch index.
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

    def getNumCells(self):
        """Get the number of cells in x, y.

        May return (0, 0) if no cells are defined.

        Returns
        -------
        result : `lsst.skymap.Index2D`
            The number of cells in x, y.
        """
        return self._numCells

    num_cells = property(getNumCells)

    def getCellBorder(self):
        return self._cellBorder

    cell_border = property(getCellBorder)

    def getCellInfo(self, index):
        """Return information for the specified cell.

        Parameters
        ----------
        index : `lsst.skymap.Index2D` or `int`
            Index of cell, as `Index2D`, or `Iterable` [`int`, `int`];
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
        if self._numCells.x == 0 or self._numCells.y == 0:
            raise IndexError("Patch does not contain cells.")
        if isinstance(index, Index2D):
            _index = index
        else:
            if isinstance(index, numbers.Number):
                _index = self.getCellIndexPair(index)
            else:
                _index = Index2D(*index)
        if (not 0 <= _index.x < self._numCells.x) \
                or (not 0 <= _index.y < self._numCells.y):
            raise IndexError("Cell index %s is not in range [0-%d, 0-%d]" %
                             (_index, self._numCells.x - 1, self._numCells.y - 1))
        # We offset the index by numCellsInPatchBorder because the cells
        # start outside the inner dimensions.
        # The cells are defined relative to the patch bounding box (within the tract).
        patchInnerBBox = self.getInnerBBox()
        innerMin = Point2I(*[(_index[i] - self._numCellsInPatchBorder)*self._cellInnerDimensions[i]
                             + patchInnerBBox.getBegin()[i]
                             for i in range(2)])

        innerBBox = Box2I(innerMin, self._cellInnerDimensions)
        outerBBox = Box2I(innerBBox)
        outerBBox.grow(self._cellBorder)

        return CellInfo(
            index=_index,
            innerBBox=innerBBox,
            outerBBox=outerBBox,
            sequentialIndex=self.getSequentialCellIndexFromPair(_index),
            tractWcs=self._wcs
        )

    def getCellInnerDimensions(self):
        """Get dimensions of inner region of the cells (all are the same)
        """
        return self._cellInnerDimensions

    cell_inner_dimensions = property(getCellInnerDimensions)

    def getSequentialCellIndex(self, cellInfo):
        """Return a single integer that uniquely identifies
        the given cell within this patch.

        Parameters
        ----------
        cellInfo : `lsst.skymap.CellInfo`

        Returns
        -------
        sequentialIndex : `int`

        Raises
        ------
        IndexError
            If index is out of range.
        """
        index = cellInfo.getIndex()
        return self.getSequentialCellIndexFromPair(index)

    def getSequentialCellIndexFromPair(self, index):
        """Return a single integer that uniquely identifies
        the given cell within this patch.

        Parameters
        ----------
        index : `lsst.skymap.Index2D`

        Returns
        -------
        sequentialIndex : `int`

        Raises
        ------
        IndexError
            If index is out of range.
        """
        if isinstance(index, Index2D):
            _index = index
        else:
            _index = Index2D(*index)
        nx, ny = self.getNumCells()
        return nx*_index.y + _index.x

    def getCellIndexPair(self, sequentialIndex):
        """Convert a sequential index into an index pair.

        Parameters
        ----------
        sequentialIndex : `int`

        Returns
        -------
        x, y : `lsst.skymap.Index2D`

        Raises
        ------
        IndexError
            If index is out of range.
        """
        if self._numCells.x == 0 or self._numCells.y == 0:
            raise IndexError("Patch does not contain cells.")

        nx, ny = self.getNumCells()
        x = sequentialIndex % nx
        y = sequentialIndex // nx
        return Index2D(x=x, y=y)

    def __iter__(self):
        xNum, yNum = self.getNumCells()
        for y in range(yNum):
            for x in range(xNum):
                yield self.getCellInfo(Index2D(x=x, y=y))

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
