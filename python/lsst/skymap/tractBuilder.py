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
__all__ = ["tractBuilderRegistry",
           "BaseTractBuilderConfig", "BaseTractBuilder",
           "LegacyTractBuilderConfig", "LegacyTractBuilder",
           "CellTractBuilderConfig", "CellTractBuilder"]

import abc
import numbers
import struct
from collections.abc import Iterable

import lsst.pex.config as pexConfig
import lsst.geom as geom
from .patchInfo import PatchInfo
from .detail import Index2D


class BaseTractBuilderConfig(pexConfig.Config):
    """Configuration that is to be shared amongst all tract builders."""
    pass


class BaseTractBuilder(metaclass=abc.ABCMeta):
    """Base class for algorithms that define patches within the tract.

    Parameters
    ----------
    config : `lsst.pexConfig.Config`
        Input for configuring the algorithm
    """
    def __init__(self, config):
        self.config = config

    def setupPatches(self, minBBox, wcs):
        """Set up the patches of a particular size in a tract.

        We grow the tract bounding box to hold an exact multiple of
        the desired size (patchInnerDimensions or
        numCellsPerPatchInner*cellInnerDimensions), while keeping
        the center roughly the same.  We return the final tract
        bounding box, and the number of patches in each dimension
        (as an Index2D).

        Parameters
        ----------
        minBBox : `lsst.geom.Box2I`
            Minimum bounding box for tract.
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs object.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            final bounding box, number of patches.
        numPatches : `lsst.skymap.Index2D`
        """
        bbox = geom.Box2I(minBBox)
        bboxMin = bbox.getMin()
        bboxDim = bbox.getDimensions()
        numPatchesList = [0, 0]
        for i, innerDim in enumerate(self._patchInnerDimensions):
            num = (bboxDim[i] + innerDim - 1) // innerDim  # round up
            deltaDim = (innerDim*num) - bboxDim[i]
            if deltaDim > 0:
                bboxDim[i] = innerDim * num
                bboxMin[i] -= deltaDim // 2
            numPatchesList[i] = num
        numPatches = Index2D(*numPatchesList)
        bbox = geom.Box2I(bboxMin, bboxDim)
        self._numPatches = numPatches
        # The final tract BBox starts at zero.
        self._tractBBox = geom.Box2I(geom.Point2I(0, 0), bbox.getDimensions())
        self._initialized = True

        return bbox, numPatches

    def getPatchBorder(self):
        return self._patchBorder

    @abc.abstractmethod
    def getPatchInfo(self, index, tractWcs):
        """Return information for the specified patch.

        Parameters
        ----------
        index : `lsst.skymap.Index2D` or `~collections.abc.Iterable` of 2 `int`
            Index of patch, as Index2D or pair of ints;
            or a sequential index as returned by getSequentialPatchIndex;
            negative values are not supported.
        tractWcs : `lsst.afw.geom.SkyWcs`
            WCS associated with the tract.

        Returns
        -------
        result : `lsst.skymap.PatchInfo`
            The patch info for that index.

        Raises
        ------
        IndexError
            Raised if index is out of range.
        """
        raise NotImplementedError("Must be implemented by a subclass")

    def getPatchInnerDimensions(self):
        """Get dimensions of inner region of the patches (all are the same)
        """
        return self._patchInnerDimensions

    def getSequentialPatchIndex(self, patchInfo):
        """Return a single integer that uniquely identifies
        the given patch within this tract.

        Parameters
        ----------
        patchInfo : `lsst.skymap.PatchInfo`

        Returns
        -------
        sequentialIndex : `int`
        """
        index = patchInfo.getIndex()
        return self.getSequentialPatchIndexFromPair(index)

    def getSequentialPatchIndexFromPair(self, index):
        """Return a single integer that uniquely identifies
        the patch index within the tract.

        Parameters
        ----------
        index : `lsst.skymap.Index2D` or `~collections.abc.Iterable` of 2 `int`

        Returns
        -------
        sequentialIndex : `int`
        """
        if isinstance(index, Index2D):
            _index = index
        else:
            if not isinstance(index, Iterable):
                raise ValueError("Input index is not an iterable.")
            if len(index) != 2:
                raise ValueError("Input index does not have two values.")
            _index = Index2D(*index)
        nx, ny = self._numPatches
        return nx*_index.y + _index.x

    def getPatchIndexPair(self, sequentialIndex):
        """Convert sequential index into patch index (x,y) pair.

        Parameters
        ----------
        sequentialIndex : `int`

        Returns
        -------
        x, y : `lsst.skymap.Index2D`
        """
        nx, ny = self._numPatches
        x = sequentialIndex % nx
        y = sequentialIndex // nx
        return Index2D(x=x, y=y)

    @abc.abstractmethod
    def getPackedConfig(self, config):
        """Get a packed config suitable for using in a sha1.

        Parameters
        ----------
        config : `lsst.skymap.BaseTractBuilderConfig`

        Returns
        -------
        configPacked : `bytes`
        """
        raise NotImplementedError("Must be implemented by a subclass")


class LegacyTractBuilderConfig(BaseTractBuilderConfig):
    patchInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of patches (x,y pixels)",
        dtype=int,
        length=2,
        default=(4000, 4000),
    )
    patchBorder = pexConfig.Field(
        doc="border between patch inner and outer bbox (pixels)",
        dtype=int,
        default=100,
    )


class LegacyTractBuilder(BaseTractBuilder):
    ConfigClass = LegacyTractBuilderConfig

    def __init__(self, config):
        super().__init__(config)

        self._patchInnerDimensions = geom.Extent2I(*(val
                                                     for val in config.patchInnerDimensions))
        self._patchBorder = config.patchBorder
        self._initialized = False

    def getPatchInfo(self, index, tractWcs):
        # This should always be initialized
        if not self._initialized:
            raise RuntimeError("Programmer error; this should always be initialized.")
        if isinstance(index, Index2D):
            _index = index
        else:
            if isinstance(index, numbers.Number):
                _index = self.getPatchIndexPair(index)
            else:
                _index = Index2D(*index)
        if (not 0 <= _index.x < self._numPatches.x) \
           or (not 0 <= _index.y < self._numPatches.y):
            raise IndexError("Patch index %s is not in range [0-%d, 0-%d]" %
                             (_index, self._numPatches.x - 1, self._numPatches.y - 1))
        innerMin = geom.Point2I(*[_index[i] * self._patchInnerDimensions[i] for i in range(2)])
        innerBBox = geom.Box2I(innerMin, self._patchInnerDimensions)
        if not self._tractBBox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" %
                (_index, innerBBox, self._tractBBox))
        outerBBox = geom.Box2I(innerBBox)
        outerBBox.grow(self.getPatchBorder())
        outerBBox.clip(self._tractBBox)
        return PatchInfo(
            index=_index,
            innerBBox=innerBBox,
            outerBBox=outerBBox,
            sequentialIndex=self.getSequentialPatchIndexFromPair(_index),
            tractWcs=tractWcs
        )

    def getPackedConfig(self, config):
        subConfig = config.tractBuilder[config.tractBuilder.name]
        configPacked = struct.pack(
            "<iiidd3sd",
            subConfig.patchInnerDimensions[0],
            subConfig.patchInnerDimensions[1],
            subConfig.patchBorder,
            config.tractOverlap,
            config.pixelScale,
            config.projection.encode('ascii'),
            config.rotation
        )

        return configPacked


class CellTractBuilderConfig(BaseTractBuilderConfig):
    cellInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of cells (x,y pixels)",
        dtype=int,
        length=2,
        default=(150, 150),
    )
    cellBorder = pexConfig.Field(
        doc="Border between cell inner and outer bbox (pixels)",
        dtype=int,
        default=50,
    )
    numCellsPerPatchInner = pexConfig.Field(
        doc="Number of cells per inner patch.",
        dtype=int,
        default=20,
    )
    numCellsInPatchBorder = pexConfig.Field(
        doc="Number of cells in the patch border (outside the inner patch region).",
        dtype=int,
        default=1,
    )

    def validate(self):
        if len(self.cellInnerDimensions) != 2:
            raise ValueError("cellInnerDimensions must be 2 ints.")

        if self.cellInnerDimensions[0] != self.cellInnerDimensions[1]:
            raise ValueError("cellInnerDimensions must be equal (for square cells).")


class CellTractBuilder(BaseTractBuilder):
    ConfigClass = CellTractBuilderConfig

    def __init__(self, config):
        super().__init__(config)

        self._cellInnerDimensions = geom.Extent2I(*(val
                                                    for val in config.cellInnerDimensions))
        self._cellBorder = config.cellBorder
        self._numCellsPerPatchInner = config.numCellsPerPatchInner
        self._numCellsInPatchBorder = config.numCellsInPatchBorder
        self._patchInnerDimensions = geom.Extent2I(*(val*self._numCellsPerPatchInner
                                                     for val in config.cellInnerDimensions))
        # The patch border is the number of cells in the border + the cell
        # border.
        self._patchBorder = config.numCellsInPatchBorder*config.cellInnerDimensions[0] + self._cellBorder
        self._initialized = False

    def getPatchInfo(self, index, tractWcs):
        # This should always be initialized
        if not self._initialized:
            raise RuntimeError("Programmer error; this should always be initialized.")
        if isinstance(index, Index2D):
            _index = index
        else:
            if isinstance(index, numbers.Number):
                _index = self.getPatchIndexPair(index)
            else:
                _index = Index2D(*index)
        if (not 0 <= _index.x < self._numPatches.x) \
           or (not 0 <= _index.y < self._numPatches.y):
            raise IndexError("Patch index %s is not in range [0-%d, 0-%d]" %
                             (_index, self._numPatches.x - 1, self._numPatches.y - 1))
        innerMin = geom.Point2I(*[_index[i]*self._patchInnerDimensions[i] for i in range(2)])
        innerBBox = geom.Box2I(innerMin, self._patchInnerDimensions)
        if not self._tractBBox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" %
                (_index, innerBBox, self._tractBBox))
        outerBBox = geom.Box2I(innerBBox)
        outerBBox.grow(self.getPatchBorder())
        # We do not clip the patch for cell-based tracts.
        return PatchInfo(
            index=_index,
            innerBBox=innerBBox,
            outerBBox=outerBBox,
            sequentialIndex=self.getSequentialPatchIndexFromPair(_index),
            tractWcs=tractWcs,
            cellInnerDimensions=self._cellInnerDimensions,
            cellBorder=self._cellBorder,
            numCellsPerPatchInner=self._numCellsPerPatchInner,
            numCellsInPatchBorder=self._numCellsInPatchBorder
        )

    def getPackedConfig(self, config):
        subConfig = config.tractBuilder[config.tractBuilder.name]
        configPacked = struct.pack(
            "<iiiiidd3sd",
            subConfig.cellInnerDimensions[0],
            subConfig.cellInnerDimensions[1],
            subConfig.cellBorder,
            subConfig.numCellsPerPatchInner,
            subConfig.numCellsInPatchBorder,
            config.tractOverlap,
            config.pixelScale,
            config.projection.encode('ascii'),
            config.rotation
        )

        return configPacked


tractBuilderRegistry = pexConfig.makeRegistry(
    doc="A registry of Tract Builders (subclasses of BaseTractBuilder)",
)

tractBuilderRegistry.register("legacy", LegacyTractBuilder)
tractBuilderRegistry.register("cells", CellTractBuilder)
