# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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

from lsst.pex.config import ListField
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo

__all__ = ["DiscreteSkyMap"]

class DiscreteSkyMapConfig(BaseSkyMap.ConfigClass):
    """Configuration for the DiscreteSkyMap"""
    raList = ListField(dtype=float, default=[], doc="Right Ascensions of tracts (ICRS, degrees)")
    decList = ListField(dtype=float, default=[], doc="Declinations of tracts (ICRS, degrees)")
    radiusList = ListField(dtype=float, default=[], doc="Radii of tracts (degrees)")

    def validate(self):
        super(DiscreteSkyMapConfig, self).validate()
        if len(self.radiusList) != len(self.raList):
            raise ValueError("Number of radii (%d) and RAs (%d) do not match" %
                             (len(self.radiusList), len(self.raList)))
        if len(self.radiusList) != len(self.decList):
            raise ValueError("Number of radii (%d) and Decs (%d) do not match" %
                             (len(self.radiusList), len(self.decList)))

class DiscreteTractInfo(TractInfo):
    """Tract for DiscreteSkyMap

    A tract is placed at the nominated coordinates, with the nominated radius.
    The tracts are square (i.e., the radius is really a half-size).
    """
    def __init__(self, ident, patchInnerDimensions, patchBorder, ctrCoord, tractOverlap, wcs, radius):
        vertexList = []
        self._radius = radius
        super(DiscreteTractInfo, self).__init__(ident, patchInnerDimensions, patchBorder, ctrCoord,
                                                vertexList, tractOverlap, wcs)

    def _minimumBoundingBox(self, wcs):
        """The minimum bounding box is calculated using the nominated radius"""
        bbox = afwGeom.Box2D()
        for i in range(4):
            coord = self._ctrCoord.clone()
            coord.offset(i * 90 * afwGeom.degrees, self._radius)
            pixPos = wcs.skyToPixel(coord)
            bbox.include(pixPos)
        return bbox


class DiscreteSkyMap(BaseSkyMap):
    """Discrete sky map pixelization.

    We put a square Tract at each of the nominated coordinates.
    """
    ConfigClass = DiscreteSkyMapConfig
    _version = (1, 0) # for pickle
    numAngles = 4 # Number of angles for vertices

    def __init__(self, config=None, version=0):
        """Constructor

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        @param[in] version: software version of this class, to retain compatibility with old instances
        """
        super(DiscreteSkyMap, self).__init__(config)
        self._version = version
        self._numTracts = len(self.config.radiusList)
        self._tractCache = [None] * self._numTracts
        self._tractInfo = None # We shouldn't need this; we will generate tracts on demand

    def __reduce__(self):
        """To support pickling"""
        return (self.__class__, (self.config, self._version))

    def __getitem__(self, index):
        """Get the TractInfo for a particular index

        The tract is returned from a cache, if available, otherwise generated
        on the fly.
        """
        if index < 0 or index > self._numTracts:
            raise IndexError("Index out of range: %d vs %d" % (index, self._numTracts))
        if self._tractCache[index] is not None:
            return self._tractCache[index]
        center = IcrsCoord(self.config.raList[index] * afwGeom.degrees,
                           self.config.decList[index] * afwGeom.degrees)
        radius = self.config.radiusList[index]
        wcs = self._wcsFactory.makeWcs(crPixPos=afwGeom.Point2D(0,0), crValCoord=center)
        tract = DiscreteTractInfo(index, self.config.patchInnerDimensions, self.config.patchBorder, center,
                                  self.config.tractOverlap * afwGeom.degrees, wcs, radius * afwGeom.degrees)
        self._tractCache[index] = tract
        return tract

    def __iter__(self):
        """Iterator over tracts"""
        for i in xrange(self._numTracts):
            yield self[i]

    def __len__(self):
        """Length is number of tracts"""
        return self._numTracts


