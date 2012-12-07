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

import numpy
try:
    import healpy
except Exception, e:
    class DummyHealpy(object):
        """An object which blows up when we try to read it"""
        def __getattr__(self, name):
            raise RuntimeError("Was unable to import healpy: %s" % e)
    healpy = DummyHealpy()

from lsst.pex.config import Field
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo


def angToCoord(thetaphi):
    """Convert healpy's ang to an afw Coord

    The ang is provided as a single object, thetaphi, so the output
    of healpy functions can be directed to this function without
    additional translation.
    """
    return IcrsCoord(float(thetaphi[1])*afwGeom.radians, float(thetaphi[0] - 0.5*numpy.pi)*afwGeom.radians)

def coordToAng(coord):
    """Convert an afw Coord to a healpy ang (theta, phi)"""
    return (coord.getLatitude().asRadians() - 0.5*numpy.pi, coord.getLongitude().asRadians())

class HealpixTractInfo(TractInfo):
    """Tract for the HealpixSkyMap"""
    def __init__(self, nSide, ident, nest, patchInnerDimensions, patchBorder, ctrCoord, tractOverlap, wcs):
        """Set vertices from nside, ident, nest"""
        theta, phi = healpy.vec2ang(numpy.transpose(healpy.boundary(nSide, ident, nest=nest)))
        vertexList = [angToCoord(thetaphi) for thetaphi in zip(theta,phi)]
        super(HealpixTractInfo, self).__init__(ident, patchInnerDimensions, patchBorder, ctrCoord,
                                               vertexList, tractOverlap, wcs)


class HealpixSkyMapConfig(BaseSkyMap.ConfigClass):
    """Configuration for the HealpixSkyMap"""
    nSide = Field(dtype=int, default=0, doc="Number of sides, expressed in powers of 2")
    nest = Field(dtype=bool, default=False, doc="Use NEST ordering instead of RING?")
    def setDefaults(self):
        self.rotation = 45 # HEALPixels are oriented at 45 degrees

class HealpixSkyMap(BaseSkyMap):
    """HEALPix-based sky map pixelization.

    We put a Tract at the position of each HEALPixel.
    """
    ConfigClass = HealpixSkyMapConfig
    _version = (1, 0) # for pickle
    numAngles = 4 # Number of angles for vertices

    def __init__(self, config=None, version=0):
        """Constructor

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        @param[in] version: software version of this class, to retain compatibility with old instances
        """
        super(HealpixSkyMap, self).__init__(config)
        self._version = version
        self.nside = 1 << self.config.nSide
        self._numTracts = healpy.nside2npix(self._nside)
        self._tractCache = {}
        self._tractInfo = None # We shouldn't be using this; we will generate tracts on demand

    def __reduce__(self):
        """To support pickling"""
        return (self.__class__, (self.config, self._version))

    def findTract(self, coord):
        """Find the tract whose inner region includes the coord."""
        theta, phi = coordToAng(coord.toIcrs())
        index = healpy.ang2pix(self.nside, theta, phi, nest=self.config.nest)
        return self[index]

    def __getitem__(self, index):
        """Get the TractInfo for a particular index

        The tract is returned from a cache, if available, otherwise generated
        on the fly.
        """
        if index < 0 or index > self._numTracts:
            raise IndexError("Index out of range: %d vs %d" % (index, self._numTracts))
        if index in self._tractCache:
            return self._tractCache[index]
        center = angToCoord(healpy.pix2ang(self.nside, index, nest=self.config.nest))
        wcs = self._wcsFactory.makeWcs(crPixPos=afwGeom.Point2D(0,0), crValCoord=center)
        tract = HealpixTractInfo(self.nside, index, self.config.nest, self.config.patchInnerDimensions,
                                 self.config.patchBorder, center, self.config.tractOverlap*afwGeom.degrees,
                                 wcs)
        self._tractCache[index] = tract
        return tract

    def __iter__(self):
        """Iterator over tracts"""
        for i in xrange(self._numTracts):
            yield self[i]

    def __len__(self):
        """Length is number of tracts"""
        return self._numTracts


