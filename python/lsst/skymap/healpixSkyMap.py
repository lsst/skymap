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

__all__ = ['HealpixSkyMapConfig', 'HealpixSkyMap']

from deprecated.sphinx import deprecated
import struct
import numpy

import hpgeom

from lsst.pex.config import Field
import lsst.geom as geom
from .cachingSkyMap import CachingSkyMap
from .tractInfo import TractInfo


# TODO: Remove with DM-44799
@deprecated("angToCoord has been deprecated and will be removed after v28.",
            category=FutureWarning, version=28)
def angToCoord(thetaphi):
    """Convert hpgeom ang to an lsst.geom.SpherePoint

    The angle is provided as a single object, thetaphi, so the output
    of hpgeom functions can be directed to this function without
    additional translation.
    """
    return geom.SpherePoint(float(thetaphi[1]), float(thetaphi[0] - 0.5*numpy.pi), geom.radians)


# TODO: Remove with DM-44799
@deprecated("coordToAnge has been deprecated and will be removed after v28.",
            category=FutureWarning, version=28)
def coordToAng(coord):
    """Convert an lsst.geom.SpherePoint to a hpgeom ang (theta, phi)

    The Healpix convention is that 0 <= theta <= pi, 0 <= phi < 2pi.
    """
    return (coord.getLatitude().asRadians() + 0.5*numpy.pi, coord.getLongitude().asRadians())


# TODO: Remove with DM-44799
@deprecated("HealpixTractInfo has been deprecated and will be removed after v28.",
            category=FutureWarning, version=28)
class HealpixTractInfo(TractInfo):
    """Tract for the HealpixSkyMap"""

    def __init__(self, nSide, ident, nest, tractBuilder, ctrCoord, tractOverlap, wcs):
        """Set vertices from nside, ident, nest"""
        theta, phi = hpgeom.boundaries(nSide, ident, nest=nest, lonlat=False)
        vertexList = [angToCoord(thetaphi) for thetaphi in zip(theta, phi)]
        super(HealpixTractInfo, self).__init__(ident, tractBuilder, ctrCoord,
                                               vertexList, tractOverlap, wcs)


# TODO: Remove with DM-44799
@deprecated("HealpixSkyMapConfig has been deprecated and will be removed after v28.",
            category=FutureWarning, version=28)
class HealpixSkyMapConfig(CachingSkyMap.ConfigClass):
    """Configuration for the HealpixSkyMap"""
    log2NSide = Field(dtype=int, default=0, doc="Number of sides, expressed in powers of 2")
    nest = Field(dtype=bool, default=False, doc="Use NEST ordering instead of RING?")

    def setDefaults(self):
        self.rotation = 45  # HEALPixels are oriented at 45 degrees


# TODO: Remove with DM-44799
@deprecated("HealpixSkyMap has been deprecated and will be removed after v28.",
            category=FutureWarning, version=28)
class HealpixSkyMap(CachingSkyMap):
    """HEALPix-based sky map pixelization.

    We put a Tract at the position of each HEALPixel.


    Parameters
    ----------
    config : `lsst.skymap.BaseSkyMapConfig`
        The configuration for this SkyMap.
    version : `int` or `tuple` of `int` (optional)
        Software version of this class, to retain compatibility with old
        instances.
    """
    ConfigClass = HealpixSkyMapConfig
    _version = (1, 0)  # for pickle
    numAngles = 4  # Number of angles for vertices

    def __init__(self, config, version=0):
        self._nside = 1 << config.log2NSide
        numTracts = hpgeom.nside_to_npixel(self._nside)
        super(HealpixSkyMap, self).__init__(numTracts, config, version)

    def findTract(self, coord):
        """Find the tract whose inner region includes the coord.

        Parameters
        ----------
        coord : `lsst.geom.SpherePoint`
            ICRS sky coordinate to search for.

        Returns
        -------
        tractInfo : `TractInfo`
            Info for tract whose inner region includes the coord.
        """
        theta, phi = coordToAng(coord)
        index = hpgeom.angle_to_pixel(self._nside, theta, phi, nest=self.config.nest, lonlat=False)
        return self[index]

    def generateTract(self, index):
        """Generate TractInfo for the specified tract index."""
        center = angToCoord(hpgeom.pixel_to_angle(self._nside, index, nest=self.config.nest, lonlat=False))
        wcs = self._wcsFactory.makeWcs(crPixPos=geom.Point2D(0, 0), crValCoord=center)
        return HealpixTractInfo(self._nside, index, self.config.nest, self._tractBuilder,
                                center, self.config.tractOverlap*geom.degrees,
                                wcs)

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1."""
        sha1.update(struct.pack("<i?", self.config.log2NSide, self.config.nest))
