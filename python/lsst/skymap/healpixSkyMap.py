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

import struct
import numpy

# We want to register the HealpixSkyMap, but want "healpy" to be an
# optional dependency.  However, the HealpixSkyMap requires the use
# of healpy.  Therefore, we'll only raise an exception on the healpy
# import when it comes time to using it.
try:
    import healpy
except Exception as e:
    class DummyHealpy:
        """An object which blows up when we try to read it"""

        def __getattr__(self, name, e=e):
            raise RuntimeError("Was unable to import healpy: %s" % e)
    healpy = DummyHealpy()

from lsst.pex.config import Field
import lsst.geom as geom
from .cachingSkyMap import CachingSkyMap
from .tractInfo import TractInfo


def angToCoord(thetaphi):
    """Convert healpy's ang to an lsst.geom.SpherePoint

    The ang is provided as a single object, thetaphi, so the output
    of healpy functions can be directed to this function without
    additional translation.
    """
    return geom.SpherePoint(float(thetaphi[1]), float(thetaphi[0] - 0.5*numpy.pi), geom.radians)


def coordToAng(coord):
    """Convert an lsst.geom.SpherePoint to a healpy ang (theta, phi)

    The Healpix convention is that 0 <= theta <= pi, 0 <= phi < 2pi.
    """
    return (coord.getLatitude().asRadians() + 0.5*numpy.pi, coord.getLongitude().asRadians())


class HealpixTractInfo(TractInfo):
    """Tract for the HealpixSkyMap"""

    def __init__(self, nSide, ident, nest, patchBuilder, ctrCoord, tractOverlap, wcs):
        """Set vertices from nside, ident, nest"""
        theta, phi = healpy.vec2ang(numpy.transpose(healpy.boundaries(nSide, ident, nest=nest)))
        vertexList = [angToCoord(thetaphi) for thetaphi in zip(theta, phi)]
        super(HealpixTractInfo, self).__init__(ident, patchBuilder, ctrCoord,
                                               vertexList, tractOverlap, wcs)


class HealpixSkyMapConfig(CachingSkyMap.ConfigClass):
    """Configuration for the HealpixSkyMap"""
    log2NSide = Field(dtype=int, default=0, doc="Number of sides, expressed in powers of 2")
    nest = Field(dtype=bool, default=False, doc="Use NEST ordering instead of RING?")

    def setDefaults(self):
        self.rotation = 45  # HEALPixels are oriented at 45 degrees


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
        numTracts = healpy.nside2npix(self._nside)
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
        index = healpy.ang2pix(self._nside, theta, phi, nest=self.config.nest)
        return self[index]

    def generateTract(self, index):
        """Generate TractInfo for the specified tract index."""
        center = angToCoord(healpy.pix2ang(self._nside, index, nest=self.config.nest))
        wcs = self._wcsFactory.makeWcs(crPixPos=geom.Point2D(0, 0), crValCoord=center)
        return HealpixTractInfo(self._nside, index, self.config.nest, self._patchBuilder,
                                center, self.config.tractOverlap*geom.degrees,
                                wcs)

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1."""
        sha1.update(struct.pack("<i?", self.config.log2NSide, self.config.nest))
