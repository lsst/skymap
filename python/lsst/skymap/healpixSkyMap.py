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
        """A null pattern which blows up when we try to read it"""
        def __getattr__(self, name):
            raise RuntimeError("Was unable to import healpy: %s" % e)
    healpy = DummyHealpy()

from lsst.pex.config import Field
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo

class HealpixSkyMapConfig(BaseSkyMap.ConfigClass):
    nSide = Field(dtype=int, default=0, doc="Number of sides, expressed in powers of 2")
    nest = Field(dtype=bool, default=False, doc="Use NEST ordering instead of RING?")

class HealpixSkyMap(BaseSkyMap):
    """HEALPix-based sky map pixelization.

    We put a Tract at the position of each HEALPixel.
    """
    ConfigClass = HealpixSkyMapConfig
    _version = (1, 0) # for pickle
    numAngles = 12 # Number of angles for vertices

    def __init__(self, config=None):
        """Constructor

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        """
        super(HealpixSkyMap, self).__init__(config)
        self.nside = 1 << self.config.nSide
        numPixels = healpy.nside2npix(self.nside)
        maxRadius = healpy.max_pixrad(self.nside) * afwGeom.radians
        indices = numpy.arange(numPixels)
        theta, phi = healpy.pix2ang(self.nside, indices, nest=self.config.nest)
        for i, th, ph in zip(indices, theta, phi):
            center = IcrsCoord(afwGeom.Angle(ph, afwGeom.radians),
                               afwGeom.Angle(th - 0.5*numpy.pi, afwGeom.radians))
            vertices = []
            for angle in numpy.linspace(0.0, 360.0, self.numAngles, endpoint=False):
                coord = center.clone()
                coord.offset(afwGeom.Angle(angle, afwGeom.degrees), maxRadius)
                vertices.append(coord)

            # make initial WCS; don't worry about crPixPos because TractInfo will shift it as required
            wcs = self._wcsFactory.makeWcs(crPixPos=afwGeom.Point2D(0,0), crValCoord=center)

            tract = TractInfo(id=i, patchInnerDimensions=self.config.patchInnerDimensions,
                              patchBorder=self.config.patchBorder, ctrCoord=center, vertexCoordList=vertices,
                              tractOverlap=self.config.tractOverlap*afwGeom.degrees, wcs=wcs)
            self._tractInfoList.append(tract)


    def __getstate__(self):
        """Return state, for pickling"""
        return dict(version=self._version, config=self.config)

    def __setstate__(self, state):
        """Set state, from pickling"""
        version = state["version"]
        if version >= (2, 0):
            raise RuntimeError("Version = %s >= (2,0); cannot unpickle" % (version,))
        self.__init__(state["config"])
    
    def findTract(self, coord):
        """Find the tract whose inner region includes the coord."""
        coord = coord.toIcrs()
        return self[healpy.ang2pix(self.nside, coord.getLatitude() + 0.5*numpy.pi,
                                   coord.getLongitude(), nest=self.config.nest)]
