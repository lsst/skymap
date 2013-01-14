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

import math

from lsst.pex.config import Field
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from .cachingSkyMap import CachingSkyMap
from .tractInfo import ExplicitTractInfo

__all__ = ["RingsSkyMap"]

class RingsSkyMapConfig(CachingSkyMap.ConfigClass):
    """Configuration for the RingsSkyMap"""
    numRings = Field(dtype=int, doc="Number of rings", check=lambda x: x > 0)
    raStart = Field(dtype=float, default=0.0, doc="Starting center RA for each ring (degrees)",
                    check=lambda x: x >= 0.0 and x < 360.0)


class RingsSkyMap(CachingSkyMap):
    """Rings sky map pixelization.

    We divide the sphere into N rings of Declination, plus the two polar
    caps, which sets the size of the individual tracts.  The rings are
    divided in RA into an integral number of tracts of this size; this
    division is made at the Declination closest to zero so as to ensure
    full overlap.
    """
    ConfigClass = RingsSkyMapConfig
    _version = (1, 0) # for pickle

    def __init__(self, config, version=0):
        """Constructor

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        @param[in] version: software version of this class, to retain compatibility with old instances
        """
        # We count rings from south to north
        # Note: pole caps together count for one additional ring
        self._ringSize = math.pi / (config.numRings + 1) # Size of a ring in Declination (radians)
        self._ringNums = [] # Number of tracts for each ring
        for i in range(config.numRings):
            startDec = self._ringSize*(i + 0.5) - 0.5*math.pi
            stopDec = startDec + self._ringSize
            dec = min(math.fabs(startDec), math.fabs(stopDec)) # Declination for determining division in RA
            self._ringNums.append(int(2*math.pi*math.cos(dec)/self._ringSize) + 1)
        numTracts = sum(self._ringNums) + 2
        super(RingsSkyMap, self).__init__(numTracts, config, version)

    def getRingIndices(self, index):
        """Calculate ring indices given a numerical index of a tract

        The ring indices are the ring number and the tract number within
        the ring.

        The ring number is -1 for the south polar cap and increases to the
        north.  The north polar cap has ring number = numRings.  The tract
        number is zero for either of the polar caps.
        """
        if index == 0: # South polar cap
            return -1, 0
        if index == self._numTracts - 1: # North polar cap
            return self.config.numRings, 0
        index -= 1
        ring = 0
        while ring < self.config.numRings and index > self._ringNums[ring]:
            index -= self._ringNums[ring]
            ring += 1
        return ring, index

    def generateTract(self, index):
        """Generate the TractInfo for this index"""
        ringNum, tractNum = self.getRingIndices(index)
        if ringNum == -1: # South polar cap
            ra, dec = 0, -0.5*math.pi
        elif ringNum == self.config.numRings: # North polar cap
            ra, dec = 0, 0.5*math.pi
        else:
            dec = self._ringSize*(ringNum + 1) - 0.5*math.pi
            ra = math.fmod(self.config.raStart + 2*math.pi*tractNum/self._ringNums[ringNum], 2*math.pi)

        center = IcrsCoord(ra*afwGeom.radians, dec*afwGeom.radians)
        wcs = self._wcsFactory.makeWcs(crPixPos=afwGeom.Point2D(0,0), crValCoord=center)
        return ExplicitTractInfo(index, self.config.patchInnerDimensions, self.config.patchBorder, center,
                                 0.5*self._ringSize*afwGeom.radians, self.config.tractOverlap*afwGeom.degrees,
                                 wcs)
