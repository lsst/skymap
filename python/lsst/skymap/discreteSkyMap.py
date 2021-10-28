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

__all__ = ["DiscreteSkyMap"]

import struct

from lsst.pex.config import ListField
import lsst.geom as geom
from .cachingSkyMap import CachingSkyMap
from .tractInfo import ExplicitTractInfo


class DiscreteSkyMapConfig(CachingSkyMap.ConfigClass):
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


class DiscreteSkyMap(CachingSkyMap):
    """Discrete sky map pixelization.

    We put a square Tract at each of the nominated coordinates.

    Parameters
    ----------
    config : `lsst.skyMap.BaseSkyMapConfig`
        The configuration for this SkyMap; if None use the default config.
    version : `int` or `tuple` of `int` (optional)
        Software version of this class, to retain compatibility with old
        instances.
    """
    ConfigClass = DiscreteSkyMapConfig
    _version = (1, 0)  # for pickle

    def __init__(self, config, version=0):
        numTracts = len(config.radiusList)
        super(DiscreteSkyMap, self).__init__(numTracts, config, version)

    def generateTract(self, index):
        """Generate TractInfo for the specified tract index."""
        center = geom.SpherePoint(self.config.raList[index], self.config.decList[index], geom.degrees)
        radius = self.config.radiusList[index]
        wcs = self._wcsFactory.makeWcs(crPixPos=geom.Point2D(0, 0), crValCoord=center)
        return ExplicitTractInfo(index, self._tractBuilder, center,
                                 radius*geom.degrees, self.config.tractOverlap*geom.degrees, wcs)

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1."""
        fmt = "<{}d".format(len(self.config.radiusList))
        sha1.update(struct.pack(fmt, *self.config.raList))
        sha1.update(struct.pack(fmt, *self.config.decList))
        sha1.update(struct.pack(fmt, *self.config.radiusList))
