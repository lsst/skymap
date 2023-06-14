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

__all__ = ['EquatSkyMapConfig', 'EquatSkyMap']

import struct

import lsst.pex.config as pexConfig
import lsst.geom as geom
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo


class EquatSkyMapConfig(BaseSkyMap.ConfigClass):
    numTracts = pexConfig.Field(
        doc="number of tracts; warning: TAN projection requires at least 3",
        dtype=int,
        default=4,
    )
    decRange = pexConfig.ListField(
        doc="range of declination (deg)",
        dtype=float,
        length=2,
        default=(-1.25, 1.25),
    )

    def setDefaults(self):
        self.projection = "CEA"


class EquatSkyMap(BaseSkyMap):
    """Equatorial sky map pixelization, e.g. for SDSS stripe 82 image data.

    EquatSkyMap represents an equatorial band of sky divided along declination
    into overlapping tracts.

    Parameters
    ----------
    config : `lsst.skymap.BaseSkyMapConfig` (optional)
        The configuration for this SkyMap; if None use the default config.
    """
    ConfigClass = EquatSkyMapConfig
    _version = (1, 0)  # for pickle

    def __init__(self, config=None):
        BaseSkyMap.__init__(self, config)

        decRange = tuple(geom.Angle(dr, geom.degrees) for dr in self.config.decRange)
        midDec = (decRange[0] + decRange[1]) / 2.0
        tractWidthRA = geom.Angle(360.0 / self.config.numTracts, geom.degrees)
        tractOverlap = geom.Angle(self.config.tractOverlap, geom.degrees)

        for id in range(self.config.numTracts):
            begRA = tractWidthRA * id
            endRA = begRA + tractWidthRA
            vertexCoordList = (
                geom.SpherePoint(begRA, decRange[0]),
                geom.SpherePoint(endRA, decRange[0]),
                geom.SpherePoint(endRA, decRange[1]),
                geom.SpherePoint(begRA, decRange[1]),
            )

            midRA = begRA + tractWidthRA / 2.0
            ctrCoord = geom.SpherePoint(midRA, midDec)

            # CRVal must have Dec=0 for symmetry about the equator
            crValCoord = geom.SpherePoint(midRA, geom.Angle(0.0))

            # Make initial WCS; don't worry about crPixPos because TractInfo
            # will shift it as required.
            wcs = self._wcsFactory.makeWcs(crPixPos=geom.Point2D(0, 0), crValCoord=crValCoord)

            self._tractInfoList.append(TractInfo(
                id=id,
                tractBuilder=self._tractBuilder,
                ctrCoord=ctrCoord,
                vertexCoordList=vertexCoordList,
                tractOverlap=tractOverlap,
                wcs=wcs,
            ))

    def __getstate__(self):
        """Support pickle.

        Returns
        -------
        stateDict : `dict`
            a dict containing:
            - version: a pair of ints
            - config: the config
        """
        return dict(
            version=self._version,
            config=self.config,
        )

    def __setstate__(self, stateDict):
        """Support unpickle

        Parameters
        ----------
        stateDict : `dict`
            a dict containing:
            - version: a pair of ints
            - config: the config
        """
        version = stateDict["version"]
        if version >= (2, 0):
            raise RuntimeError("Version = %s >= (2,0); cannot unpickle" % (version,))
        self.__init__(stateDict["config"])

    def getVersion(self):
        """Return version (e.g. for pickle).

        Returns
        -------
        result : `tuple` of `int`
            Version as a pair of integers.
        """
        return self._version

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1."""
        sha1.update(struct.pack("<i2d", self.config.numTracts, *self.config.decRange))
