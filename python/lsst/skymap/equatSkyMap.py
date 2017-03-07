from __future__ import division
from builtins import range
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
import lsst.pex.config as pexConfig
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
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

    EquatSkyMap represents an equatorial band of sky divided along declination into overlapping tracts.
    """
    ConfigClass = EquatSkyMapConfig
    _version = (1, 0) # for pickle

    def __init__(self, config=None):
        """Construct a EquatSkyMap

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        """
        BaseSkyMap.__init__(self, config)

        decRange = tuple(afwGeom.Angle(dr, afwGeom.degrees) for dr in self.config.decRange)
        midDec = (decRange[0] + decRange[1]) / 2.0
        tractWidthRA = afwGeom.Angle(360.0 / self.config.numTracts, afwGeom.degrees)
        tractOverlap = afwGeom.Angle(self.config.tractOverlap, afwGeom.degrees)

        for id in range(self.config.numTracts):
            begRA = tractWidthRA * id
            endRA = begRA + tractWidthRA
            vertexCoordList = (
                afwCoord.IcrsCoord(begRA, decRange[0]),
                afwCoord.IcrsCoord(endRA, decRange[0]),
                afwCoord.IcrsCoord(endRA, decRange[1]),
                afwCoord.IcrsCoord(begRA, decRange[1]),
            )

            midRA = begRA + tractWidthRA / 2.0
            ctrCoord = afwCoord.IcrsCoord(midRA, midDec)

            # CRVal must have Dec=0 for symmetry about the equator
            crValCoord = afwCoord.IcrsCoord(midRA, afwGeom.Angle(0.0))

            # make initial WCS; don't worry about crPixPos because TractInfo will shift it as required
            wcs = self._wcsFactory.makeWcs(crPixPos=afwGeom.Point2D(0, 0), crValCoord=crValCoord)

            self._tractInfoList.append(TractInfo(
                id=id,
                patchInnerDimensions=self.config.patchInnerDimensions,
                patchBorder=self.config.patchBorder,
                ctrCoord=ctrCoord,
                vertexCoordList=vertexCoordList,
                tractOverlap=tractOverlap,
                wcs=wcs,
            ))

    def __getstate__(self):
        """Support pickle

        @return a dict containing:
        - version: a pair of ints
        - config: the config
        """
        return dict(
            version=self._version,
            config=self.config,
        )

    def __setstate__(self, stateDict):
        """Support unpickle

        @param[in] stateDict: a dict containing:
        - version: a pair of ints
        - config: the config
        """
        version = stateDict["version"]
        if version >= (2, 0):
            raise runtimeError("Version = %s >= (2,0); cannot unpickle" % (version,))
        self.__init__(stateDict["config"])

    def getVersion(self):
        """Return version (e.g. for pickle)

        @return version as a pair of integers
        """
        return self._version
