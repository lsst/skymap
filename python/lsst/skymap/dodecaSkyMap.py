#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/) .
#
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

"""Consider tweaking pixel scale so the average scale is as specified, rather
than the scale at the center.
"""

__all__ = ['DodecaSkyMapConfig', 'DodecaSkyMap']

import struct

import lsst.pex.config as pexConfig
import lsst.geom as geom
from . import detail
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo


class DodecaSkyMapConfig(BaseSkyMap.ConfigClass):
    withTractsOnPoles = pexConfig.Field(
        doc="if True center a tract on each pole, else put a vertex on each pole",
        dtype=bool,
        default=False,
    )

    def setDefaults(self):
        self.tractOverlap = 3.5
        self.patchBorder = 250
        self.pixelScale = 10.0 / 50.0  # LSST plate scale is 50 um/arcsec and pixel size is 10 um
        self.patchInnerDimensions = (4000, 4000)
        self.projection = "STG"


class DodecaSkyMap(BaseSkyMap):
    """Dodecahedron-based sky map pixelization.

    DodecaSkyMap divides the sky into 12 overlapping Tracts arranged as the
    faces of a dodecahedron.

    Parameters
    ----------
    config : `lsst.skymap.BaseSkyMapConfig` (optional)
        The configuration for this SkyMap; if None use the default config.
    """
    ConfigClass = DodecaSkyMapConfig
    _version = (1, 0)  # for pickle

    def __init__(self, config=None):
        BaseSkyMap.__init__(self, config)
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles=self.config.withTractsOnPoles)

        tractOverlap = geom.Angle(self.config.tractOverlap, geom.degrees)

        for id in range(12):
            tractVec = self._dodecahedron.getFaceCtr(id)
            tractCoord = detail.coordFromVec(tractVec, defRA=geom.Angle(0))
            tractRA = tractCoord.getLongitude()
            vertexVecList = self._dodecahedron.getVertices(id)

            # make initial WCS; don't worry about crPixPos because TractInfo will shift it as required
            wcs = self._wcsFactory.makeWcs(crPixPos=geom.Point2D(0, 0), crValCoord=tractCoord)

            self._tractInfoList.append(
                TractInfo(
                    id=id,
                    patchInnerDimensions=self.config.patchInnerDimensions,
                    patchBorder=self.config.patchBorder,
                    ctrCoord=tractCoord,
                    vertexCoordList=[detail.coordFromVec(vec, defRA=tractRA) for vec in vertexVecList],
                    tractOverlap=tractOverlap,
                    wcs=wcs,
                )
            )

    def __getstate__(self):
        """Support pickle.

        Returns
        -------
        result : `dict`
            A dict containing:
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
            - version: a pair of ints
            - config: the config
        """
        version = stateDict["version"]
        if version >= (2, 0):
            raise RuntimeError("Version = %s >= (2,0); cannot unpickle" % (version,))
        self.__init__(stateDict["config"])

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
        return self[self._dodecahedron.getFaceInd(coord.getVector())]

    def getVersion(self):
        """Return version (e.g. for pickle).

        Returns
        -------
        version : `tuple` of `int`
            Version as a pair of integers.
        """
        return self._version

    def getWithTractsOnPoles(self):
        """Return True if there are tracts centered on the poles.
        """
        return self._dodecahedron.getWithFacesOnPoles()

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1."""
        sha1.update(struct.pack("<?", self.config.withTractsOnPoles))
