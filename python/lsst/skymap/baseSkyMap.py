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
"""
@todo
- Consider tweaking pixel scale so the average scale is as specified, rather than the scale at the center
"""
import hashlib
import struct

from builtins import object
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
from . import detail

__all__ = ["BaseSkyMap"]


class BaseSkyMapConfig(pexConfig.Config):
    patchInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of patches (x,y pixels)",
        dtype=int,
        length=2,
        default=(4000, 4000),
    )
    patchBorder = pexConfig.Field(
        doc="border between patch inner and outer bbox (pixels)",
        dtype=int,
        default=100,
    )
    tractOverlap = pexConfig.Field(
        doc="minimum overlap between adjacent sky tracts, on the sky (deg)",
        dtype=float,
        default=1.0,
    )
    pixelScale = pexConfig.Field(
        doc="nominal pixel scale (arcsec/pixel)",
        dtype=float,
        default=0.333
    )
    projection = pexConfig.Field(
        doc="""one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        """,
        dtype=str,
        default="STG",
    )
    rotation = pexConfig.Field(
        doc="Rotation for WCS (deg)",
        dtype=float,
        default=0,
    )


class BaseSkyMap(object):
    """A collection of overlapping Tracts that map part or all of the sky.

    See TractInfo for more information.

    BaseSkyMap is an abstract base class. Subclasses must do the following:
    @li define __init__ and have it construct the TractInfo objects and put them in _tractInfoList
    @li define __getstate__ and __setstate__ to allow pickling (the butler saves sky maps using pickle);
        see DodecaSkyMap for an example of how to do this. (Most of that code could be moved
        into this base class, but that would make it harder to handle older versions of pickle data.)
    @li define updateSha1 to add any subclass-specific state to the hash.
    """
    ConfigClass = BaseSkyMapConfig

    def __init__(self, config=None):
        """Construct a BaseSkyMap

        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        """
        if config is None:
            config = self.ConfigClass()
        config.freeze() # just to be sure, e.g. for pickling
        self.config = config
        self._tractInfoList = []
        self._wcsFactory = detail.WcsFactory(
            pixelScale=afwGeom.Angle(self.config.pixelScale, afwGeom.arcseconds),
            projection=self.config.projection,
            rotation=afwGeom.Angle(self.config.rotation, afwGeom.degrees),
        )
        self._sha1 = None

    def findTract(self, coord):
        """Find the tract whose center is nearest the specified coord.

        @param[in] coord: sky coordinate (afwCoord.Coord)
        @return TractInfo of tract whose center is nearest the specified coord

        @warning:
        - if tracts do not cover the whole sky then the returned tract may not include the coord

        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between multiple sky tract centers then one is arbitrarily chosen.
        - The default implementation is not very efficient; subclasses may wish to override.
        """
        icrsCoord = coord.toIcrs()
        distTractInfoList = []
        for tractInfo in self:
            angSep = icrsCoord.angularSeparation(tractInfo.getCtrCoord()).asDegrees()
            distTractInfoList.append((angSep, tractInfo))
        distTractInfoList.sort()
        return distTractInfoList[0][1]

    def findTractPatchList(self, coordList):
        """Find tracts and patches that overlap a region

        @param[in] coordList: list of sky coordinates (afwCoord.Coord)
        @return list of (TractInfo, list of PatchInfo) for tracts and patches that contain,
            or may contain, the specified region. The list will be empty if there is no overlap.

        @warning this uses a naive algorithm that may find some tracts and patches that do not overlap
            the region (especially if the region is not a rectangle aligned along patch x,y).
        """
        retList = []
        for tractInfo in self:
            patchList = tractInfo.findPatchList(coordList)
            if patchList:
                retList.append((tractInfo, patchList))
        return retList

    def findClosestTractPatchList(self, coordList):
        """Find closest tract and patches that overlap coordinates

        @param[in] coordList: list of sky coordinates (afwCoord.Coord)
        @return list of (TractInfo, list of PatchInfo) for tracts and patches that contain,
            or may contain, the specified region. The list will be empty if there is no overlap.
        """
        retList = []
        for coord in coordList:
            tractInfo = self.findTract(coord)
            patchList = tractInfo.findPatchList(coordList)
            if patchList and not (tractInfo, patchList) in retList:
                retList.append((tractInfo, patchList))
        return retList

    def __getitem__(self, ind):
        return self._tractInfoList[ind]

    def __iter__(self):
        return iter(self._tractInfoList)

    def __len__(self):
        return len(self._tractInfoList)

    def getSha1(self):
        """Return a SHA1 hash that uniquely identifies this SkyMap instance.

        Returns
        -------
        sha1 : bytes
            A 20-byte hash that uniquely identifies this SkyMap instance.

        Subclasses should almost always override `updateSha1()` instead of
        this function to add subclass-specific state to the hash.
        """
        if self._sha1 is None:
            sha1 = hashlib.sha1()
            sha1.update(type(self).__name__.encode('utf-8'))
            configPacked = struct.pack(
                "<iiidd3sd",
                self.config.patchInnerDimensions[0],
                self.config.patchInnerDimensions[1],
                self.config.patchBorder,
                self.config.tractOverlap,
                self.config.pixelScale,
                self.config.projection.encode('ascii'),
                self.config.rotation
            )
            sha1.update(configPacked)
            self.updateSha1(sha1)
            self._sha1 = sha1.digest()
        return self._sha1

    def updateSha1(self, sha1):
        """Add subclass-specific state or configuration options to the SHA1.

        Parameters
        ----------
        sha1 : hashlib.sha1
            A hashlib object on which `update()` can be called to add
            additional state to the hash.

        This method is conceptually "protected": it should be reimplemented by
        all subclasses, but called only by the base class implementation of
        `getSha1()`.
        """
        raise NotImplementedError()
