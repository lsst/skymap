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
todo: Consider tweaking pixel scale so the average scale is as specified,
rather than the scale at the center.
"""

__all__ = ["BaseSkyMapConfig", "BaseSkyMap"]

import hashlib
import struct

import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
from lsst.geom import SpherePoint, Angle, arcseconds, degrees
from . import detail


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
        doc="one of the FITS WCS projection codes, such as:"
        "- STG: stereographic projection"
        "- MOL: Molleweide's projection"
        "- TAN: tangent-plane projection",
        dtype=str,
        default="STG",
    )
    rotation = pexConfig.Field(
        doc="Rotation for WCS (deg)",
        dtype=float,
        default=0,
    )


class BaseSkyMap:
    """A collection of overlapping Tracts that map part or all of the sky.

    See TractInfo for more information.

    Parameters
    ----------
    config : `BaseSkyMapConfig` or None (optional)
        The configuration for this SkyMap; if None use the default config.

    Notes
    -----
    BaseSkyMap is an abstract base class. Subclasses must do the following:
    define ``__init__`` and have it construct the TractInfo objects and put
    them in ``__tractInfoList__`` define ``__getstate__`` and ``__setstate__``
    to allow pickling (the butler saves sky maps using pickle);
    see DodecaSkyMap for an example of how to do this. (Most of that code could
    be moved into this base class, but that would make it harder to handle
    older versions of pickle data.) define updateSha1 to add any
    subclass-specific state to the hash.

    All SkyMap subclasses must be conceptually immutable; they must always
    refer to the same set of mathematical tracts and patches even if the in-
    memory representation of those objects changes.
    """

    ConfigClass = BaseSkyMapConfig

    def __init__(self, config=None):
        if config is None:
            config = self.ConfigClass()
        config.freeze()  # just to be sure, e.g. for pickling
        self.config = config
        self._tractInfoList = []
        self._wcsFactory = detail.WcsFactory(
            pixelScale=Angle(self.config.pixelScale, arcseconds),
            projection=self.config.projection,
            rotation=Angle(self.config.rotation, degrees),
        )
        self._sha1 = None

    def findTract(self, coord):
        """Find the tract whose center is nearest the specified coord.

        Parameters
        ----------
        coord : `lsst.geom.SpherePoint`
            ICRS sky coordinate to search for.

        Returns
        -------
        result : `TractInfo`
            TractInfo of tract whose center is nearest the specified coord.

        Notes
        -----
        - If coord is equidistant between multiple sky tract centers then one
          is arbitrarily chosen.

        - The default implementation is not very efficient; subclasses may wish
          to override.

        **Warning:**
        If tracts do not cover the whole sky then the returned tract may not
        include the coord.
        """
        distTractInfoList = []
        for i, tractInfo in enumerate(self):
            angSep = coord.separation(tractInfo.getCtrCoord()).asDegrees()
            # include index in order to disambiguate identical angSep values
            distTractInfoList.append((angSep, i, tractInfo))
        distTractInfoList.sort()
        return distTractInfoList[0][2]

    def findTractPatchList(self, coordList):
        """Find tracts and patches that overlap a region.

        Parameters
        ----------
        coordList : `list` of `lsst.geom.SpherePoint`
            List of ICRS sky coordinates to search for.

        Returns
        -------
        reList : `list` of (`TractInfo`, `list` of `PatchInfo`)
            For tracts and patches that contain, or may contain, the specified
            region. The list will be empty if there is no overlap.

        Notes
        -----
        **warning:**
            This uses a naive algorithm that may find some tracts and patches
            that do not overlap the region (especially if the region is not a
            rectangle aligned along patch x, y).
        """
        retList = []
        for tractInfo in self:
            patchList = tractInfo.findPatchList(coordList)
            if patchList:
                retList.append((tractInfo, patchList))
        return retList

    def findClosestTractPatchList(self, coordList):
        """Find closest tract and patches that overlap coordinates.

        Parameters
        ----------
        coordList : `lsst.geom.SpherePoint`
            List of ICRS sky coordinates to search for.

        Returns
        -------
        retList : `list`
            list of (TractInfo, list of PatchInfo) for tracts and patches
            that contain, or may contain, the specified region.
            The list will be empty if there is no overlap.
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

    def __hash__(self):
        return hash(self.getSha1())

    def __eq__(self, other):
        try:
            return self.getSha1() == other.getSha1()
        except AttributeError:
            return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    def logSkyMapInfo(self, log):
        """Write information about a sky map to supplied log

        Parameters
        ----------
        log : `lsst.log.Log`
            Log object that information about skymap will be written
        """
        log.info("sky map has %s tracts" % (len(self),))
        for tractInfo in self:
            wcs = tractInfo.getWcs()
            posBox = afwGeom.Box2D(tractInfo.getBBox())
            pixelPosList = (
                posBox.getMin(),
                afwGeom.Point2D(posBox.getMaxX(), posBox.getMinY()),
                posBox.getMax(),
                afwGeom.Point2D(posBox.getMinX(), posBox.getMaxY()),
            )
            skyPosList = [wcs.pixelToSky(pos).getPosition(afwGeom.degrees) for pos in pixelPosList]
            posStrList = ["(%0.3f, %0.3f)" % tuple(skyPos) for skyPos in skyPosList]
            log.info("tract %s has corners %s (RA, Dec deg) and %s x %s patches" %
                     (tractInfo.getId(), ", ".join(posStrList),
                      tractInfo.getNumPatches()[0], tractInfo.getNumPatches()[1]))

    def getSha1(self):
        """Return a SHA1 hash that uniquely identifies this SkyMap instance.

        Returns
        -------
        sha1 : `bytes`
            A 20-byte hash that uniquely identifies this SkyMap instance.

        Notes
        -----
        Subclasses should almost always override ``updateSha1`` instead of
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
        sha1 : `hashlib.sha1`
            A hashlib object on which `update` can be called to add
            additional state to the hash.

        Notes
        -----
        This method is conceptually "protected" : it should be reimplemented by
        all subclasses, but called only by the base class implementation of
        `getSha1` .
        """
        raise NotImplementedError()

    def register(self, name, registry):
        """Add SkyMap, Tract, and Patch Dimension entries to the given Gen3
        Butler Registry.

        Parameters
        ----------
        name : `str`
            The name of the skymap.
        registry : `lsst.daf.butler.Registry`
            The registry to add to.
        """
        nxMax = 0
        nyMax = 0
        for tractInfo in self:
            nx, ny = tractInfo.getNumPatches()
            nxMax = max(nxMax, nx)
            nyMax = max(nyMax, ny)
        with registry.transaction():
            registry.addDimensionEntry(
                "skymap",
                {"skymap": name,
                 "hash": self.getSha1(),
                 "tract_max": len(self),
                 "patch_nx_max": nxMax,
                 "patch_ny_max": nyMax}
            )
            for tractInfo in self:
                region = tractInfo.getOuterSkyPolygon()
                centroid = SpherePoint(region.getCentroid())
                entry = (
                    {"skymap": name, "tract": tractInfo.getId(),
                     "region": region,
                     "ra": centroid.getRa().asDegrees(),
                     "dec": centroid.getDec().asDegrees()}
                )
                registry.addDimensionEntry("tract", entry)
                patchDataIdList = []
                for patchInfo in tractInfo:
                    cellX, cellY = patchInfo.getIndex()
                    patchDataIdList.append(
                        {"skymap": name, "tract": tractInfo.getId(),
                         "patch": tractInfo.getSequentialPatchIndex(patchInfo),
                         "cell_x": cellX, "cell_y": cellY,
                         "region": patchInfo.getOuterSkyPolygon(tractInfo.getWcs())}
                    )
                registry.addDimensionEntryList("patch", patchDataIdList)
