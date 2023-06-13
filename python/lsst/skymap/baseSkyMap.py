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
import numpy as np

import lsst.geom as geom
import lsst.pex.config as pexConfig
from lsst.geom import SpherePoint, Angle, arcseconds, degrees
from . import detail
from .tractBuilder import tractBuilderRegistry


class BaseSkyMapConfig(pexConfig.Config):
    tractBuilder = tractBuilderRegistry.makeField(
        doc="Algorithm for creating patches within the tract.",
        default="legacy"
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

    # Backwards compatibility
    # We can't use the @property decorator because it makes pexConfig sad.
    def getPatchInnerDimensions(self):
        """Get the patch inner dimensions, for backwards compatibility.

        This value is only used with the ``legacy`` tract builder,
        and is ignored otherwise.  In general, the config should be
        accessed directly with config.tractBuilder["legacy"].patchInnerDimensions.

        Returns
        -------
        innerDimensions : `list` [`int`, `int`]
        """
        return self.tractBuilder["legacy"].patchInnerDimensions

    def setPatchInnerDimensions(self, value):
        """Set the patch inner dimensions, for backwards compatibility.

        This value is only used with the ``legacy`` tract builder,
        and is ignored otherwise.  In general, the config should be
        accessed directly with config.tractBuilder["legacy"].patchInnerDimensions.

        Parameters
        ----------
        value : `list` [`int`, `int`]
        """
        self.tractBuilder["legacy"].patchInnerDimensions = value

    patchInnerDimensions = property(getPatchInnerDimensions, setPatchInnerDimensions)

    def getPatchBorder(self):
        """Get the patch border, for backwards compatibility.

        This value is only used with the ``legacy`` tract builder,
        and is ignored otherwise.  In general, the config should be
        accessed directly with config.tractBuilder["legacy"].patchBorder.

        Returns
        -------
        border: `int`
        """
        return self.tractBuilder["legacy"].patchBorder

    def setPatchBorder(self, value):
        """Set the patch border, for backwards compatibility.

        This value is only used with the ``legacy`` tract builder,
        and is ignored otherwise.  In general, the config should be
        accessed directly with config.tractBuilder["legacy"].patchBorder.

        Parameters
        ----------
        border: `int`
        """
        self.tractBuilder["legacy"].patchBorder = value

    patchBorder = property(getPatchBorder, setPatchBorder)


class BaseSkyMap:
    """A collection of overlapping Tracts that map part or all of the sky.

    See TractInfo for more information.

    Parameters
    ----------
    config : `BaseSkyMapConfig` or `None` (optional)
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
        self._tractBuilder = config.tractBuilder.apply()

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

        .. warning::

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

    def findTractIdArray(self, ra, dec, degrees=False):
        """Find array of tract IDs with vectorized operations (where supported).

        If a given sky map does not support vectorized operations, then a loop
        over findTract will be called.

        Parameters
        ----------
        ra : `numpy.ndarray`
            Array of Right Ascension.  Units are radians unless
            degrees=True.
        dec : `numpy.ndarray`
            Array of Declination.  Units are radians unless
            degrees=True.
        degrees : `bool`, optional
            Input ra, dec arrays are degrees if `True`.

        Returns
        -------
        tractId : `numpy.ndarray`
            Array of tract IDs

        Notes
        -----
        - If coord is equidistant between multiple sky tract centers then one
          is arbitrarily chosen.

        .. warning::

           If tracts do not cover the whole sky then the returned tract may not
           include the given ra/dec.
        """
        units = geom.degrees if degrees else geom.radians
        coords = [geom.SpherePoint(r*units, d*units) for r, d in zip(np.atleast_1d(ra),
                                                                     np.atleast_1d(dec))]

        tractId = np.array([self.findTract(coord).getId() for coord in coords])

        return tractId

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
        .. warning::

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
        log : `logging.Logger`
            Log object that information about skymap will be written.
        """
        log.info("sky map has %s tracts" % (len(self),))
        for tractInfo in self:
            wcs = tractInfo.getWcs()
            posBox = geom.Box2D(tractInfo.getBBox())
            pixelPosList = (
                posBox.getMin(),
                geom.Point2D(posBox.getMaxX(), posBox.getMinY()),
                posBox.getMax(),
                geom.Point2D(posBox.getMinX(), posBox.getMaxY()),
            )
            skyPosList = [wcs.pixelToSky(pos).getPosition(geom.degrees) for pos in pixelPosList]
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
            configPacked = self._tractBuilder.getPackedConfig(self.config)
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

    SKYMAP_RUN_COLLECTION_NAME = "skymaps"

    SKYMAP_DATASET_TYPE_NAME = "skyMap"

    def register(self, name, butler):
        """Add skymap, tract, and patch Dimension entries to the given Gen3
        Butler.

        Parameters
        ----------
        name : `str`
            The name of the skymap.
        butler : `lsst.daf.butler.Butler`
            The butler to add to.

        Raises
        ------
        lsst.daf.butler.registry.ConflictingDefinitionError
            Raised if a different skymap exists with the same name.

        Notes
        -----
        Registering the same skymap multiple times (with the exact same
        definition) is safe, but inefficient; most of the work of computing
        the rows to be inserted must be done first in order to check for
        consistency between the new skymap and any existing one.

        Re-registering a skymap with different tract and/or patch definitions
        but the same summary information may not be detected as a conflict but
        will never result in updating the skymap; there is intentionally no
        way to modify a registered skymap (aside from manual administrative
        operations on the database), as it is hard to guarantee that this can
        be done without affecting reproducibility.
        """
        numPatches = [tractInfo.getNumPatches() for tractInfo in self]
        nxMax = max(nn[0] for nn in numPatches)
        nyMax = max(nn[1] for nn in numPatches)

        skyMapRecord = {
            "skymap": name,
            "hash": self.getSha1(),
            "tract_max": len(self),
            "patch_nx_max": nxMax,
            "patch_ny_max": nyMax,
        }
        butler.registry.registerRun(self.SKYMAP_RUN_COLLECTION_NAME)
        # Kind of crazy that we've got three different capitalizations of
        # "skymap" here, but that's what the various conventions (or at least
        # precedents) dictate.
        from lsst.daf.butler import DatasetType
        from lsst.daf.butler.registry import ConflictingDefinitionError
        datasetType = DatasetType(
            name=self.SKYMAP_DATASET_TYPE_NAME,
            dimensions=["skymap"],
            storageClass="SkyMap",
            universe=butler.dimensions
        )
        butler.registry.registerDatasetType(datasetType)
        with butler.transaction():
            try:
                inserted = butler.registry.syncDimensionData("skymap", skyMapRecord)
            except ConflictingDefinitionError as err:
                raise ConflictingDefinitionError(
                    f"SkyMap with hash {self.getSha1().hex()} is already registered with a different name."
                ) from err
            if inserted:
                for tractInfo in self:
                    tractId = tractInfo.getId()
                    tractRegion = tractInfo.getOuterSkyPolygon()
                    centroid = SpherePoint(tractRegion.getCentroid())
                    tractWcs = tractInfo.getWcs()
                    tractRecord = dict(
                        skymap=name,
                        tract=tractId,
                        region=tractRegion,
                        ra=centroid.getRa().asDegrees(),
                        dec=centroid.getDec().asDegrees(),
                    )
                    butler.registry.insertDimensionData("tract", tractRecord)

                    patchRecords = []
                    for patchInfo in tractInfo:
                        xx, yy = patchInfo.getIndex()
                        patchRecords.append(
                            dict(
                                skymap=name,
                                tract=tractId,
                                patch=tractInfo.getSequentialPatchIndex(patchInfo),
                                cell_x=xx,
                                cell_y=yy,
                                region=patchInfo.getOuterSkyPolygon(tractWcs),
                            )
                        )
                    butler.registry.insertDimensionData("patch", *patchRecords)

                butler.put(self, datasetType, {"skymap": name}, run=self.SKYMAP_RUN_COLLECTION_NAME)
