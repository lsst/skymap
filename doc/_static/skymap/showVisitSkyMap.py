#!/usr/bin/env python
#
# This file is part of skymap.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import logging
import os

import matplotlib
import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord
from lsst.afw.cameraGeom import FIELD_ANGLE, FOCAL_PLANE, DetectorType
from lsst.daf.butler import Butler
from lsst.geom import Angle, Box2D, Point2D, SpherePoint, degrees
from lsst.sphgeom import ConvexPolygon
from matplotlib.legend import Legend

logger = logging.getLogger("lsst.skymap.bin.showVisitSkyMap")


def bboxToRaDec(bbox, wcs):
    """Get the corners of a BBox and convert them to lists of RA and Dec."""
    sphPoints = wcs.pixelToSky(Box2D(bbox).getCorners())
    ra = [float(sph.getRa().asDegrees()) for sph in sphPoints]
    dec = [float(sph.getDec().asDegrees()) for sph in sphPoints]
    return ra, dec


def getValueAtPercentile(values, percentile=0.5):
    """Return a value a fraction of the way between the min and max values in a
    list.

    Parameters
    ----------
    values : `list` [`float`]
        The list of values under consideration.
    percentile : `float`, optional
        The percentile (expressed as a number between 0.0 and 1.0) at which
        to determine the value in `values`.

    Returns
    -------
    result : `float`
        The value at the given percentile of ``values``.
    """
    m = min(values)
    interval = max(values) - m
    return m + percentile * interval


def get_cmap(n, name="gist_rainbow"):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color.

    Uses ``gist_rainbow`` by default: vivid, non-cyclic, so first and last
    colors are visually distinct.
    """
    return matplotlib.colormaps[name].resampled(n)


def queryImageDatasets(butler, whereStr, imageDatasetType=None):
    """Query image datasets with support for current and legacy names."""
    datasetTypes = (
        [imageDatasetType]
        if imageDatasetType is not None
        else [
            "visit_image",
            "preliminary_visit_image",
            "calexp",
        ]
    )
    for datasetType in datasetTypes:
        logger.info("Querying image dataset type: %s", datasetType)
        dataRefs = list(butler.registry.queryDatasets(datasetType, where=whereStr).expanded())
        if len(dataRefs) > 0:
            return datasetType, dataRefs
    logger.warning("No visit-level data found for any of the tested dataset types; unable to plot visits.")
    return None, []


def getVisitSummaryForVisit(butler, visit, visitSummaryDatasetType=None):
    """Fetch visit summary for a visit, supporting legacy and newer names."""
    datasetTypes = (
        [visitSummaryDatasetType]
        if visitSummaryDatasetType is not None
        else [
            "visit_summary",
            "preliminary_visit_summary",
            "visitSummary",
        ]
    )
    for datasetType in datasetTypes:
        try:
            return butler.get(datasetType, visit=visit), datasetType
        except LookupError:
            pass
    raise LookupError(f"Visit summary for visit {visit!r} not found in any of {datasetTypes!r}.")


def main(
    repo,
    collections,
    skymapName=None,
    tracts=None,
    patches=None,
    visits=None,
    physicalFilters=None,
    bands=None,
    ccds=None,
    ccdKey="detector",
    showPatch=False,
    showPatchSelectedTractsOnly=False,
    saveFile=None,
    showCcds=False,
    showCcdsAll=False,
    plotFailsOnly=False,
    visitVetoFile=None,
    minOverlapFraction=None,
    trimToTracts=False,
    trimToOverlappingTracts=False,
    doUnscaledLimitRatio=False,
    forceScaledLimitRatio=False,
    imageDatasetType=None,
    visitSummaryDatasetType=None,
    dpi=150,
    maxVisitForLegend=20,
):
    if minOverlapFraction is not None and tracts is None:
        raise RuntimeError("Must specify --tracts if --minOverlapFraction is set")
    if dpi <= 0:
        raise RuntimeError("--dpi must be > 0")
    if maxVisitForLegend < 0:
        raise RuntimeError("--maxVisitForLegend must be >= 0")
    logger.info("Instantiating butler for repo '%s' with collections = %s", repo, collections)
    butler = Butler.from_config(repo, collections=collections)
    cameraDataset = butler.find_dataset("camera")
    if cameraDataset is None:
        raise RuntimeError("Could not find required dataset type: camera")
    instrument = str(cameraDataset.dataId["instrument"])
    detectorSkipList = []
    # Make a guess at the skymapName if not provided
    if skymapName is None:
        if instrument == "HSC":
            skymapName = "hsc_rings_v1"
            detectorSkipList = [9]  # detector 9 has long been dead for HSC
        elif instrument == "LSSTCam-imSim":
            skymapName = "DC2"
        elif instrument == "LSSTComCamSim":
            skymapName = "ops_rehersal_prep_2k_v1"
        elif instrument == "LATISS":
            skymapName = "latiss_v1"
        elif instrument == "DECam":
            skymapName = "decam_rings_v1"
        elif instrument == "LSSTComCam":
            skymapName = "lsst_cells_v1"
        elif instrument == "LSSTCam":
            skymapName = "lsst_cells_v2"
        else:
            raise RuntimeError(
                f"Unknown skymapName for instrument: {instrument}. Must specify --skymapName on command line."
            )

    logger.info("Using instrument = '%s' and skymapName = '%s'", instrument, skymapName)
    camera = butler.get("camera", instrument=instrument)
    skymap = butler.get("skyMap", instrument=instrument, skymap=skymapName)

    whereStr = ""
    if tracts is not None:
        tractStr = makeWhereInStr("tract", tracts, int)
        whereStr += tractStr

    if patches is not None:
        patchStr = makeWhereInStr("patch", patches, int)
        whereStr += " AND " + patchStr

    if visits is not None:
        visitStr = makeWhereInStr("exposure", visits, int)
        if len(whereStr) < 1:
            whereStr += visitStr
        else:
            whereStr += " AND " + visitStr

    if physicalFilters is not None:
        physicalFilterStr = makeWhereInStr("physical_filter", physicalFilters, str)
        whereStr += " AND " + physicalFilterStr if len(whereStr) else " " + physicalFilterStr

    if bands is not None:
        bandStr = makeWhereInStr("band", bands, str)
        whereStr += " AND " + bandStr if len(whereStr) else " " + bandStr

    if len(whereStr) > 1:
        whereStr = f"instrument='{instrument}' AND skymap='{skymapName}' AND {whereStr}"
    else:
        whereStr = f"instrument='{instrument}' AND skymap='{skymapName}'"
    logger.info("Querying the butler with the following dataId where clause: %s", whereStr)

    imageDatasetTypeUsed, imageDataRefs = queryImageDatasets(
        butler, whereStr, imageDatasetType=imageDatasetType
    )
    logger.info("Using image dataset type: %s", imageDatasetTypeUsed)

    failedDataIds = set()
    if plotFailsOnly:
        processedDataIds = set((ref.dataId["visit"], ref.dataId["detector"]) for ref in imageDataRefs)
        selectedTracts = set(tracts) if tracts is not None else None
        logger.info("Querying raw datasets to find failed detectors...")
        rawDataRefs = list(butler.registry.queryDatasets("raw", where=whereStr).expanded())
        logger.info("Found %d raw datasets", len(rawDataRefs))
        for ref in rawDataRefs:
            exposure = ref.dataId["exposure"]
            detector = ref.dataId["detector"]
            if (exposure, detector) not in processedDataIds:
                if selectedTracts is not None:
                    raCorners, decCorners = getDetRaDecCorners(
                        ccdKey,
                        detector,
                        exposure,
                        visitSummary=None,
                        butler=butler,
                        imageDatasetType="raw",
                        doLogWarn=False,
                    )
                    if raCorners is None or decCorners is None:
                        continue
                    finiteCornerPairs = [
                        (ra, dec)
                        for ra, dec in zip(raCorners, decCorners)
                        if np.isfinite(ra) and np.isfinite(dec)
                    ]
                    if len(finiteCornerPairs) == 0:
                        continue
                    rawRas = [ra for ra, _ in finiteCornerPairs]
                    rawDecs = [dec for _, dec in finiteCornerPairs]
                    rawTracts = set(skymap.findTractIdArray(rawRas, rawDecs, degrees=True))
                    if selectedTracts.isdisjoint(rawTracts):
                        continue
                failedDataIds.add((exposure, detector))
        logger.info(
            "Found %d failed (raw but unprocessed) visit-detector pairs",
            len(failedDataIds),
        )

    visitSummaryDatasetTypeUsed = visitSummaryDatasetType

    visitVetoList = []
    if visitVetoFile is not None:
        with open(visitVetoFile) as f:
            content = f.readlines()
            visitVetoList = [int(visit.strip()) for visit in content]

    if visits is None:
        visits = []
        if plotFailsOnly:
            for failVisit, _ in failedDataIds:
                if failVisit not in visits and failVisit not in visitVetoList:
                    visits.append(failVisit)
        else:
            for dataRef in imageDataRefs:
                visit = dataRef.dataId.visit.id
                if visit not in visits and visit not in visitVetoList:
                    visits.append(visit)
        visits.sort()
        logger.info("List of visits (N=%d) satisfying where and veto clauses: %s", len(visits), visits)
    else:
        if len(visitVetoList) > 1:
            visitListTemp = visits.copy()
            for visit in visitListTemp:
                if visit in visitVetoList:
                    visits.remove(visit)
            logger.info("List of visits (N=%d) excluding veto list: %s}", len(visits), visits)
        logger.info("List of visits (N=%d): %s", len(visits), visits)

    if len(visits) > 0:
        try:
            _, visitSummaryDatasetTypeUsed = getVisitSummaryForVisit(
                butler, visits[0], visitSummaryDatasetType=visitSummaryDatasetType
            )
            logger.info("Using visit summary dataset type: %s", visitSummaryDatasetTypeUsed)
        except LookupError:
            if visitSummaryDatasetType is None:
                logger.info(
                    "No visit summary dataset type found for auto-detection; "
                    "will fall back to detector-level WCS/bbox lookups."
                )
            else:
                logger.info(
                    "Configured visit summary dataset type '%s' was not found for sampled visit; "
                    "will fall back to detector-level WCS/bbox lookups.",
                    visitSummaryDatasetType,
                )

    ccdIdList = []
    for ccd in camera:
        ccdId = ccd.getId()
        if (
            (ccds is None or ccdId in ccds)
            and ccd.getType() == DetectorType.SCIENCE
            and ccdId not in detectorSkipList
        ):
            ccdIdList.append(ccdId)
    ccdIdList.sort()
    nDetTot = len(ccdIdList)
    missingVisitSummaryRows = {}
    nonFiniteDetectorCorners = {}

    visitIncludeList = []
    # Determine the fraction of detectors that overlap any tract under
    # consideration. If this fraction does not exceed minOverlapFraction,
    # skip this visit.
    if minOverlapFraction is not None:
        for i_v, visit in enumerate(visits):
            ccdOverlapList = []
            try:
                visitSummary, _ = getVisitSummaryForVisit(
                    butler, visit, visitSummaryDatasetType=visitSummaryDatasetTypeUsed
                )
            except LookupError as e:
                logger.warning("%s  Will try to get wcs from %s.", e, imageDatasetTypeUsed)
                visitSummary = None
            if tracts is not None:
                for tract in tracts:
                    tractInfo = skymap[tract]
                    sphCorners = tractInfo.wcs.pixelToSky(Box2D(tractInfo.bbox).getCorners())
                    tractConvexHull = ConvexPolygon.convexHull([coord.getVector() for coord in sphCorners])
                    for ccdId in ccdIdList:
                        if ccdId not in ccdOverlapList:
                            raCorners, decCorners = getDetRaDecCorners(
                                ccdKey,
                                ccdId,
                                visit,
                                visitSummary=visitSummary,
                                butler=butler,
                                imageDatasetType=imageDatasetTypeUsed,
                                doLogWarn=False,
                            )
                            if raCorners is not None and decCorners is not None:
                                finiteCornerPairs = [
                                    (ra, dec)
                                    for ra, dec in zip(raCorners, decCorners)
                                    if np.isfinite(ra) and np.isfinite(dec)
                                ]
                                if len(finiteCornerPairs) < 3:
                                    logger.debug(
                                        "visit %d det %d (tract-overlap path): only %d finite "
                                        "corner(s); raw ra=%s dec=%s",
                                        visit,
                                        ccdId,
                                        len(finiteCornerPairs),
                                        raCorners,
                                        decCorners,
                                    )
                                    continue
                                detSphCorners = []
                                for ra, dec in finiteCornerPairs:
                                    pt = SpherePoint(Angle(ra, degrees), Angle(dec, degrees))
                                    detSphCorners.append(pt)
                                try:
                                    detConvexHull = ConvexPolygon.convexHull(
                                        [coord.getVector() for coord in detSphCorners]
                                    )
                                except ValueError as e:
                                    logger.debug(
                                        "visit %d det %d (tract-overlap path): hull ValueError (%s); "
                                        "corners ra=%s dec=%s",
                                        visit,
                                        ccdId,
                                        e,
                                        raCorners,
                                        decCorners,
                                    )
                                    continue
                                if tractConvexHull.contains(detConvexHull):
                                    ccdOverlapList.append(ccdId)

                    if len(ccdOverlapList) / nDetTot >= minOverlapFraction:
                        break
                if len(ccdOverlapList) / nDetTot < minOverlapFraction:
                    logger.info(
                        "Fraction of detectors overlapping any tract for visit %d (%.2f) < "
                        "minimum required (%.2f).  Skipping visit...",
                        visit,
                        len(ccdOverlapList) / nDetTot,
                        minOverlapFraction,
                    )
                else:
                    if visit not in visitIncludeList:
                        visitIncludeList.append(visit)
    else:
        visitIncludeList = visits

    # Draw the CCDs.
    ras, decs = [], []
    ccdBBoxesPlotted = []  # Bounding boxes for CCDs that were actually drawn
    ccdLabelsToPlot = []  # Store CCD label info to plot after limits finalized
    cmap = get_cmap(len(visitIncludeList))
    alphaEdge = 0.7
    finalVisitList = []
    finalVisitColorIndices = []
    includedBands = []
    includedPhysicalFilters = []
    for i_v, visit in enumerate(visitIncludeList):
        print("Working on visit %d [%d of %d]" % (visit, i_v + 1, len(visitIncludeList)), end="\r")
        inLegend = False
        color = cmap(i_v)
        fillKwargs = {"fill": False, "alpha": alphaEdge, "facecolor": None, "edgecolor": color, "lw": 0.6}
        try:
            visitSummary, _ = getVisitSummaryForVisit(
                butler, visit, visitSummaryDatasetType=visitSummaryDatasetTypeUsed
            )
        except Exception as e:
            logger.warning("%s  Will try to get wcs from %s.", e, imageDatasetTypeUsed)
            visitSummary = None

        band, physicalFilter = getBand(visitSummary=visitSummary, butler=butler, visit=visit)
        if band not in includedBands:
            includedBands.append(band)
        if physicalFilter not in includedPhysicalFilters:
            includedPhysicalFilters.append(physicalFilter)

        for ccdId in ccdIdList:
            if plotFailsOnly and (visit, ccdId) not in failedDataIds:
                continue
            raCorners, decCorners = getDetRaDecCorners(
                ccdKey,
                ccdId,
                visit,
                visitSummary=None if plotFailsOnly else visitSummary,
                butler=butler,
                imageDatasetType="raw" if plotFailsOnly else imageDatasetTypeUsed,
                missingVisitSummaryRows=missingVisitSummaryRows,
            )
            if raCorners is not None and decCorners is not None:
                cornerPairs = list(zip(raCorners, decCorners))
                finiteCornerPairs = [
                    (ra, dec) for ra, dec in cornerPairs if np.isfinite(ra) and np.isfinite(dec)
                ]
                if len(finiteCornerPairs) < len(cornerPairs):
                    nonFiniteDetectorCorners.setdefault(visit, set()).add(ccdId)
                    # Skip plotting malformed polygons for this detector.
                    continue

                ras += raCorners
                decs += decCorners
                if not inLegend and len(visitIncludeList) <= maxVisitForLegend:
                    plt.fill(raCorners, decCorners, label=str(visit), **fillKwargs)
                    inLegend = True
                else:
                    plt.fill(raCorners, decCorners, **fillKwargs)
                plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge / 4, color=color, edgecolor=color)
                if visit not in finalVisitList:
                    finalVisitList.append(visit)
                    finalVisitColorIndices.append(i_v)
                # Collect bboxes and CCD label info for later
                if showCcds or showCcdsAll:
                    ccdCenterRa = getValueAtPercentile(raCorners)
                    ccdCenterDec = getValueAtPercentile(decCorners)
                    overlapFrac = 0.2
                    deltaRa = max(raCorners) - min(raCorners)
                    deltaDec = max(decCorners) - min(decCorners)
                    minPoint = Point2D(
                        min(raCorners) + overlapFrac * deltaRa, min(decCorners) + overlapFrac * deltaDec
                    )
                    maxPoint = Point2D(
                        max(raCorners) - overlapFrac * deltaRa, max(decCorners) - overlapFrac * deltaDec
                    )
                    bboxDouble = Box2D(minPoint, maxPoint)
                    ccdLabelsToPlot.append(
                        {
                            "ra": ccdCenterRa,
                            "dec": ccdCenterDec,
                            "ccdId": ccdId,
                            "bbox": bboxDouble,
                        }
                    )

        if visit in missingVisitSummaryRows:
            missingDetectors = sorted(missingVisitSummaryRows[visit])
            logger.warning(
                "visit summary table for visit %d missing detectors: %s",
                visit,
                missingDetectors,
            )

        if visit in nonFiniteDetectorCorners:
            badDetectors = sorted(nonFiniteDetectorCorners[visit])
            logger.warning(
                "Non-finite detector corners for visit %d (N=%d): %s",
                visit,
                len(badDetectors),
                badDetectors,
            )

    logger.info(
        "Final list of visits (N=%d) satisfying where and minOverlapFraction clauses: %s",
        len(finalVisitList),
        finalVisitList,
    )

    raToDecLimitRatio = None
    midRa = None
    midDec = None
    midSkyCoord = None
    if len(ras) > 0:
        finiteCoordPairs = [(ra, dec) for ra, dec in zip(ras, decs) if np.isfinite(ra) and np.isfinite(dec)]
        droppedCoordCount = len(ras) - len(finiteCoordPairs)
        if droppedCoordCount > 0:
            logger.warning(
                "Dropping %d non-finite detector corner coordinates before tract lookup",
                droppedCoordCount,
            )
        if len(finiteCoordPairs) == 0:
            if tracts is not None:
                logger.info(
                    "No finite detector corners found, but --tracts list was provided, so will go ahead and "
                    "plot the empty tracts."
                )
                tractList = tracts
                trimToTracts = True
            else:
                raise RuntimeError("No finite detector corners available for tract lookup")
        else:
            ras, decs = zip(*finiteCoordPairs)
            ras = list(ras)
            decs = list(decs)
            tractList = list(set(skymap.findTractIdArray(ras, decs, degrees=True)))
            minVisitRa, maxVisitRa = min(ras), max(ras)
            minVisitDec, maxVisitDec = min(decs), max(decs)
            raVisitDiff = maxVisitRa - minVisitRa
            decVisitDiff = maxVisitDec - minVisitDec
            midVisitRa = minVisitRa + 0.5 * raVisitDiff
            midVisitDec = minVisitDec + 0.5 * decVisitDiff
            midRa = np.atleast_1d((midVisitRa * units.deg).to(units.radian).value).astype(np.float64)
            midDec = np.atleast_1d((midVisitDec * units.deg).to(units.radian).value).astype(np.float64)
            midSkyCoord = SkyCoord(midVisitRa * units.deg, midVisitDec * units.deg)
    else:
        if tracts is not None:
            logger.info(
                "No detectors were found, but --tracts list was provided, so will go ahead and "
                "plot the empty tracts."
            )
            tractList = tracts
            trimToTracts = True
        else:
            raise RuntimeError(
                "No data to plot (if you want to plot empty tracts, include them as "
                "a blank-space separated list to the --tracts option)."
            )
    tractList, invalidTracts = sanitizeTractList(skymap, tractList)
    if len(invalidTracts) > 0:
        logger.warning("Ignoring invalid tract ids: %s", invalidTracts)
    if len(tractList) == 0:
        raise RuntimeError("No valid tract ids found for plotting")
    logger.info("List of tracts overlapping data:  %s", tractList)

    # Determine which tracts to use for plot axis limits.
    # trimToOverlappingTracts always uses all tracts overlapping the visits;
    # trimToTracts uses user tracts when available, else falls back
    # to overlapping tracts; otherwise plot limits come from visit footprints.
    if trimToOverlappingTracts:
        tractLimitsForPlotting = tractList
    elif trimToTracts:
        tractLimitsForPlotting = tracts if tracts is not None else tractList
    else:
        tractLimitsForPlotting = None  # use visit footprints for plot limits

    # Compute the tract limits dict only if needed for plot limits.
    if tractLimitsForPlotting is not None:
        tractLimitsDict = getTractLimitsDict(skymap, tractLimitsForPlotting)
    else:
        tractLimitsDict = None

    if forceScaledLimitRatio:
        doUnscaledLimitRatio = False
    else:
        # Roughly compute radius in degrees of a single detector.  If RA/Dec
        # coverage is more than 30 times the detector radius, and the RA/Dec
        # limit ratio is greater than raDecScaleThresh, don't try to scale to
        # detector coords.
        radiusMm = camera.computeMaxFocalPlaneRadius()
        fpRadiusPt = Point2D(radiusMm, radiusMm)
        focalPlaneToFieldAngle = camera.getTransformMap().getTransform(FOCAL_PLANE, FIELD_ANGLE)
        fpRadiusDeg = np.rad2deg(focalPlaneToFieldAngle.applyForward(fpRadiusPt))[0]
        detectorRadiusDeg = fpRadiusDeg / np.sqrt(len(camera))

        if tractLimitsDict is not None:
            xLimMin, xLimMax, yLimMin, yLimMax = getMinMaxLimits(tractLimitsDict)
            xDelta0 = xLimMax - xLimMin
            yDelta0 = yLimMax - yLimMin
        else:
            xDelta0 = raVisitDiff
            yDelta0 = decVisitDiff
            yLimMin = minVisitDec
            yLimMax = maxVisitDec
        raDecScaleThresh = 1.5  # This is a best guess with current testing.
        if (
            (xDelta0 / yDelta0 > raDecScaleThresh or yDelta0 / xDelta0 > raDecScaleThresh)
            and max(xDelta0, yDelta0) > 70 * detectorRadiusDeg
            and yLimMin < 75.0
            and yLimMax > -75.0
        ):
            logger.info(
                "Sky coverage is large (and not too close to a pole), so not scaling to detector coords."
            )
            doUnscaledLimitRatio = True

    if not doUnscaledLimitRatio and midSkyCoord is not None:
        # Find a detector that contains the mid point in RA/Dec (or the closest
        # one) to set the plot aspect ratio.
        minDistToMidCoord = 1e12
        minSepVisit = None
        minSepCcdId = None
        for i_v, visit in enumerate(visits):
            try:
                visitSummary, _ = getVisitSummaryForVisit(
                    butler, visit, visitSummaryDatasetType=visitSummaryDatasetTypeUsed
                )
            except Exception as e:
                logger.warning("%s  Will try to get wcs from %s.", e, imageDatasetTypeUsed)
                visitSummary = None
            for ccdId in ccdIdList:
                raCorners, decCorners = getDetRaDecCorners(
                    ccdKey,
                    ccdId,
                    visit,
                    visitSummary=visitSummary,
                    butler=butler,
                    imageDatasetType=imageDatasetTypeUsed,
                    doLogWarn=False,
                )
                if raCorners is not None and decCorners is not None:
                    finiteCornerPairs = [
                        (ra, dec)
                        for ra, dec in zip(raCorners, decCorners)
                        if np.isfinite(ra) and np.isfinite(dec)
                    ]
                    if len(finiteCornerPairs) < 3:
                        logger.debug(
                            "visit %d det %d (midpoint path): only %d finite corner(s); raw ra=%s dec=%s",
                            visit,
                            ccdId,
                            len(finiteCornerPairs),
                            raCorners,
                            decCorners,
                        )
                        continue
                    detSphCorners = []
                    for ra, dec in finiteCornerPairs:
                        pt = SpherePoint(Angle(ra, degrees), Angle(dec, degrees))
                        detSphCorners.append(pt)
                        ptSkyCoord = SkyCoord(ra * units.deg, dec * units.deg)
                        separation = (midSkyCoord.separation(ptSkyCoord)).degree
                        if separation < minDistToMidCoord:
                            minSepVisit = visit
                            minSepCcdId = ccdId
                            minDistToMidCoord = separation
                    try:
                        detConvexHull = ConvexPolygon([coord.getVector() for coord in detSphCorners])
                    except ValueError as e:
                        logger.debug(
                            "visit %d det %d (midpoint path): hull ValueError (%s); corners ra=%s dec=%s",
                            visit,
                            ccdId,
                            e,
                            raCorners,
                            decCorners,
                        )
                        continue
                    if detConvexHull.contains(midRa, midDec) and raToDecLimitRatio is None:
                        logger.info(
                            "visit/det overlapping plot coord mid point in RA/Dec: %d %d", visit, ccdId
                        )
                        raToDecLimitRatio = (max(raCorners) - min(raCorners)) / (
                            max(decCorners) - min(decCorners)
                        )
                        det = camera[ccdId]
                        width = det.getBBox().getWidth()
                        height = det.getBBox().getHeight()
                        if raToDecLimitRatio > 1.0:
                            raToDecLimitRatio /= max(height / width, width / height)
                        else:
                            if raToDecLimitRatio < 1.0:
                                raToDecLimitRatio *= max(height / width, width / height)
                        break
            if raToDecLimitRatio is not None:
                break

        if raToDecLimitRatio is None and minSepVisit is not None:
            canComputeNearestDetRatio = True
            try:
                visitSummary, _ = getVisitSummaryForVisit(
                    butler, minSepVisit, visitSummaryDatasetType=visitSummaryDatasetTypeUsed
                )
            except Exception as e:
                logger.warning("%s  Will try to get wcs from %s.", e, imageDatasetTypeUsed)
                visitSummary = None
            raCorners, decCorners = getDetRaDecCorners(
                ccdKey,
                minSepCcdId,
                minSepVisit,
                visitSummary=visitSummary,
                butler=butler,
                imageDatasetType=imageDatasetTypeUsed,
                doLogWarn=False,
            )
            if raCorners is None or decCorners is None:
                logger.warning(
                    "Could not determine detector corners for visit/det nearest plot midpoint: %d %d",
                    minSepVisit,
                    minSepCcdId,
                )
                canComputeNearestDetRatio = False
            else:
                finiteCornerPairs = [
                    (ra, dec)
                    for ra, dec in zip(raCorners, decCorners)
                    if np.isfinite(ra) and np.isfinite(dec)
                ]
                if len(finiteCornerPairs) < 3:
                    logger.warning(
                        "Insufficient finite detector corners for visit/det nearest plot midpoint: %d %d",
                        minSepVisit,
                        minSepCcdId,
                    )
                    canComputeNearestDetRatio = False
                else:
                    raCorners = [ra for ra, _ in finiteCornerPairs]
                    decCorners = [dec for _, dec in finiteCornerPairs]
            if canComputeNearestDetRatio:
                logger.info(
                    "visit/det closest to plot coord mid point in RA/Dec (none actually overlap it): %d %d",
                    minSepVisit,
                    minSepCcdId,
                )
                raToDecLimitRatio = (max(raCorners) - min(raCorners)) / (max(decCorners) - min(decCorners))
                det = camera[minSepCcdId]
                width = det.getBBox().getWidth()
                height = det.getBBox().getHeight()
                if raToDecLimitRatio > 1.0:
                    raToDecLimitRatio /= max(height / width, width / height)
                else:
                    if raToDecLimitRatio < 1.0:
                        raToDecLimitRatio *= max(height / width, width / height)

    elif not doUnscaledLimitRatio:
        logger.info(
            "Skipping midpoint-based detector aspect-ratio scaling: "
            "no detector-derived midpoint is available."
        )

    if plotFailsOnly and not doUnscaledLimitRatio and raToDecLimitRatio is None:
        # In fail-only mode we can lack midpoint detector context; default to
        # 1:1 sky-coordinate limits for predictable geometry.
        raToDecLimitRatio = 1.0

    if tractLimitsDict is not None:
        xlim, ylim = derivePlotLimits(tractLimitsDict, raToDecLimitRatio=raToDecLimitRatio, buffFrac=0.04)
    else:
        visitLimitsDict = {"allVisits": {"ras": [minVisitRa, maxVisitRa], "decs": [minVisitDec, maxVisitDec]}}
        xlim, ylim = derivePlotLimits(visitLimitsDict, raToDecLimitRatio=raToDecLimitRatio, buffFrac=0.04)

    if doUnscaledLimitRatio:
        boxAspectRatio = abs((ylim[1] - ylim[0]) / (xlim[1] - xlim[0]))
    else:
        boxAspectRatio = 1.0

    # Plot deferred CCD labels after plot limits are finalized.
    ccdLabelEdgeBuffer = 0.03  # fraction of plot width/height to use as inset
    ccdBufRa = ccdLabelEdgeBuffer * abs(xlim[1] - xlim[0])
    ccdBufDec = ccdLabelEdgeBuffer * abs(ylim[1] - ylim[0])
    for ccdLabel in ccdLabelsToPlot:
        # Keep center inside buffered bounds (RA xlim may be inverted).
        if (
            min(xlim) + ccdBufRa < ccdLabel["ra"] < max(xlim) - ccdBufRa
            and min(ylim) + ccdBufDec < ccdLabel["dec"] < max(ylim) - ccdBufDec
        ):
            if showCcdsAll:
                plt.text(
                    ccdLabel["ra"],
                    ccdLabel["dec"],
                    str(ccdLabel["ccdId"]),
                    fontsize=6,
                    ha="center",
                    va="center",
                    color="darkblue",
                )
            else:
                overlaps = [ccdLabel["bbox"].overlaps(otherBbox) for otherBbox in ccdBBoxesPlotted]
                if not any(overlaps):
                    plt.text(
                        ccdLabel["ra"],
                        ccdLabel["dec"],
                        str(ccdLabel["ccdId"]),
                        fontsize=6,
                        ha="center",
                        va="center",
                        color="darkblue",
                    )
                    ccdBBoxesPlotted.append(ccdLabel["bbox"])

    # Draw the skymap.
    alpha0 = 1.0
    tractHandleList = []
    tractStrList = []
    if tracts is not None:
        tractOutlineList = list(set(tracts + tractList))
    else:
        tractOutlineList = tractList
    tractOutlineList.sort()
    selectedPatchTracts = set(tracts) if tracts is not None else None
    if showPatchSelectedTractsOnly and selectedPatchTracts is None:
        logger.info(
            "--showPatchSelectedTractsOnly was set without --tracts; showing patches for all plotted tracts."
        )
    logger.info("List of tract outlines being plotted: %s", tractOutlineList)
    for i_t, tract in enumerate(tractOutlineList):
        alpha = max(0.1, alpha0 - i_t * 1.0 / len(tractOutlineList))
        tractInfo = skymap[tract]
        tCenter = tractInfo.ctr_coord
        tCenterRa = tCenter.getRa().asDegrees()
        tCenterDec = tCenter.getDec().asDegrees()
        fracDeltaX = 0.02 * abs((xlim[1] - xlim[0]))
        fracDeltaY = 0.02 * abs((ylim[1] - ylim[0]))
        if (
            xlim[1] + fracDeltaX < tCenterRa < xlim[0] - fracDeltaX
            and ylim[0] + fracDeltaY < tCenterDec < ylim[1] - fracDeltaY
        ):
            if len(tractOutlineList) > 1 or not showPatch:
                if not showPatch:
                    plt.text(tCenterRa, tCenterDec, tract, fontsize=7, alpha=alpha, ha="center", va="center")
                else:
                    plt.text(
                        tCenterRa,
                        tCenterDec,
                        tract,
                        fontsize=7,
                        alpha=1,
                        color="white",
                        path_effects=[pathEffects.withStroke(linewidth=3, foreground="black")],
                        fontweight=500,
                        ha="center",
                        va="center",
                        zorder=5,
                    )
        ra, dec = bboxToRaDec(tractInfo.bbox, tractInfo.getWcs())
        plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed", alpha=alpha)
        tractArtist = matplotlib.patches.Patch(fill=False, edgecolor="k", linestyle="dashed", alpha=alpha)
        tractHandleList.append(tractArtist)
        tractStrList.append(str(tract))
        if showPatch and (
            not showPatchSelectedTractsOnly or selectedPatchTracts is None or tract in selectedPatchTracts
        ):
            patchColor = "k"
            for patch in tractInfo:
                ra, dec = bboxToRaDec(patch.getInnerBBox(), tractInfo.getWcs())
                plt.fill(ra, dec, fill=False, edgecolor=patchColor, lw=0.5, linestyle="dotted", alpha=alpha)
                if (
                    xlim[1] + fracDeltaX < getValueAtPercentile(ra) < xlim[0] - fracDeltaX
                    and ylim[0] + fracDeltaY < getValueAtPercentile(dec) < ylim[1] - fracDeltaY
                ):
                    plt.text(
                        getValueAtPercentile(ra),
                        getValueAtPercentile(dec),
                        str(patch.sequential_index),
                        fontsize=5,
                        color=patchColor,
                        ha="center",
                        va="center",
                        alpha=alpha,
                    )

    # Add labels and save.
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_box_aspect(boxAspectRatio)
    if abs(xlim[1] > 99.99):
        ax.tick_params("x", labelrotation=45, pad=0, labelsize=8)
    else:
        ax.tick_params("x", labelrotation=0, pad=0, labelsize=8)
    ax.tick_params("y", labelsize=8)
    ax.set_xlabel("RA (deg)", fontsize=9)
    ax.set_ylabel("Dec (deg)", fontsize=9)

    if len(visitIncludeList) > maxVisitForLegend:
        nVisits = len(finalVisitList)
        nz = matplotlib.colors.Normalize(vmin=0, vmax=len(visitIncludeList) - 1)
        cax, _ = matplotlib.colorbar.make_axes(plt.gca(), pad=0.03)
        cax.tick_params(labelsize=7)
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=nz, alpha=alphaEdge)
        nTicks = min(8, nVisits)
        tickPositions = [finalVisitColorIndices[int(round(x))] for x in np.linspace(0, nVisits - 1, nTicks)]
        cb.set_ticks(tickPositions)
        cb.set_ticklabels([str(finalVisitList[int(round(x))]) for x in np.linspace(0, nVisits - 1, nTicks)])
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
        cb.set_label("visit", rotation=-90, labelpad=13, fontsize=9)
        tractLegend = Legend(
            ax,
            tractHandleList,
            tractStrList,
            loc="upper right",
            fancybox=True,
            shadow=True,
            fontsize=5,
            title_fontsize=6,
            title="tracts",
        )
        ax.add_artist(tractLegend)
    else:
        if len(visitIncludeList) > 0:
            handles, labels = ax.get_legend_handles_labels()
            if len(handles) > 0:
                xBboxAnchor = min(1.25, max(1.03, boxAspectRatio * 1.15))
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(xBboxAnchor, 0.5),
                    fancybox=True,
                    shadow=True,
                    fontsize=6,
                    title_fontsize=6,
                    title="visits",
                )
        # Create the second legend and add the artist manually.
        tractLegend = Legend(
            ax,
            tractHandleList,
            tractStrList,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            fancybox=True,
            shadow=True,
            fontsize=6,
            title_fontsize=6,
            title="tracts",
        )
        ax.add_artist(tractLegend)

    titleStr = repo + "\n" + collections[0]
    if len(collections) > 1:
        for collection in collections[1:]:
            titleStr += "\n" + collection
    if plotFailsOnly:
        titleStr += "\nFailed calibration: N={}".format(len(failedDataIds))
    else:
        titleStr += "\nnVisit: {}".format(str(len(finalVisitList)))
    if minOverlapFraction is not None:
        titleStr += " (minOverlapFraction = {:.2f})".format(minOverlapFraction)
    if len(includedBands) > 0:
        titleStr += ", bands: {}".format(str(includedBands).translate({ord(i): None for i in "[]'"}))
    if len(includedPhysicalFilters) > 0:
        if len(includedPhysicalFilters[0]) > 9:
            titleStr += "\n"
        else:
            titleStr += ","
        titleStr += " physical filters: {}".format(
            str(includedPhysicalFilters).translate({ord(i): None for i in "[]'"})
        )
    ax.set_title("{}".format(titleStr), fontsize=8)

    fig = plt.gcf()
    if boxAspectRatio > 1.0:
        minInches = max(4.0, 0.3 * abs(xlim[1] - xlim[0]))
        xInches = minInches
        yInches = min(120.0, boxAspectRatio * minInches)
        fig.set_size_inches(xInches, yInches)
    if boxAspectRatio < 1.0:
        minInches = max(4.0, 0.3 * abs(ylim[1] - ylim[0]))
        xInches = min(120.0, minInches / boxAspectRatio)
        yInches = minInches
        fig.set_size_inches(xInches, yInches)
    if saveFile is not None:
        if plotFailsOnly and "fail" not in saveFile:
            fileRoot, fileExt = os.path.splitext(saveFile)
            saveFile = f"{fileRoot}_failed{fileExt}" if fileExt else f"{saveFile}_failed"
        logger.info("Saving file in: %s", saveFile)
        fig.savefig(saveFile, bbox_inches="tight", dpi=dpi)
    else:
        fig.show()


def makeWhereInStr(parameterName, parameterList, parameterType):
    """Create the string to be used in the where clause for registry lookup."""
    typeStr = "'" if parameterType is str else ""
    whereInStr = parameterName + " IN (" + typeStr + str(parameterList[0]) + typeStr
    if len(parameterList) > 1:
        for param in parameterList[1:]:
            whereInStr += ", " + typeStr + str(param) + typeStr
    whereInStr += ")"

    return whereInStr


def sanitizeTractList(skymap, tractList):
    """Split tract ids into valid and invalid entries for the given skymap."""
    validTracts = []
    invalidTracts = []
    for tract in tractList:
        tractInt = int(tract)
        try:
            _ = skymap[tractInt]
        except IndexError:
            invalidTracts.append(tractInt)
            continue
        if tractInt not in validTracts:
            validTracts.append(tractInt)
    validTracts.sort()
    invalidTracts.sort()
    return validTracts, invalidTracts


def getTractLimitsDict(skymap, tractList):
    """Return a dict containing tract limits needed for outline plotting.

    Parameters
    ----------
    skymap : `lsst.skymap.BaseSkyMap`
        The sky map used for this dataset. Used to obtain tract
        parameters.
    tractList : `list` [`int`]
        The list of tract ids (as integers) for which to determine the
        limits.

    Returns
    -------
    tractLimitsDict : `dict` [`dict`]
        A dictionary keyed on tract id. Each entry includes a `dict`
        including the tract RA corners, Dec corners, and the tract center,
        all in units of degrees. These are used for plotting the tract
        outlines.
    """
    tractLimitsDict = {}
    for tract in tractList:
        tractInfo = skymap[tract]
        tractBbox = tractInfo.outer_sky_polygon.getBoundingBox()
        tractCenter = tractBbox.getCenter()
        tractRa0 = (tractCenter[0] - tractBbox.getWidth() / 2).asDegrees()
        tractRa1 = (tractCenter[0] + tractBbox.getWidth() / 2).asDegrees()
        tractDec0 = (tractCenter[1] - tractBbox.getHeight() / 2).asDegrees()
        tractDec1 = (tractCenter[1] + tractBbox.getHeight() / 2).asDegrees()
        tractLimitsDict[tract] = {
            "ras": [tractRa0, tractRa1, tractRa1, tractRa0, tractRa0],
            "decs": [tractDec0, tractDec0, tractDec1, tractDec1, tractDec0],
            "center": [tractCenter[0].asDegrees(), tractCenter[1].asDegrees()],
        }

    return tractLimitsDict


def getMinMaxLimits(limitsDict):
    """Derive the min and max axis limits of points in limitsDict.

    Parameters
    ----------
    limitsDict : `dict` [`dict`]
        A dictionary keyed on any id. Each entry includes a `dict`
        keyed on "ras" and "decs" including (at least the minimum
        and maximum) RA and Dec values in units of degrees.

    Returns
    -------
    xLimMin, xLimMax, yLimMin, yLimMax : `tuple` [`float`]
        The min and max values for the x and y-axis limits, respectively.
    """
    xLimMin, yLimMin = 1e12, 1e12
    xLimMax, yLimMax = -1e12, -1e12
    for limitId, limits in limitsDict.items():
        xLimMin = min(xLimMin, min(limits["ras"]))
        xLimMax = max(xLimMax, max(limits["ras"]))
        yLimMin = min(yLimMin, min(limits["decs"]))
        yLimMax = max(yLimMax, max(limits["decs"]))

    return xLimMin, xLimMax, yLimMin, yLimMax


def derivePlotLimits(limitsDict, raToDecLimitRatio=1.0, buffFrac=0.0):
    """Derive the axis limits to encompass all points in limitsDict.

    Parameters
    ----------
    limitsDict : `dict` [`dict`]
        A dictionary keyed on any id. Each entry includes a `dict`
        keyed on "ras" and "decs" including (at least the minimum
        and maximum) RA and Dec values in units of degrees.
    raToDecLimitRatio : `float`, optional
        The aspect ratio between RA and Dec to set the plot limits to.  This
        is to namely to set this ratio to that of the focal plane (i.e. such
        that a square detector appears as a square), but any aspect ratio can,
        in principle, be requested.

    Returns
    -------
    xlim, ylim : `tuple` [`float`]
        Two tuples containing the derived min and max values for the x and
        y-axis limits (in degrees), respectively.
    """
    xLimMin, xLimMax, yLimMin, yLimMax = getMinMaxLimits(limitsDict)

    xDelta0 = xLimMax - xLimMin
    yDelta0 = yLimMax - yLimMin
    if raToDecLimitRatio is None:
        padFrac = 0.05
        xlim = xLimMax + padFrac * xDelta0, xLimMin - padFrac * xDelta0
        ylim = yLimMin - padFrac * yDelta0, yLimMax + padFrac * yDelta0
        return xlim, ylim

    if raToDecLimitRatio == 1.0:
        if xDelta0 > yDelta0:
            xLimMin -= buffFrac * yDelta0
            xLimMax += buffFrac * yDelta0
        else:
            yLimMin -= buffFrac * yDelta0
            yLimMax += buffFrac * yDelta0
        xLimMin, xLimMax, yLimMin, yLimMax = setLimitsToEqualRatio(xLimMin, xLimMax, yLimMin, yLimMax)
    else:
        xLimMin -= buffFrac * xDelta0
        xLimMax += buffFrac * xDelta0
        yLimMin -= buffFrac * yDelta0
        yLimMax += buffFrac * yDelta0
        xLimMin, xLimMax, yLimMin, yLimMax = setLimitsToEqualRatio(xLimMin, xLimMax, yLimMin, yLimMax)
        xDelta = xLimMax - xLimMin
        yDelta = yLimMax - yLimMin
        if raToDecLimitRatio > 1.0:
            if yDelta0 > xDelta:
                xMid = xLimMin + 0.5 * (xDelta)
                xLimMin = xMid - 0.5 * yDelta * raToDecLimitRatio
                xLimMax = xMid + 0.5 * yDelta * raToDecLimitRatio
            else:
                yMid = yLimMin + 0.5 * (yDelta)
                yLimMin = yMid - 0.5 * xDelta / raToDecLimitRatio
                yLimMax = yMid + 0.5 * xDelta / raToDecLimitRatio
        else:
            if xDelta0 > yDelta0:
                yMid = yLimMin + 0.5 * (yDelta)
                yLimMin = yMid - 0.5 * xDelta / raToDecLimitRatio
                yLimMax = yMid + 0.5 * xDelta / raToDecLimitRatio
            else:
                xMid = xLimMin + 0.5 * (xDelta)
                xLimMin = xMid - 0.5 * yDelta * raToDecLimitRatio
                xLimMax = xMid + 0.5 * yDelta * raToDecLimitRatio
    xlim = xLimMax, xLimMin
    ylim = yLimMin, yLimMax
    return xlim, ylim


def setLimitsToEqualRatio(xMin, xMax, yMin, yMax):
    """For a given set of x/y min/max, redefine to have equal aspect ratio.

    The limits are extended on both ends such that the central value is
    preserved.

    Parameters
    ----------
    xMin, xMax, yMin, yMax : `float`
        The min/max values of the x/y ranges for which to match in dynamic
        range while preserving the central values.

    Returns
    -------
    xMin, xMax, yMin, yMax : `float`
        The adjusted min/max values of the x/y ranges with equal aspect ratios.
    """
    xDelta = xMax - xMin
    yDelta = yMax - yMin
    deltaDiff = yDelta - xDelta
    if deltaDiff > 0:
        xMin -= 0.5 * deltaDiff
        xMax += 0.5 * deltaDiff
    elif deltaDiff < 0:
        yMin -= 0.5 * abs(deltaDiff)
        yMax += 0.5 * abs(deltaDiff)
    return xMin, xMax, yMin, yMax


def getDetRaDecCorners(
    ccdKey,
    ccdId,
    visit,
    visitSummary=None,
    butler=None,
    imageDatasetType="calexp",
    doLogWarn=True,
    missingVisitSummaryRows=None,
):
    """Compute the RA/Dec corners lists for a given detector in a visit."""
    raCorners, decCorners = None, None
    if visitSummary is not None:
        row = visitSummary.find(ccdId)
        if row is None:
            if doLogWarn:
                if missingVisitSummaryRows is not None:
                    missingVisitSummaryRows.setdefault(visit, set()).add(ccdId)
                else:
                    logger.warning(
                        "No row found for %d in visit summary table for visit %d. Skipping and continuing...",
                        ccdId,
                        visit,
                    )
        else:
            raCorners = list(row["raCorners"])
            decCorners = list(row["decCorners"])
    else:
        if butler is None:
            raise RuntimeError("A butler instance is required when visitSummary is not provided")
        try:
            if imageDatasetType == "raw":
                dataId = {"exposure": visit, ccdKey: ccdId}
                # Raw exposures may not have a bbox component.
                exposure = butler.get(imageDatasetType, dataId)
                wcs = exposure.getWcs()
                bbox = exposure.getBBox()
            else:
                dataId = {"visit": visit, ccdKey: ccdId}
                try:
                    wcs = butler.get(f"{imageDatasetType}.wcs", dataId)
                    bbox = butler.get(f"{imageDatasetType}.bbox", dataId)
                except LookupError:
                    # Some dataset types are only available as full exposures.
                    exposure = butler.get(imageDatasetType, dataId)
                    wcs = exposure.getWcs()
                    bbox = exposure.getBBox()
            if wcs is None or bbox is None:
                if doLogWarn:
                    logger.warning(
                        "WCS or BBox is None for datasetType=%s visit=%s %s=%s. Skipping and continuing...",
                        imageDatasetType,
                        visit,
                        ccdKey,
                        ccdId,
                    )
                return None, None
            raCorners, decCorners = bboxToRaDec(bbox, wcs)
        except LookupError as e:
            if doLogWarn:
                logger.warning("%s Skipping and continuing...", e)
        except Exception as e:
            if doLogWarn:
                logger.warning("%s Skipping and continuing...", e)

    return raCorners, decCorners


def getBand(visitSummary=None, butler=None, visit=None):
    """Determine band and physical filter for given visit.

    Parameters
    ----------
    visitSummary : `lsst.afw.table.ExposureCatalog` or `None`, optional
        The visitSummary table for the visit for which to determine the band.
    butler : `lsst.daf.butler.Butler` or `None`, optional
        The butler from which to look up the Dimension Records. Only needed
        if ``visitSummary`` is `None`.
    visit : `int` or `None, optional
        The visit number for which to determine the band. Only needed
        if ``visitSummary`` is `None`.

    Returns
    -------
    band, physicalFilter : `str`
        The band and physical filter for the given visit.
    """
    if visitSummary is not None:
        band = visitSummary[0]["band"]
        physicalFilter = visitSummary[0]["physical_filter"]
    else:
        if butler is None:
            raise RuntimeError("A butler instance is required when visitSummary is not provided")
        record = list(butler.registry.queryDimensionRecords("band", visit=visit))[0]
        band = record.name
        record = list(butler.registry.queryDimensionRecords("physical_filter", visit=visit))[0]
        physicalFilter = record.name
    return band, physicalFilter


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "repo", type=str, help="URI or path to an existing data repository root or configuration file"
    )
    parser.add_argument(
        "--collections",
        type=str,
        nargs="+",
        help="Blank-space separated list of collection names for butler instantiation",
        metavar=("COLLECTION1", "COLLECTION2"),
        required=True,
    )
    parser.add_argument("--skymapName", default=None, help="Name of the skymap for the collection")
    parser.add_argument(
        "--tracts",
        type=int,
        nargs="+",
        default=None,
        help=("Blank-space separated list of tracts to constrain search for visit overlap"),
        metavar=("TRACT1", "TRACT2"),
    )
    parser.add_argument(
        "--patches",
        type=int,
        nargs="+",
        default=None,
        help=("Blank-space separated list of patches to constrain search for visit overlap"),
        metavar=("PATCH1", "PATCH2"),
    )
    parser.add_argument(
        "--visits",
        type=int,
        nargs="+",
        default=None,
        help="Blank-space separated list of visits to include",
        metavar=("VISIT1", "VISIT2"),
    )
    parser.add_argument(
        "--physicalFilters",
        type=str,
        nargs="+",
        default=None,
        help=("Blank-space separated list of physical filter names to constrain search for visits"),
        metavar=("PHYSICAL_FILTER1", "PHYSICAL_FILTER2"),
    )
    parser.add_argument(
        "--bands",
        type=str,
        nargs="+",
        default=None,
        help=("Blank-space separated list of canonical band names to constrain search for visits"),
        metavar=("BAND1", "BAND2"),
    )
    parser.add_argument(
        "-c",
        "--ccds",
        nargs="+",
        type=int,
        default=None,
        help="Blank-space separated list of CCDs to show",
        metavar=("CCD1", "CCD2"),
    )
    parser.add_argument(
        "-p", "--showPatch", action="store_true", default=False, help="Show the patch boundaries"
    )
    parser.add_argument(
        "--showPatchSelectedTractsOnly",
        action="store_true",
        default=False,
        help="When showing patches, only draw patch boundaries for the tracts passed to --tracts.",
    )
    parser.add_argument(
        "--saveFile", type=str, default="showVisitSkyMap.png", help="Filename to write the plot to"
    )
    parser.add_argument("--ccdKey", default="detector", help="Data ID name of the CCD key")
    parser.add_argument(
        "--showCcds",
        action="store_true",
        default=False,
        help="Show ccd ID numbers on output image, suppressing overlapping labels.",
    )
    parser.add_argument(
        "--showCcdsAll",
        action="store_true",
        default=False,
        help="Show all ccd ID numbers on output image (without suppressing overlaps). "
        "Takes precedence over --showCcds when both are set.",
    )
    parser.add_argument(
        "--plotFailsOnly",
        action="store_true",
        default=False,
        help=(
            "Plot only detector outlines for visits that have raw data but no "
            "processed image (i.e. failed calibration). Outlines are drawn "
            "from the raw WCS. Output filename is suffixed with '_failed' automatically."
        ),
    )
    parser.add_argument(
        "--visitVetoFile",
        type=str,
        default=None,
        help="Full path to single-column file containing a list of visits to veto",
    )
    parser.add_argument(
        "--minOverlapFraction",
        type=float,
        default=None,
        help="Minimum fraction of detectors that overlap any tract for visit to be included",
    )
    parser.add_argument(
        "--trimToTracts",
        action="store_true",
        default=False,
        help=(
            "Set plot limits based on the extent of the tracts specified via --tracts "
            "(or all overlapping tracts if --tracts is not given), rather than the visit footprints."
        ),
    )
    parser.add_argument(
        "--trimToOverlappingTracts",
        action="store_true",
        default=False,
        help=(
            "Set plot limits based on the extent of all tracts overlapping the plotted visits, "
            "regardless of --tracts. Takes precedence over --trimToTracts when both are set."
        ),
    )
    parser.add_argument(
        "--doUnscaledLimitRatio",
        action="store_true",
        default=False,
        help="Let axis limits get set by sky coordinate range without scaling to focal "
        "plane based projection (ignored if --forceScaledLimitRatio is passed).",
    )
    parser.add_argument(
        "--forceScaledLimitRatio",
        action="store_true",
        default=False,
        help="Force the axis limit scaling to focal plane based projection (takes "
        "precedence over --doUnscaledLimitRatio.",
    )
    parser.add_argument(
        "--imageDatasetType",
        type=str,
        default=None,
        help=(
            "Image dataset type used for visit discovery and WCS/bbox fallback; "
            "defaults to commonly used detector storage types."
        ),
    )
    parser.add_argument(
        "--visitSummaryDatasetType",
        type=str,
        default=None,
        help=("Visit summary dataset type to use; defaults to commonly used visit summary types."),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="DPI used when writing output via --saveFile.",
    )
    parser.add_argument(
        "--maxVisitForLegend",
        type=int,
        default=20,
        help="Maximum number of visits to include in the legend before switching to a colorbar.",
    )
    parser.add_argument(
        "--logLevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO).",
    )
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.logLevel), format="%(levelname)s:%(name)s: %(message)s")
    main(
        args.repo,
        args.collections,
        skymapName=args.skymapName,
        tracts=args.tracts,
        patches=args.patches,
        visits=args.visits,
        physicalFilters=args.physicalFilters,
        bands=args.bands,
        ccds=args.ccds,
        ccdKey=args.ccdKey,
        showPatch=args.showPatch,
        showPatchSelectedTractsOnly=args.showPatchSelectedTractsOnly,
        saveFile=args.saveFile,
        showCcds=args.showCcds,
        showCcdsAll=args.showCcdsAll,
        plotFailsOnly=args.plotFailsOnly,
        visitVetoFile=args.visitVetoFile,
        minOverlapFraction=args.minOverlapFraction,
        trimToTracts=args.trimToTracts,
        trimToOverlappingTracts=args.trimToOverlappingTracts,
        doUnscaledLimitRatio=args.doUnscaledLimitRatio,
        forceScaledLimitRatio=args.forceScaledLimitRatio,
        imageDatasetType=args.imageDatasetType,
        visitSummaryDatasetType=args.visitSummaryDatasetType,
        dpi=args.dpi,
        maxVisitForLegend=args.maxVisitForLegend,
    )
