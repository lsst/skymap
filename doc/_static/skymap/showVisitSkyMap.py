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

import sys
import argparse
import logging
from astropy import units
from astropy.coordinates import SkyCoord
import matplotlib
import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import numpy as np

import lsst.afw.cameraGeom as cameraGeom
import lsst.daf.butler as dafButler
import lsst.geom as geom
import lsst.sphgeom as sphgeom

# from lsst.utils.plotting import (
#     get_multiband_plot_colors,
#     get_multiband_plot_symbols,
#     get_multiband_plot_linestyles
# )

from lsst.utils.plotting import publication_plots, make_figure, set_rubin_plotstyle

# publication_plots.set_rubin_plotstyle()

bands_dict = publication_plots.get_band_dicts()
# bandColorDict = bands_dict["colors"]
# Force to unpushed updated color scheme.
bandColorDict ={"u": "#48A8D4",
                "g": "#31DE1F",
                "r": "#B52626",
                "i": "#2915A4",
                "z": "#AD03EA",
                "y": "#2D0201",
}
bandSymbolDict = bands_dict["symbols"]
bandLinestyleDict = bands_dict["line_styles"]
set_rubin_plotstyle()

logger = logging.getLogger("lsst.skymap.bin.showVisitSkyMap")


def bboxToRaDec(bbox, wcs):
    """Get the corners of a BBox and convert them to lists of RA and Dec."""
    sphPoints = wcs.pixelToSky(geom.Box2D(bbox).getCorners())
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
    return m + percentile*interval


def get_cmap(n, name="hsv"):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.
    """
    return matplotlib.colormaps[name].resampled(n)


def main(repo, collections, skymapName=None, tracts=None, visits=None, patches=None,
         physicalFilters=None, bands=None, ccds=None, ccdKey="detector", showPatch=False,
         saveFile=None, showCcds=False, showCcdsAll=True, visitVetoFile=None,
         minOverlapFraction=None, trimToTracts=False, doUnscaledLimitRatio=False,
         forceScaledLimitRatio=False, plotFailsOnly=False):
    if minOverlapFraction is not None and tracts is None:
        raise RuntimeError("Must specify --tracts if --minOverlapFraction is set")
    logger.info("Making butler for collections = %s in repo %s", collections, repo)
    butler = dafButler.Butler(repo, collections=collections)
    instrument = butler.find_dataset("camera").dataId["instrument"]
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
            skymapName = "lsst_cells_v1"
        else:
            logger.error("Unknown skymapName for instrument: %s.  Must specify --skymapName on command line.",
                         instrument)
    logger.info("instrument = %s skymapName = %s", instrument, skymapName)
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

    if "Trifid-Lagoon" in collections[0]:
        whereStr += " AND detector NOT IN (120,121,122,78,79,80,0,20,27,65,123,161,168,188,158,68,169,187,19) AND detector<189 AND exposure NOT IN (2025050400538)"
        detectorSkipList = [120,121,122,78,79,80,0,20,27,65,
                            123,161,168,188,158,68,169,187,19]
        # whereStr += " AND detector NOT IN (0,20,27,65,123,161,168,188,78,79,80,121,122,169,187,120,158,30,68,1,19,170,186,117,155,33,71,2,18,165,185,124,160,28,64,3,23,125,118,119) AND detector<189 AND exposure NOT IN (2025050400538)"
        # detectorSkipList = [0,20,27,65,123,161,168,188,78,79,80,
        #                     121,122,169,187,120,158,30,68,
        #                     1,19,170,186,117,155,33,71,
        #                     2,18,165,185,124,160,28,64,
        #                     3,23,125,118,119]

    elif "202508" in collections[0]:
        whereStr += " AND detector NOT IN (0,20,65,161,188,168,123,27,1,19,68,158,187,169,120,30) AND detector<189"
        detectorSkipList = [0,20,65,161,188,168,123,27,1,19,68,158,187,169,120,30]
    elif instrument == "LSSTCam":
        whereStr += " AND detector NOT IN (120,121,122,78,79,80,0,20,27,65,123,161,168,188) AND detector<189"
        detectorSkipList = [120,121,122,78,79,80,0,20,27,65,123,161,168,188]

    if len(whereStr) > 1:
        whereStr = "instrument=\'" + instrument + "\' AND skymap=\'" + skymapName + "\' AND " + whereStr
    else:
        whereStr = "instrument=\'" + instrument + "\' AND skymap=\'" + skymapName + "\'"

    logger.info("Applying the following where clause in dataId search: %s", whereStr)

    getIdsFromIsr = False  # True
    isrVisits = []
    isrDetectors = None
    if getIdsFromIsr:
        isrDataRefs = list(butler.registry.queryDatasets("post_isr_image", where=whereStr).expanded())
        print("Number of post_isr_image datasets in {}: {}".format(collections[0], len(isrDataRefs)))
        isrDataIds = []
        isrDetectors = []
        for isrDataRef in isrDataRefs:
            exposure = isrDataRef.dataId["exposure"]
            if exposure not in isrVisits:
                isrVisits.append(exposure)
            detector = isrDataRef.dataId["detector"]
            if detector not in isrDetectors:
                isrDetectors.append(detector)
            isrDataIds.append((exposure, detector))

    if len(isrVisits) > 0 :
        isrVisitStr = makeWhereInStr("exposure", isrVisits, int)
        if len(whereStr) < 1:
            whereStr += isrVisitStr
        else:
            whereStr += " AND " + isrVisitStr
        visits = isrVisits

    visitVetoList = []
    if visitVetoFile is not None:
        with open(visitVetoFile) as f:
            content = f.readlines()
            visitVetoList = [int(visit.strip()) for visit in content]

    # imageTypes = ["science", "acq"]
    # imageTypes = ["science"]
    imageTypes = None
    if imageTypes is not None:
        imageTypeStr = makeWhereInStr("exposure.observation_type", imageTypes, str)
        rawWhereStr = whereStr + " AND " + imageTypeStr if len(whereStr) else " " + imageTypeStr
    else:
        rawWhereStr = whereStr
    # rawWhereStr += " AND exposure.exposure_time>25.0"
    if len(isrVisits) == 0:
        rawDataRefs = list(butler.registry.queryDatasets("raw", where=rawWhereStr).expanded())
        logger.info("Limiting search to query: %s", rawWhereStr)
        rawVisits = []
        for rawDataRef in rawDataRefs:
            rawVisit = rawDataRef.dataId.exposure.id
            if rawVisit not in rawVisits and rawVisit not in visitVetoList:
                rawVisits.append(rawVisit)
    else:
        rawVisits = isrVisits
        rawDataRefs = isrDataRefs

    rawVisits.sort()
    logger.info("Number of raw visits satisfying where and veto clauses: %d (nDatId: %d)",
                len(rawVisits), len(rawDataRefs))

    # Guess at data type (OG calexp vs. new initial_pvi vs. newer preliminary_visit_image)
    datasetTypesToTry = ["preliminary_visit_image", "calexp", "initial_pvi"]
    for datasetType in datasetTypesToTry:
        dataRefs = list(butler.registry.queryDatasets(datasetType, where=whereStr).expanded())
        logger.info("Number of %s dataset = %d", datasetType, len(dataRefs))
        if len(dataRefs) > 0:
            break
    if len(dataRefs) == 0 and tracts == None:
        raise RuntimeError("No visits of any of {} found.".format(datasetTypesToTry))
    logger.info("Determined that the visit-level image datasetType is %s", datasetType)
    if visits is None:
        visits = []
        for dataRef in dataRefs:
            visit = dataRef.dataId.visit.id
            if visit not in visits and visit not in visitVetoList and visit in rawVisits:
                visits.append(visit)
        visits.sort()
        logger.info(
            "List of visits (N=%d) satisfying where and veto clauses: %s", len(visits), visits
        )
    else:
        if len(visitVetoList) > 1:
            visitListTemp = visits.copy()
            for visit in visitListTemp:
                if visit in visitVetoList or visit not in rawVisits:
                    visits.remove(visit)
            logger.info("List of visits (N=%d) excluding veto list: %s", len(visits), visits)
        logger.info("List of visits (N=%d): %s", len(visits), visits)

    if datasetType == "preliminary_visit_image":
        visitSummaryStr = "preliminary_visit_summary"  # "visit_summary"
    else:
        visitSummaryStr = "visitSummary"

    if isrDetectors is None:
        ccdIdList = []
        for ccd in camera:
            ccdId = ccd.getId()
            if ((ccds is None or ccdId in ccds) and ccd.getType() == cameraGeom.DetectorType.SCIENCE
                and ccdId not in detectorSkipList):
                ccdIdList.append(ccdId)
        ccdIdList.sort()
    else:
        ccdIdList = isrDetectors
    print("ccdIdList: {}".format(ccdIdList))
    print("detectorSkipList: {}".format(detectorSkipList))
    nDetTot = len(ccdIdList)

    visitIncludeList = []
    # Determine the fraction of detectors that overlap any tract under
    # consideration. If this fraction does not exceed minOverlapFraction,
    # skip this visit.
    if minOverlapFraction is not None:
        nDataIdsFound = 0
        for i_v, visit in enumerate(visits):
            ccdOverlapList = []
            try:
                visitSummary = butler.get(visitSummaryStr, visit=visit)
            except LookupError as e:
                logger.warning("%s  Will try to get wcs from exposure.", e)
                visitSummary = None
            if tracts is not None:
                for tract in tracts:
                    tractInfo = skymap[tract]
                    sphCorners = tractInfo.wcs.pixelToSky(geom.Box2D(tractInfo.bbox).getCorners())
                    tractConvexHull = sphgeom.ConvexPolygon.convexHull(
                        [coord.getVector() for coord in sphCorners])
                    for ccdId in ccdIdList:
                        if ccdId not in ccdOverlapList:
                            raCorners, decCorners = getDetRaDecCorners(
                                ccdKey, ccdId, datasetType, visit, visitSummary=visitSummary,
                                butler=butler, doLogWarn=False)
                            if raCorners is not None and decCorners is not None:
                                nDataIdsFound += 1
                                detSphCorners = []
                                for ra, dec in zip(raCorners, decCorners):
                                    pt = geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                                          geom.Angle(dec, geom.degrees))
                                    detSphCorners.append(pt)
                                detConvexHull = sphgeom.ConvexPolygon.convexHull(
                                    [coord.getVector() for coord in detSphCorners])
                                if tractConvexHull.contains(detConvexHull):
                                    ccdOverlapList.append(ccdId)

                    if len(ccdOverlapList)/nDetTot >= minOverlapFraction:
                        break
                if len(ccdOverlapList)/nDetTot < minOverlapFraction:
                    logger.info("Fraction of detectors overlaping any tract for visit %d (%.2f) < "
                                "minimum required (%.2f).  Skipping visit...",
                                visit, len(ccdOverlapList)/nDetTot, minOverlapFraction)
                else:
                    if visit not in visitIncludeList:
                        visitIncludeList.append(visit)
    else:
        visitIncludeList = visits

    print("len(visitIncludeList) = {}".format(len(visitIncludeList)))
    # for finalVisit in visitIncludeList:
    #    print("LSSTCam {}".format(finalVisit))
    # raise("Remove me....")

    # Draw the CCDs.
    ras, decs = [], []
    rasRaw, decsRaw = [], []
    bboxesPlotted = []
    cmap = get_cmap(len(visitIncludeList))
    alphaEdge = 0.7
    maxVisitForLegend = 150  # 25  # 46  # 20
    finalVisitList = []
    includedBands = []
    includedPhysicalFilters = []
    nDataIdFound = 0
    nRawDataId = 0
    rawDataIdList = []
    dataIdFoundList = []
    noRawDataIdList = []
    failedDataIdList = []
    removededDataIdList = []
    failedVisitList = []
    for i_v, visit in enumerate(visitIncludeList):
        print("Working on visit %d [%d of %d]" % (visit, i_v + 1, len(visitIncludeList)), end="\r")
        try:
            visitSummary = butler.get(visitSummaryStr, visit=visit)
        except Exception as e:
            logger.warning("%s  Will try to get wcs from exposure.", e)
            visitSummary = None

        band, physicalFilter = getBand(visitSummary=visitSummary, butler=butler, visit=visit)
        if band not in includedBands:
            includedBands.append(band)
        if physicalFilter not in includedPhysicalFilters:
            includedPhysicalFilters.append(physicalFilter)
        # color = cmap(i_v)
        color = bandColorDict[band]
        linestyle = bandLinestyleDict[band]
        fillKwargs = {"fill": False, "alpha": alphaEdge, "facecolor": None, "edgecolor": color, "lw": 0.6}
        inLegend = False
        iDataId = 0
        halfWay = int(maxVisitForLegend/2) - 1

        doPlotRawOutlines = False  # True
        for ccdId in ccdIdList:
            if getIdsFromIsr:
                if (visit, ccdId) not in isrDataIds:
                    # logger.info("Skipping {} {}".format(visit, ccdId))
                    continue

            # Plot raw WCSs
            raCornersRaw, decCornersRaw = getDetRaDecCorners(
                ccdKey, ccdId, "raw", visit, visitSummary=visitSummary, butler=butler)
            if raCornersRaw is not None and decCornersRaw is not None:
                nRawDataId += 1
                rawDataIdList.append({"visit": visit, "detector": ccdId})
                rasRaw += raCornersRaw
                decsRaw += decCornersRaw
                if doPlotRawOutlines:
                    if not plotFailsOnly:
                        plt.fill(raCornersRaw, decCornersRaw, fill=False, alpha=alphaEdge, color=color,
                                 edgecolor=color, ls=":", lw=0.8)
                # Always plot gray outline of raws when plotting failures only.
                if plotFailsOnly:
                    plt.fill(raCornersRaw, decCornersRaw, fill=True, alpha=0.4, color="lightgray",
                             edgecolor="lightgray", ls=":", lw=0.8)
                    plt.fill(raCornersRaw, decCornersRaw, fill=False, alpha=1.0, color="lightgray",
                             edgecolor="lightgray", ls=":", lw=0.8)
                if visit not in finalVisitList:
                    finalVisitList.append(visit)

            # dataIdExists = butler.exists(datasetType, visit=visit, detector=ccdId)
            dataIdFailed = False
            row = []
            if visitSummary is not None:
                row = visitSummary[visitSummary["id"] == ccdId]
                if len(row) > 0:
                    ra = row["ra"][0]
                    dec = row["dec"][0]
                    psfSigma = row["psfSigma"][0]
                    zeroPoint = row["zeroPoint"][0]
                    effTime = row["effTime"][0]
                    if (~np.isfinite(ra) or ~np.isfinite(dec) or ~np.isfinite(psfSigma)
                        or ~np.isfinite(zeroPoint) or ~np.isfinite(effTime)):
                        dataIdFailed = True
                else:
                    dataIdFailed = True
            else:
                try:
                    exp = butler.get(datasetType, visit=visit, detector=ccdId)
                    wcs = exp.getWcs()
                    photoCalib = exp.getPhotoCalib()
                    if wcs is None or photoCalib is None:
                        dataIdFailed = True
                except Exception as e:
                    print(e)
                    dataIdFailed = True

            if not dataIdFailed and not plotFailsOnly:
                iDataId += 1
                raCorners, decCorners = getDetRaDecCorners(
                    ccdKey, ccdId, datasetType, visit, visitSummary=visitSummary, butler=butler)
                if raCorners is not None and decCorners is not None:
                    nDataIdFound += 1
                    dataIdFoundList.append({"visit": visit, "detector": ccdId})
                    ras += raCorners
                    decs += decCorners
                    if not inLegend:
                        if i_v < halfWay or i_v > len(visitIncludeList) - halfWay:
                            plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4,
                                     color=color,
                                     edgecolor=color, ls=linestyle, lw=0.8, label="{} {}".
                                     format(str(visit), band))
                        elif i_v == halfWay:
                            plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4, color=color,
                                     edgecolor=color, ls=linestyle, lw=0.8, label="[...]")
                        else:
                            plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4, color=color,
                                     edgecolor=color, ls=linestyle, lw=0.8)
                    else:
                        plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4, color=color,
                                 edgecolor=color, ls=linestyle, lw=0.8)
                    inLegend = True

                    # if not inLegend and len(visitIncludeList) <= maxVisitForLegend:
                    #     plt.fill(raCorners, decCorners, label="{} {}".format(str(visit), band),
                    #              **fillKwargs)
                    #     inLegend = True
                    # else:
                    #     plt.fill(raCorners, decCorners, **fillKwargs)
                    # plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4, color=color,
                    #          edgecolor=color)
                    if visit not in finalVisitList:
                        finalVisitList.append(visit)
                    # add CCD serial numbers
                    if showCcds or showCcdsAll and not plotFailsOnly:
                        overlapFrac = 0.2
                        deltaRa = max(raCorners) - min(raCorners)
                        deltaDec = max(decCorners) - min(decCorners)
                        minPoint = geom.Point2D(min(raCorners) + overlapFrac*deltaRa,
                                                min(decCorners) + overlapFrac*deltaDec)
                        maxPoint = geom.Point2D(max(raCorners) - overlapFrac*deltaRa,
                                                max(decCorners) - overlapFrac*deltaDec)
                        # Use doubles in Box2D to check overlap
                        bboxDouble = geom.Box2D(minPoint, maxPoint)
                        overlaps = [not bboxDouble.overlaps(otherBbox) for otherBbox in bboxesPlotted]
                        if showCcdsAll:
                            plt.text(getValueAtPercentile(raCorners), getValueAtPercentile(decCorners),
                                     str(ccdId), fontsize=6, ha="center", va="center",
                                     color=color, alpha=1.0)
                            bboxesPlotted.append(bboxDouble)
                        else:
                            if all(overlaps):
                                plt.text(getValueAtPercentile(raCorners), getValueAtPercentile(decCorners),
                                         str(ccdId), fontsize=6, ha="center", va="center",
                                         color=color, alpha=1.0)
                                # str(ccdId), fontsize=6, ha="center", va="center", color="darkblue")
                                bboxesPlotted.append(bboxDouble)
                else:
                    # if visit not in noRawVisitList:
                    #     noRawVisitList.append(visit)
                    noRawDataIdList.append({"visit": visit, "detector": ccdId})
            else:
                # raCorners, decCorners = getDetRaDecCorners(
                #     ccdKey, ccdId, datasetType, visit, visitSummary=visitSummary, butler=butler)
                # if raCorners == None or decCorners == None:
                if dataIdFailed:
                    if visit not in failedVisitList:
                        failedVisitList.append(visit)
                    failedDataIdList.append({"visit": visit, "detector": ccdId})

    logger.info("Final list of visits (nVisit=%d nDataId=%d) satisfying where and minOverlapFraction "
                "clauses: %s", len(finalVisitList), nDataIdFound, finalVisitList)

    if plotFailsOnly:
        nNoIcExp = 0
        maxFailedVisitForLegend = 25
        logger.warning("Only plotting failed detectors...")
        visitSummary = None
        rasFailed, decsFailed = [], []
        # Plot raw WCS of the failed (i.e. raw exists but calexp does not) dataIds only.
        # for rawDataId in rawDataIdList:
        for iFailed, failedDataId in enumerate(failedDataIdList):
            visit = failedDataId["visit"]
            ccdId = failedDataId["detector"]
            band, physicalFilter = getBand(visitSummary=None, butler=butler, visit=visit)
            # color = "teal"  # cmap(i_v)
            color = bandColorDict[band]
            # marker = bandSymbolDict[band]
            linestyle = bandLinestyleDict[band]
            print("Working on failed visit %d detector %d band %s [%d of %d]" % (
                visit, ccdId, band, iFailed + 1, len(failedDataIdList)), end="\r")
            # Plot raw WCSs
            raCornersFailed, decCornersFailed = getDetRaDecCorners(
                ccdKey, ccdId, "raw", visit, butler=butler)
                # ccdKey, ccdId, "raw", visit, visitSummary=visitSummary, butler=butler)
            # print(visit, ccdId, raCornersFailed, decCornersFailed)
            if raCornersFailed is not None and decCornersFailed is not None:
                rasFailed += raCornersFailed
                decsFailed += decCornersFailed
                # See if icExp exists to see if it failed there, noting that it
                # may actually have been removed/pruned due to PSF funkyness
                # rather than a true calibration fail.
                icExpExists = butler.exists("icExp", visit=visit, detector=ccdId)
                if icExpExists:
                    halfWay = int(maxFailedVisitForLegend/2) - 1
                    if iFailed < halfWay or iFailed > len(failedDataIdList) - halfWay:
                            plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                     color=color, edgecolor=color, ls=linestyle, lw=0.8, label="{} {} {}".
                                     format(str(visit), ccdId, band))
                    elif iFailed == halfWay:
                        plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                 color=color, edgecolor=color, ls=linestyle, lw=0.8, label="[...]")
                    else:
                        plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                 color=color, edgecolor=color, ls=linestyle, lw=0.8)
                else:
                    nNoIcExp += 1
                    if iFailed < maxFailedVisitForLegend:
                        plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                 color=color, edgecolor="darkblue", lw=1.2, zorder=100, label="{} {} {}".
                                 format(str(visit), ccdId, band))
                    elif iFailed == maxFailedVisitForLegend:
                        plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                 color=color, edgecolor="darkblue", lw=1.2, zorder=100, label="etc...")
                    else:
                        plt.fill(raCornersFailed, decCornersFailed, fill=True, alpha=alphaEdge/4,
                                 color=color, edgecolor="darkblue", lw=1.2, zorder=100)

            # add CCD serial numbers
            if showCcds or showCcdsAll:
                overlapFrac = 0.2
                deltaRa = max(raCornersFailed) - min(raCornersFailed)
                deltaDec = max(decCornersFailed) - min(decCornersFailed)
                minPoint = geom.Point2D(min(raCornersFailed) + overlapFrac*deltaRa,
                                        min(decCornersFailed) + overlapFrac*deltaDec)
                maxPoint = geom.Point2D(max(raCornersFailed) - overlapFrac*deltaRa,
                                        max(decCornersFailed) - overlapFrac*deltaDec)
                # Use doubles in Box2D to check overlap
                bboxDouble = geom.Box2D(minPoint, maxPoint)
                overlaps = [not bboxDouble.overlaps(otherBbox) for otherBbox in bboxesPlotted]
                if showCcdsAll:
                    plt.text(getValueAtPercentile(raCornersFailed),
                             getValueAtPercentile(decCornersFailed),
                             str(ccdId), fontsize=6, ha="center", va="center", color=color, alpha=1.0)
                    bboxesPlotted.append(bboxDouble)
                else:
                    if all(overlaps):
                        plt.text(getValueAtPercentile(raCornersFailed),
                                 getValueAtPercentile(decCornersFailed),
                                 str(ccdId), fontsize=6, ha="center", va="center", color=color, alpha=1.0)
                        # str(ccdId), fontsize=6, ha="center", va="center", color="darkblue")
                        bboxesPlotted.append(bboxDouble)
        finalVisitList = failedVisitList
        ras = rasRaw  # rasFailed
        decs = decsRaw  # decsFailed

        if len(failedDataIdList)%5 == 0 and len(ras) > 0:
            logger.info("Saving interim plot at N = {} detectors added.".format(len(failedDataIdList)))
            savePlot(repo, butler, collections, skymap, camera, tracts, visits, ras, decs, plotFailsOnly,
                     finalVisitList, nDataIdFound, nRawDataId, ccdIdList, ccdKey, datasetType,
                     failedDataIdList, minOverlapFraction, includedBands, includedPhysicalFilters,
                     trimToTracts, showPatch, forceScaledLimitRatio, doUnscaledLimitRatio, visitSummaryStr,
                     visitIncludeList, cmap, alphaEdge, maxVisitForLegend, saveFile)

    savePlot(repo, butler, collections, skymap, camera, tracts, visits, ras, decs, plotFailsOnly,
             finalVisitList, nDataIdFound, nRawDataId, ccdIdList, ccdKey, datasetType,
             failedDataIdList, minOverlapFraction, includedBands, includedPhysicalFilters,
             trimToTracts, showPatch, forceScaledLimitRatio, doUnscaledLimitRatio, visitSummaryStr,
             visitIncludeList, cmap, alphaEdge, maxVisitForLegend, saveFile)


def savePlot(repo, butler, collections, skymap, camera, tracts, visits, ras, decs, plotFailsOnly,
             finalVisitList, nDataIdFound, nRawDataId, ccdIdList, ccdKey, datasetType,
             failedDataIdList, minOverlapFraction, includedBands, includedPhysicalFilters,
             trimToTracts, showPatch, forceScaledLimitRatio, doUnscaledLimitRatio, visitSummaryStr,
             visitIncludeList, cmap, alphaEdge, maxVisitForLegend, saveFile):

    raToDecLimitRatio = None
    if len(ras) > 0:
        tractList = list(set(skymap.findTractIdArray(ras, decs, degrees=True)))
        minVisitRa, maxVisitRa = min(ras), max(ras)
        minVisitDec, maxVisitDec = min(decs), max(decs)
        raVisitDiff = maxVisitRa - minVisitRa
        decVisitDiff = maxVisitDec - minVisitDec
        midVisitRa = minVisitRa + 0.5*raVisitDiff
        midVisitDec = minVisitDec + 0.5*decVisitDiff
        midRa = np.atleast_1d((midVisitRa*units.deg).to(units.radian).value).astype(np.float64)
        midDec = np.atleast_1d((midVisitDec*units.deg).to(units.radian).value).astype(np.float64)
        midSkyCoord = SkyCoord(midVisitRa*units.deg, midVisitDec*units.deg)
    else:
        if tracts is not None:
            logger.info("No exposures were found, but --tracts list was provided, so will go ahead and "
                        "plot the empty tracts.")
            tractList = tracts
            trimToTracts = True
        else:
            logger.warning("No data to plot (if you want to plot empty tracts, include them as "
                           "a blank-space separated list to the --tracts option.")
            return None
            # raise RuntimeError("No data to plot (if you want to plot empty tracts, include them as "
            #                    "a blank-space separated list to the --tracts option.")
    tractList.sort()
    logger.info("List of tracts overlapping data:  %s", tractList)
    tractLimitsDict = getTractLimitsDict(skymap, tractList)

    if forceScaledLimitRatio:
        doUnscaledLimitRatio = False
    else:
        # Roughly compute radius in degrees of a single detector.  If RA/Dec
        # coverage is more than 30 times the detector radius, and the RA/Dec
        # limit ratio is greater than raDecScaleThresh, don't try to scale to
        # detector coords.
        radiusMm = camera.computeMaxFocalPlaneRadius()
        fpRadiusPt = geom.Point2D(radiusMm, radiusMm)
        focalPlaneToFieldAngle = camera.getTransformMap().getTransform(
            cameraGeom.FOCAL_PLANE, cameraGeom.FIELD_ANGLE
        )
        fpRadiusDeg = np.rad2deg(focalPlaneToFieldAngle.applyForward(fpRadiusPt))[0]
        detectorRadiusDeg = fpRadiusDeg/np.sqrt(len(camera))

        if trimToTracts:
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
                (xDelta0/yDelta0 > raDecScaleThresh or yDelta0/xDelta0 > raDecScaleThresh)
                and max(xDelta0, yDelta0) > 70*detectorRadiusDeg
                and yLimMin < 75.0 and yLimMax > -75.0
        ):
            logger.info(
                "Sky coverage is large (and not too close to a pole), so not scaling to detector coords."
            )
            doUnscaledLimitRatio = True

    if not doUnscaledLimitRatio:
        # Find a detector that contains the mid point in RA/Dec (or the closest
        # one) to set the plot aspect ratio.
        minDistToMidCoord = 1e12
        minSepVisit = None
        minSepCcdId = None
        for i_v, visit in enumerate(visits):
            try:
                visitSummary = butler.get(visitSummaryStr, visit=visit)
            except Exception as e:
                logger.warning("%s  Will try to get wcs from exposure.", e)
                visitSummary = None
            for ccdId in ccdIdList:
                raCorners, decCorners = getDetRaDecCorners(
                    ccdKey, ccdId, datasetType, visit, visitSummary=visitSummary,
                    butler=butler, doLogWarn=False
                )
                if raCorners is not None and decCorners is not None:
                    detSphCorners = []
                    for ra, dec in zip(raCorners, decCorners):
                        pt = geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                              geom.Angle(dec, geom.degrees))
                        detSphCorners.append(pt)
                        ptSkyCoord = SkyCoord(ra*units.deg, dec*units.deg)
                        separation = (midSkyCoord.separation(ptSkyCoord)).degree
                        if separation < minDistToMidCoord:
                            minSepVisit = visit
                            minSepCcdId = ccdId
                            minDistToMidCoord = separation
                    detConvexHull = sphgeom.ConvexPolygon(
                        [coord.getVector() for coord in detSphCorners])
                    if detConvexHull.contains(midRa, midDec) and raToDecLimitRatio is None:
                        logger.info(
                            "visit/det overlapping plot coord mid point in RA/Dec: %d %d", visit, ccdId
                        )
                        raToDecLimitRatio = (
                            (max(raCorners) - min(raCorners))/(max(decCorners) - min(decCorners))
                        )
                        det = camera[ccdId]
                        width = det.getBBox().getWidth()
                        height = det.getBBox().getHeight()
                        if raToDecLimitRatio > 1.0:
                            raToDecLimitRatio /= max(height/width, width/height)
                        else:
                            if raToDecLimitRatio < 1.0:
                                raToDecLimitRatio *= max(height/width, width/height)
                        break
            if raToDecLimitRatio is not None:
                break

        if raToDecLimitRatio is None and minSepVisit is not None:
            try:
                visitSummary = butler.get(visitSummaryStr, visit=minSepVisit)
            except Exception as e:
                logger.warning("%s  Will try to get wcs from exposure.", e)
                visitSummary = None
            raCorners, decCorners = getDetRaDecCorners(
                ccdKey, minSepCcdId, datasetType, minSepVisit, visitSummary=visitSummary,
                butler=butler, doLogWarn=False
            )
            for ra, dec in zip(raCorners, decCorners):
                pt = geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                      geom.Angle(dec, geom.degrees))
                detSphCorners.append(pt)
            detConvexHull = sphgeom.ConvexPolygon([coord.getVector() for coord in detSphCorners])
            logger.info(
                "visit/det closest to plot coord mid point in RA/Dec (none actually overlap it): %d %d",
                minSepVisit, minSepCcdId
            )
            raToDecLimitRatio = (max(raCorners) - min(raCorners))/(max(decCorners) - min(decCorners))
            det = camera[minSepCcdId]
            width = det.getBBox().getWidth()
            height = det.getBBox().getHeight()
            if raToDecLimitRatio > 1.0:
                raToDecLimitRatio /= max(height/width, width/height)
            else:
                if raToDecLimitRatio < 1.0:
                    raToDecLimitRatio *= max(height/width, width/height)

    if trimToTracts is True:
        xlim, ylim = derivePlotLimits(tractLimitsDict, raToDecLimitRatio=raToDecLimitRatio, buffFrac=0.04)
    else:
        visitLimitsDict = {"allVisits": {"ras": [minVisitRa, maxVisitRa],
                                         "decs": [minVisitDec, maxVisitDec]}}
        xlim, ylim = derivePlotLimits(visitLimitsDict, raToDecLimitRatio=raToDecLimitRatio, buffFrac=0.04)
    # xlim = 167.5, 163.5
    # ylim = -61.5, -59.5

    if doUnscaledLimitRatio:
        boxAspectRatio = abs((ylim[1] - ylim[0])/(xlim[1] - xlim[0]))
    else:
        boxAspectRatio = 1.0

    # Draw the skymap.
    alpha0 = 1.0
    tractHandleList = []
    tractStrList = []
    if tracts is not None:
        tractOutlineList = list(set(tracts + tractList))
    else:
        tractOutlineList = tractList
    tractOutlineList.sort()
    logger.info("List of tract outlines that would get plotted...: %s", tractOutlineList)
    plotTractOutlines = True  # False
    if not plotTractOutlines:
        tractOutlineList = []
        logger.info("...but not plotting them on this run!")
    for i_t, tract in enumerate(tractOutlineList):
        alpha = max(0.1, alpha0 - i_t*1.0/len(tractOutlineList))
        tractInfo = skymap[tract]
        tCenter = tractInfo.ctr_coord
        tCenterRa = tCenter.getRa().asDegrees()
        tCenterDec = tCenter.getDec().asDegrees()
        fracDeltaX = 0.02*abs((xlim[1] - xlim[0]))
        fracDeltaY = 0.02*abs((ylim[1] - ylim[0]))
        if (xlim[1] + fracDeltaX < tCenterRa < xlim[0] - fracDeltaX
                and ylim[0] + fracDeltaY < tCenterDec < ylim[1] - fracDeltaY):
            if len(tractOutlineList) > 1 or not showPatch:
                if not showPatch:
                    plt.text(tCenterRa, tCenterDec, tract, fontsize=7, alpha=alpha, ha="center", va="center")
                else:
                    plt.text(tCenterRa, tCenterDec, tract, fontsize=7, alpha=1, color="white",
                             path_effects=[pathEffects.withStroke(linewidth=3, foreground="black")],
                             fontweight=500, ha="center", va="center", zorder=5)
        ra, dec = bboxToRaDec(tractInfo.bbox, tractInfo.getWcs())
        plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed", alpha=alpha)
        tractArtist = matplotlib.patches.Patch(fill=False, edgecolor="k", linestyle="dashed", alpha=alpha)
        tractHandleList.append(tractArtist)
        tractStrList.append(str(tract))
        if showPatch:
            patchColor = "k"
            for patch in tractInfo:
                patchFontSize = 5 if patch.sequential_index < 1000 else 4
                ra, dec = bboxToRaDec(patch.getInnerBBox(), tractInfo.getWcs())
                plt.fill(ra, dec, fill=False, edgecolor=patchColor, lw=0.5, linestyle=(0, (5, 6)),
                         alpha=alpha)
                if (xlim[1] + fracDeltaX < getValueAtPercentile(ra) < xlim[0] - fracDeltaX
                        and ylim[0] + fracDeltaY < getValueAtPercentile(dec) < ylim[1] - fracDeltaY):
                    plt.text(getValueAtPercentile(ra), getValueAtPercentile(dec),
                             str(patch.sequential_index), fontsize=patchFontSize, color=patchColor,
                             ha="center", va="center", alpha=alpha)

    # Add labels and save.
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_box_aspect(boxAspectRatio)
    if abs(xlim[1] > 99.99):
        ax.tick_params("x", labelrotation=45, pad=5, labelsize=8)
    else:
        ax.tick_params("x", labelrotation=0, pad=5, labelsize=8)
    ax.tick_params("y", labelsize=8)
    ax.set_xlabel("R.A. (deg)", fontsize=9)
    ax.set_ylabel("Dec. (deg)", fontsize=9)

    visitScaleOffset = None
    # if 0:
    if len(visitIncludeList) > maxVisitForLegend and not plotFailsOnly:
        nz = matplotlib.colors.Normalize()
        colorBarScale = finalVisitList
        if max(finalVisitList) > 9999999:
            visitScaleOffset = min(finalVisitList)
            colorBarScale = [visit - visitScaleOffset for visit in finalVisitList]
        nz.autoscale(colorBarScale)
        cax, _ = matplotlib.colorbar.make_axes(plt.gca(), pad=0.03)
        cax.tick_params(labelsize=7)
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=nz, alpha=alphaEdge,
                                              format=lambda x, _: f"{x:.0f}")
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
        colorBarLabel = "visit number"
        if visitScaleOffset is not None:
            colorBarLabel += " - {:d}".format(visitScaleOffset)
        cb.set_label(colorBarLabel, rotation=-90, labelpad=13, fontsize=9)
        tractLegend = Legend(ax, tractHandleList, tractStrList, loc="upper right", fancybox=True,
                             shadow=True, fontsize=5, title_fontsize=6, title="tracts")
        ax.add_artist(tractLegend)
    else:
        doPlotLegend = True  # False
        if doPlotLegend:
            if len(visitIncludeList) > 0:
                labelStr = "visits" if not plotFailsOnly else "failed visits"
                xBboxAnchor = min(1.25, max(1.03, boxAspectRatio*1.15))
                ax.legend(loc="center left", bbox_to_anchor=(xBboxAnchor, 0.5), fancybox=True,
                          shadow=True, fontsize=6, title_fontsize=6, title=labelStr)
            # Create the second legend and add the artist manually.
            tractLegend = Legend(ax, tractHandleList, tractStrList, loc="center left",
                                 bbox_to_anchor=(1.0, 0.5),
                                 fancybox=True, shadow=True, fontsize=6, title_fontsize=6, title="tracts")
            ax.add_artist(tractLegend)

    titleStr = repo + "\n" + collections[0]
    if len(collections) > 1:
        for collection in collections[1:]:
            titleStr += "\n" + collection
    if not plotFailsOnly:
        titleStr += "\nnVisit: {} nDataId: {} (nRaw: {})".format(
            str(len(finalVisitList)), str(nDataIdFound), nRawDataId)
    else:
        # titleStr += "\nNO icExp: N={} NO calexp: N={}".format(nNoIcExp, len(failedDataIdList) - nNoIcExp)
        titleStr += "\nFailed calibrateImage: N={}".format(len(failedDataIdList))
    logger.info("DataIds that failed calibrateImage (N=%d): {%s}", len(failedDataIdList), failedDataIdList)
    if minOverlapFraction is not None:
        titleStr += " (minOvlpFrac = {:.2f})".format(minOverlapFraction)
    if len(includedBands) > 0:
        titleStr += "\nbands: {}".format(str(includedBands).translate({ord(i): None for i in "[]\'"}))
    if len(includedPhysicalFilters) > 0:
        if len(includedPhysicalFilters[0]) > 9:
            titleStr += "\n"
        titleStr += "  physical filters: {}".format(str(includedPhysicalFilters).translate(
            {ord(i): None for i in "[]\'"}))

    doPlotTitle = True  # False
    if doPlotTitle:
        ax.set_title("{}".format(titleStr), fontsize=8)

    fig = plt.gcf()
    # fig = make_figure()

    if boxAspectRatio > 1.0:
        minInches = max(4.0, 0.3*abs(xlim[1] - xlim[0]))
        xInches = minInches
        yInches = min(120.0, boxAspectRatio*minInches)
        fig.set_size_inches(xInches, yInches)
    if boxAspectRatio < 1.0:
        minInches = max(4.0, 0.3*abs(ylim[1] - ylim[0]))
        xInches = min(120.0, minInches/boxAspectRatio)
        yInches = minInches
        fig.set_size_inches(xInches, yInches)
    if saveFile is not None:
        if plotFailsOnly:
            if "fail" not in saveFile:
                saveFile = saveFile[0:-4] + "_failed.png"
        logger.info("Saving file in: %s", saveFile)
        fig.savefig(saveFile, bbox_inches="tight", dpi=150)
    else:
        fig.show()


def makeWhereInStr(parameterName, parameterList, parameterType):
    """Create the string to be used in the where clause for registry lookup.
    """
    typeStr = "\'" if parameterType is str else ""
    whereInStr = parameterName + " IN (" + typeStr + str(parameterList[0]) + typeStr
    if len(parameterList) > 1:
        for param in parameterList[1:]:
            whereInStr += ", " + typeStr + str(param) + typeStr
    whereInStr += ")"

    return whereInStr


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
        xlim = xLimMax + padFrac*xDelta0, xLimMin - padFrac*xDelta0
        ylim = yLimMin - padFrac*yDelta0, yLimMax + padFrac*yDelta0
        return xlim, ylim

    if raToDecLimitRatio == 1.0:
        if xDelta0 > yDelta0:
            xLimMin -= buffFrac*yDelta0
            xLimMax += buffFrac*yDelta0
        else:
            yLimMin -= buffFrac*yDelta0
            yLimMax += buffFrac*yDelta0
        xLimMin, xLimMax, yLimMin, yLimMax = setLimitsToEqualRatio(xLimMin, xLimMax, yLimMin, yLimMax)
    else:
        xLimMin -= buffFrac*xDelta0
        xLimMax += buffFrac*xDelta0
        yLimMin -= buffFrac*yDelta0
        yLimMax += buffFrac*yDelta0
        xLimMin, xLimMax, yLimMin, yLimMax = setLimitsToEqualRatio(xLimMin, xLimMax, yLimMin, yLimMax)
        xDelta = xLimMax - xLimMin
        yDelta = yLimMax - yLimMin
        if raToDecLimitRatio > 1.0:
            if yDelta0 > xDelta:
                xMid = xLimMin + 0.5*(xDelta)
                xLimMin = xMid - 0.5*yDelta*raToDecLimitRatio
                xLimMax = xMid + 0.5*yDelta*raToDecLimitRatio
            else:
                yMid = yLimMin + 0.5*(yDelta)
                yLimMin = yMid - 0.5*xDelta/raToDecLimitRatio
                yLimMax = yMid + 0.5*xDelta/raToDecLimitRatio
        else:
            if xDelta0 > yDelta0:
                yMid = yLimMin + 0.5*(yDelta)
                yLimMin = yMid - 0.5*xDelta/raToDecLimitRatio
                yLimMax = yMid + 0.5*xDelta/raToDecLimitRatio
            else:
                xMid = xLimMin + 0.5*(xDelta)
                xLimMin = xMid - 0.5*yDelta*raToDecLimitRatio
                xLimMax = xMid + 0.5*yDelta*raToDecLimitRatio
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
        range while perserving the central values.

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


def getDetRaDecCorners(ccdKey, ccdId, datasetType, visit, visitSummary=None, butler=None, doLogWarn=True):
    """Compute the RA/Dec corners lists for a given detector in a visit.
    """
    raCorners, decCorners = None, None
    if visitSummary is not None and datasetType != "raw":
        row = visitSummary.find(ccdId)
        if row is None:
            if doLogWarn:
                logger.warning("No row found for %d in visitSummary of visit %d. "
                            "Skipping and continuing...", ccdId, visit)
        else:
            raCorners = list(row["raCorners"])
            decCorners = list(row["decCorners"])
            if any(~np.isfinite(raCorners)) or any(~np.isfinite(decCorners)):
                logger.warning("Got nans for RA/Dec corners for %d in visitSummary of visit %d. "
                               "Skipping and continuing...", ccdId, visit)
                raCorners, decCorners = None, None
    else:
        try:
            if datasetType == "raw" or datasetType == "preliminary_visit_image":
                if datasetType == "raw":
                    dataId = {"exposure": visit, ccdKey: ccdId}
                else:
                    dataId = {"visit": visit, ccdKey: ccdId}
                exp = butler.get(datasetType, dataId)
                wcs = exp.getWcs()
                bbox = exp.getBBox()
            elif datasetType == "xpreliminary_visit_image":
                wcs = None
            else:
                dataId = {"visit": visit, ccdKey: ccdId}
                wcs = butler.get(datasetType + ".wcs", dataId)
                bbox = butler.get(datasetType + ".bbox", dataId)
            if wcs is not None:
                raCorners, decCorners = bboxToRaDec(bbox, wcs)
            else:
                logger.warning("WCS is None on %s detector %d of visit %d. "
                               "Skipping and continuing...", datasetType, ccdId, visit)
                raCorners, decCorners = None, None
        # except LookupError as e:
        except Exception as e:
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
        record = list(butler.registry.queryDimensionRecords("band", visit=visit))[0]
        band = record.name
        record = list(butler.registry.queryDimensionRecords("physical_filter", visit=visit))[0]
        physicalFilter = record.name
    return band, physicalFilter


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=str,
                        help="URI or path to an existing data repository root or configuration file")
    parser.add_argument("--collections", type=str, nargs="+",
                        help="Blank-space separated list of collection names for butler instantiation",
                        metavar=("COLLECTION1", "COLLECTION2"), required=True)
    parser.add_argument("--skymapName", default=None, help="Name of the skymap for the collection")
    parser.add_argument("--tracts", type=int, nargs="+", default=None,
                        help=("Blank-space separated list of tracts to constrain search for "
                              "visit overlap"), metavar=("TRACT1", "TRACT2"))
    parser.add_argument("--patches", type=int, nargs="+", default=None,
                        help=("Blank-space separated list of patches to constrain search for "
                              "visit overlap"), metavar=("PATCH1", "PATCH2"))
    parser.add_argument("--visits", type=int, nargs="+", default=None,
                        help="Blank-space separated list of visits to include",
                        metavar=("VISIT1", "VISIT2"))
    parser.add_argument("--physicalFilters", type=str, nargs="+", default=None,
                        help=("Blank-space separated list of physical filter names to constrain search for "
                              "visits"), metavar=("PHYSICAL_FILTER1", "PHYSICAL_FILTER2"))
    parser.add_argument("--bands", type=str, nargs="+", default=None,
                        help=("Blank-space separated list of canonical band names to constrin search for "
                              "visits"), metavar=("BAND1", "BAND2"))
    parser.add_argument("-c", "--ccds", nargs="+", type=int, default=None,
                        help="Blank-space separated list of CCDs to show", metavar=("CCD1", "CCD2"))
    parser.add_argument("-p", "--showPatch", action="store_true", default=False,
                        help="Show the patch boundaries")
    parser.add_argument("--saveFile", type=str, default="showVisitSkyMap.png",
                        help="Filename to write the plot to")
    parser.add_argument("--ccdKey", default="detector", help="Data ID name of the CCD key")
    parser.add_argument("--showCcds", action="store_true", default=False,
                        help="Show ccd ID numbers on output image.  Will try to avoid overlapping labels "
                        "by selectively skipping some.")
    parser.add_argument("--showCcdsAll", action="store_true", default=False,
                        help="Show all ccd ID numbers on output image")
    parser.add_argument("--visitVetoFile", type=str, default=None,
                        help="Full path to single-column file containing a list of visits to veto")
    parser.add_argument("--minOverlapFraction", type=float, default=None,
                        help="Minimum fraction of detectors that overlap any tract for visit to be included")
    parser.add_argument("--trimToTracts", action="store_true", default=False,
                        help="Set plot limits based on extent of visits (as opposed to tracts) plotted?")
    parser.add_argument("--doUnscaledLimitRatio", action="store_true", default=False,
                        help="Let axis limits get set by sky coordinate range without scaling to focal "
                        "plane based projection (ignored if --forceScaledLimitRatio is passed).")
    parser.add_argument("--forceScaledLimitRatio", action="store_true", default=False,
                        help="Force the axis limit scaling to focal plane based projection (takes "
                        "precedence over --doUnscaledLimitRatio.")
    parser.add_argument("--plotFailsOnly", action="store_true", default=False,
                        help="Only plot the raw WCS outlines for failed (i.e. no calexp) detectors.")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    main(args.repo, args.collections, skymapName=args.skymapName, tracts=args.tracts,
         patches=args.patches, visits=args.visits, physicalFilters=args.physicalFilters,
         bands=args.bands, ccds=args.ccds, ccdKey=args.ccdKey,
         showPatch=args.showPatch, saveFile=args.saveFile,
         showCcds=args.showCcds, showCcdsAll=args.showCcdsAll,
         visitVetoFile=args.visitVetoFile, minOverlapFraction=args.minOverlapFraction,
         trimToTracts=args.trimToTracts, doUnscaledLimitRatio=args.doUnscaledLimitRatio,
         forceScaledLimitRatio=args.forceScaledLimitRatio, plotFailsOnly=args.plotFailsOnly)
