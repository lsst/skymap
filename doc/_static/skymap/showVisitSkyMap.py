#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#


import argparse
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord
from matplotlib.legend import Legend

import lsst.afw.cameraGeom as cameraGeom
import lsst.daf.butler as dafButler
import lsst.geom as geom
import lsst.log as lsstLog
import lsst.sphgeom as sphgeom

logger = lsstLog.Log.getLogger("showVisitSkyMap.py")


def bboxToRaDec(bbox, wcs):
    """Get the corners of a BBox and convert them to lists of RA and Dec."""
    sphPoints = wcs.pixelToSky(geom.Box2D(bbox).getCorners())
    ra = [float(sph.getRa().asDegrees()) for sph in sphPoints]
    dec = [float(sph.getDec().asDegrees()) for sph in sphPoints]
    return ra, dec


def percent(values, p=0.5):
    """Return a value a faction of the way between the min and max values in a
    list.
    """
    m = min(values)
    interval = max(values) - m
    return m + p*interval


def get_cmap(n, name="hsv"):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.
    """
    return matplotlib.colormaps[name].resampled(n)


def main(repo, collections, skymapName=None, tracts=None, visits=None, physicalFilters=None, bands=None,
         ccds=None, ccdKey="detector", showPatch=False, saveFile=None, showCcds=False, visitVetoFile=None,
         minOverlapFraction=None, trimToTracts=False):
    if minOverlapFraction is not None and tracts is None:
        raise RuntimeError("Must specify --tracts if --minOverlapFraction is set")
    logger.info("Making butler for collections = {} in repo {}".format(collections, repo))
    butler = dafButler.Butler(repo, collections=collections)
    instrument = butler.registry.findDataset("camera").dataId["instrument"]
    detectorSkipList = []
    # Make a guess at the skymapName if not provided
    if skymapName is None:
        if instrument == "HSC":
            skymapName = "hsc_rings_v1"
            detectorSkipList = [9]  # detector 9 has long been dead for HSC
        elif instrument == "LSSTCam-imSim":
            skymapName = "DC2"
        elif instrument == "LATISS":
            skymapName = "latiss_v1"
        elif instrument == "DECam":
            skymapName = "decam_rings_v1"
        else:
            logger.error("Unknown skymapName for instrument: %s.  Must specify --skymapName on command line.",
                         instrument)
    logger.info("instrument = {} skymapName = {}".format(instrument, skymapName))
    camera = butler.get("camera", instrument=instrument)
    skymap = butler.get("skyMap", instrument=instrument, skymap=skymapName)

    whereStr = ""
    if tracts is not None:
        tractStr = makeWhereInStr("tract", tracts, int)
        whereStr += tractStr

    if physicalFilters is not None:
        physicalFilterStr = makeWhereInStr("physical_filter", physicalFilters, str)
        whereStr += " AND " + physicalFilterStr if len(whereStr) else " " + physicalFilterStr

    if bands is not None:
        bandStr = makeWhereInStr("band", bands, str)
        whereStr += " AND " + bandStr if len(whereStr) else " " + bandStr

    if len(whereStr) > 1:
        whereStr = "instrument=\'" + instrument + "\' AND skymap=\'" + skymapName + "\' AND " + whereStr
    else:
        whereStr = "instrument=\'" + instrument + "\' AND skymap=\'" + skymapName + "\'"
    logger.info("Applying the following where clause in dataId search: {} ".format(whereStr))

    visitVetoList = []
    if visitVetoFile is not None:
        with open(visitVetoFile) as f:
            content = f.readlines()
            visitVetoList = [int(visit.strip()) for visit in content]

    if visits is None:
        dataRefs = list(butler.registry.queryDatasets("calexp", where=whereStr).expanded())
        visits = []
        for dataRef in dataRefs:
            visit = dataRef.dataId.visit.id
            if visit not in visits and visit not in visitVetoList:
                visits.append(visit)
        visits.sort()
        logger.info("List of visits (N={}) satisfying where and veto clauses: {}".format(len(visits),
                                                                                         visits))
    else:
        if len(visitVetoList) > 1:
            visitListTemp = visits.copy()
            for visit in visitListTemp:
                if visit in visitVetoList:
                    visits.remove(visit)
        logger.info("List of visits (N={}) excluding veto list: {}".format(len(visits), visits))

    ccdIdList = []
    for ccd in camera:
        ccdId = ccd.getId()
        if ((ccds is None or ccdId in ccds) and ccd.getType() == cameraGeom.DetectorType.SCIENCE
                and ccdId not in detectorSkipList):
            ccdIdList.append(ccdId)
    ccdIdList.sort()
    nDetTot = len(ccdIdList)

    visitIncludeList = []
    # Determine the fraction of detectors that overlap any tract under
    # consideration. If this fraction does not exceed minOverlapFraction,
    # skip this visit.
    if minOverlapFraction is not None:
        for i_v, visit in enumerate(visits):
            ccdOverlapList = []
            try:
                visitSummary = butler.get("visitSummary", visit=visit)
            except LookupError as e:
                logger.warn("%s  Will try to get wcs from calexp.", e)
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
                                ccdKey, ccdId, visit, visitSummary=visitSummary, butler=butler,
                                doLogWarn=False)
                            if raCorners is not None and decCorners is not None:
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

    # draw the CCDs
    ras, decs = [], []
    bboxesPlotted = []
    cmap = get_cmap(len(visitIncludeList))
    alphaEdge = 0.7
    maxVisitForLegend = 20
    finalVisitList = []
    includedBands = []
    includedPhysicalFilters = []
    for i_v, visit in enumerate(visitIncludeList):
        print("Working on visit %d [%d of %d]" % (visit, i_v + 1, len(visitIncludeList)), end="\r")
        inLegend = False
        color = cmap(i_v)
        fillKwargs = {"fill": False, "alpha": alphaEdge, "facecolor": None, "edgecolor": color, "lw": 0.6}
        try:
            visitSummary = butler.get("visitSummary", visit=visit)
        except LookupError as e:
            logger.warn("%s  Will try to get wcs from calexp.", e)
            visitSummary = None

        band, physicalFilter = getBand(visitSummary=visitSummary, butler=butler, visit=visit)
        if band not in includedBands:
            includedBands.append(band)
        if physicalFilter not in includedPhysicalFilters:
            includedPhysicalFilters.append(physicalFilter)

        for ccdId in ccdIdList:
            raCorners, decCorners = getDetRaDecCorners(
                ccdKey, ccdId, visit, visitSummary=visitSummary, butler=butler)
            if raCorners is not None and decCorners is not None:
                ras += raCorners
                decs += decCorners
                if not inLegend and len(visitIncludeList) <= maxVisitForLegend:
                    plt.fill(raCorners, decCorners, label=str(visit), **fillKwargs)
                    inLegend = True
                else:
                    plt.fill(raCorners, decCorners, **fillKwargs)
                plt.fill(raCorners, decCorners, fill=True, alpha=alphaEdge/4, color=color,
                         edgecolor=color)
                if visit not in finalVisitList:
                    finalVisitList.append(visit)
                # add CCD serial numbers
                if showCcds:
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
                    if all(overlaps):
                        plt.text(percent(raCorners), percent(decCorners), str(ccdId), fontsize=6,
                                 ha="center", va="center", color="darkblue")
                        bboxesPlotted.append(bboxDouble)

    logger.info("Final list of visits (N={}) satisfying where and minOverlapFraction clauses: {}"
                .format(len(finalVisitList), finalVisitList))
    tractList = list(set(skymap.findTractIdArray(ras, decs, degrees=True)))
    tractList.sort()
    logger.info("List of tracts overlapping data:  {}".format(tractList))
    tractLimitsDict = getTractLimitsDict(skymap, tractList)

    minVisitRa, maxVisitRa = min(ras), max(ras)
    minVisitDec, maxVisitDec = min(decs), max(decs)
    raVisitDiff = maxVisitRa - minVisitRa
    decVisitDiff = maxVisitDec - minVisitDec
    midVisitRa = minVisitRa + 0.5*raVisitDiff
    midVisitDec = minVisitDec + 0.5*decVisitDiff
    midRa = np.atleast_1d((midVisitRa*units.deg).to(units.radian).value).astype(np.float64)
    midDec = np.atleast_1d((midVisitDec*units.deg).to(units.radian).value).astype(np.float64)
    midSkyCoord = SkyCoord(midVisitRa*units.deg, midVisitDec*units.deg)

    # Find a detector that contains the mid point in RA/Dec (or the closest
    # one) to set the plot aspect ratio.
    minDistToMidCood = 1e12
    minSepVisit = None
    minSepCcdId = None
    raToDecLimitRatio = None
    for i_v, visit in enumerate(visits):
        try:
            visitSummary = butler.get("visitSummary", visit=visit)
        except LookupError as e:
            logger.warn("%s  Will try to get wcs from calexp.", e)
            visitSummary = None
        for ccdId in ccdIdList:
            raCorners, decCorners = getDetRaDecCorners(
                ccdKey, ccdId, visit, visitSummary=visitSummary, butler=butler, doLogWarn=False)
            if raCorners is not None and decCorners is not None:
                detSphCorners = []
                for ra, dec in zip(raCorners, decCorners):
                    pt = geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                          geom.Angle(dec, geom.degrees))
                    detSphCorners.append(pt)
                    ptSkyCoord = SkyCoord(ra*units.deg, dec*units.deg)
                    separation = (midSkyCoord.separation(ptSkyCoord)).degree
                    if separation < minDistToMidCood:
                        minSepVisit = visit
                        minSepCcdId = ccdId
                        minDistToMidCood = separation
                detConvexHull = sphgeom.ConvexPolygon(
                    [coord.getVector() for coord in detSphCorners])
                if detConvexHull.contains(midRa, midDec):
                    logger.info("visit/det overlapping plot coord mid point in RA/Dec: %d %d", visit, ccdId)
                    raToDecLimitRatio = (max(raCorners) - min(raCorners))/(max(decCorners) - min(decCorners))
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

    if raToDecLimitRatio is None:
        try:
            visitSummary = butler.get("visitSummary", visit=minSepVisit)
        except LookupError as e:
            logger.warn("%s  Will try to get wcs from calexp.", e)
            visitSummary = None
        raCorners, decCorners = getDetRaDecCorners(
            ccdKey, minSepCcdId, minSepVisit, visitSummary=visitSummary, butler=butler, doLogWarn=False)
        for ra, dec in zip(raCorners, decCorners):
            pt = geom.SpherePoint(geom.Angle(ra, geom.degrees),
                                  geom.Angle(dec, geom.degrees))
            detSphCorners.append(pt)
        detConvexHull = sphgeom.ConvexPolygon([coord.getVector() for coord in detSphCorners])
        logger.info("visit/det closest to plot coord mid point in RA/Dec (none actually overlat it): %d %d",
                    minSepVisit, minSepCcdId)
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
        visitLimitsDict = {"allVisits": {"ras": [minVisitRa, maxVisitRa], "decs": [minVisitDec, maxVisitDec]}}
        xlim, ylim = derivePlotLimits(visitLimitsDict, raToDecLimitRatio=raToDecLimitRatio, buffFrac=0.04)

    # draw the skymap
    alpha0 = 1.0
    tractHandleList = []
    tractStrList = []
    if tracts is not None:
        tractOutlineList = list(set(tracts + tractList))
    else:
        tractOutlineList = tractList
    tractOutlineList.sort()
    logger.info("List of tract outlines being plotted: {}".format(tractOutlineList))
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
                plt.text(tCenterRa, tCenterDec, tract, fontsize=7, alpha=alpha, ha="center", va="center")
        ra, dec = bboxToRaDec(tractInfo.bbox, tractInfo.getWcs())
        plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed", alpha=alpha)
        tractArtist = mpatches.Patch(fill=False, edgecolor="k", linestyle="dashed", alpha=alpha)
        tractHandleList.append(tractArtist)
        tractStrList.append(str(tract))
        if showPatch:
            patchColor = "k"
            for patch in tractInfo:
                ra, dec = bboxToRaDec(patch.getInnerBBox(), tractInfo.getWcs())
                plt.fill(ra, dec, fill=False, edgecolor=patchColor, lw=0.5, linestyle=(0, (5, 6)),
                         alpha=alpha)
                if (xlim[1] + fracDeltaX < percent(ra) < xlim[0] - fracDeltaX
                        and ylim[0] + fracDeltaY < percent(dec) < ylim[1] - fracDeltaY):
                    plt.text(percent(ra), percent(dec), str(patch.sequential_index), fontsize=5,
                             color=patchColor, ha="center", va="center", alpha=alpha)

    # add labels and save
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_box_aspect(1)
    if abs(xlim[1] > 99.99):
        ax.tick_params("x", labelrotation=45, pad=0, labelsize=8)
    else:
        ax.tick_params("x", labelrotation=0, pad=0, labelsize=8)
    ax.tick_params("y", labelsize=8)
    ax.set_xlabel("RA (deg)", fontsize=9)
    ax.set_ylabel("Dec (deg)", fontsize=9)

    visitScaleOffset = None
    if len(visitIncludeList) > maxVisitForLegend:
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
        tractLegend = Legend(ax, tractHandleList, tractStrList, loc="upper right",
                             fancybox=True, shadow=True, fontsize=5, title_fontsize=6, title="tracts")
        ax.add_artist(tractLegend)
    else:
        if len(visitIncludeList) > 0:
            ax.legend(loc="center left", bbox_to_anchor=(1.15, 0.5), fancybox=True, shadow=True,
                      fontsize=6, title_fontsize=6, title="visits")
        # Create the second legend and add the artist manually.
        tractLegend = Legend(ax, tractHandleList, tractStrList, loc="center left", bbox_to_anchor=(1.0, 0.5),
                             fancybox=True, shadow=True, fontsize=6, title_fontsize=6, title="tracts")
        ax.add_artist(tractLegend)

    titleStr = repo + "\n" + collections[0]
    if len(collections) > 1:
        for collection in collections[1:]:
            titleStr += "\n" + collection
    titleStr += "\nnVisit: {}".format(str(len(finalVisitList)))
    if minOverlapFraction is not None:
        titleStr += " (minOvlpFrac = {:.2f})".format(minOverlapFraction)
    if len(includedBands) > 0:
        titleStr += "  bands: {}".format(str(includedBands).translate({ord(i): None for i in "[]\'"}))
    if len(includedPhysicalFilters) > 0:
        if len(includedPhysicalFilters[0]) > 9:
            titleStr += "\n"
        titleStr += "  physical filters: {}".format(str(includedPhysicalFilters).translate(
            {ord(i): None for i in "[]\'"}))
    ax.set_title("{}".format(titleStr), fontsize=8)

    fig = plt.gcf()
    if saveFile is not None:
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
    xLimMin, yLimMin = 1e12, 1e12
    xLimMax, yLimMax = -1e12, -1e12
    for limitId, limits in limitsDict.items():
        xLimMin = min(xLimMin, min(limits["ras"]))
        xLimMax = max(xLimMax, max(limits["ras"]))
        yLimMin = min(yLimMin, min(limits["decs"]))
        yLimMax = max(yLimMax, max(limits["decs"]))
    xDelta0 = xLimMax - xLimMin
    yDelta0 = yLimMax - yLimMin

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


def getDetRaDecCorners(ccdKey, ccdId, visit, visitSummary=None, butler=None, doLogWarn=True):
    """Compute the RA/Dec corners lists for a given detector in a visit.
    """
    raCorners, decCorners = None, None
    if visitSummary is not None:
        row = visitSummary.find(ccdId)
        if row is None:
            if doLogWarn:
                logger.warn("No row found for %d in visitSummary of visit %d. "
                            "Skipping and continuing...", ccdId, visit)
        else:
            raCorners = list(row["raCorners"])
            decCorners = list(row["decCorners"])
    else:
        try:
            dataId = {"visit": visit, ccdKey: ccdId}
            wcs = butler.get("calexp.wcs", dataId)
            bbox = butler.get("calexp.bbox", dataId)
            raCorners, decCorners = bboxToRaDec(bbox, wcs)
        except LookupError as e:
            logger.warn("%s Skipping and continuing...", e)

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
                        help=("Blank-space separated list of tract outlines to constrain search for "
                              "visit overlap"), metavar=("TRACT1", "TRACT2"))
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
                        help="Show ccd ID numbers on output image")
    parser.add_argument("--visitVetoFile", type=str, default=None,
                        help="Full path to single-column file containing a list of visits to veto")
    parser.add_argument("--minOverlapFraction", type=float, default=None,
                        help="Minimum fraction of detectors that overlap any tract for visit to be included")
    parser.add_argument("--trimToTracts", action="store_true", default=False,
                        help="Set plot limits based on extent of visits (as opposed to tracts) plotted?")
    args = parser.parse_args()
    main(args.repo, args.collections, skymapName=args.skymapName, tracts=args.tracts, visits=args.visits,
         physicalFilters=args.physicalFilters, bands=args.bands, ccds=args.ccds, ccdKey=args.ccdKey,
         showPatch=args.showPatch, saveFile=args.saveFile, showCcds=args.showCcds,
         visitVetoFile=args.visitVetoFile, minOverlapFraction=args.minOverlapFraction,
         trimToTracts=args.trimToTracts)
