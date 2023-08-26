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
         ccds=None, ccdKey="ccd", showPatch=False, saveFile=None, showCcds=False, visitVetoFile=None):
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
        logger.info("List of visits (N={}) satisfying where clause: {}".format(len(visits), visits))

    # draw the CCDs
    ras, decs = [], []
    bboxesPlotted = []
    cmap = get_cmap(len(visits))
    alphaEdge = 0.4
    maxVisitForLegend = 20
    visitFoundList = []
    for i_v, visit in enumerate(visits):
        logger.info("Working on visit %d [%d of %d]" % (visit, i_v + 1, len(visits)))
        inLegend = False
        color = cmap(i_v)
        try:
            visitSummary = butler.get("visitSummary", visit=visit)
            expTime = visitSummary[0].getVisitInfo().exposureTime
            if expTime < minExpTime:
                logger.warn("Skipping visit %d with expTime of %.1f (< minExpTime of %.1f)",
                            visit, expTime, minExpTime)
                continue
            visitFoundList.append(visit)
        except LoookupError as e:
            logger.warn("%s  Will try to get wcs from calexp.", e)
            visitSummary = None
        for ccd in camera:
            bbox = ccd.getBBox()
            ccdId = int(ccd.getId())

            if ((ccds is None or ccdId in ccds) and ccd.getType() == cameraGeom.DetectorType.SCIENCE
                    and ccdId not in detectorSkipList):
                dataId = {"visit": visit, ccdKey: ccdId}
                if visitSummary is None:
                    try:
                        wcs = butler.get("calexp.wcs", dataId)
                        ra, dec = bboxToRaDec(bbox, wcs)
                        if visit not in visitFoundList:
                            visitFoundList.append(visit)
                    except LookupError as e:
                        logger.warn("%s Skipping and continuing...", e)
                        continue
                else:
                    row = visitSummary.find(ccdId)
                    if row is None:
                        logger.warn("No row found for %d in visitSummary of visit %d. "
                                    "Skipping and continuing...", ccdId, visit)
                        continue
                    ra = list(row["raCorners"])
                    dec = list(row["decCorners"])
                ras += ra
                decs += dec
                if not inLegend and len(visitFoundList) <= maxVisitForLegend:
                    plt.fill(ra, dec, fill=True, alpha=alphaEdge, facecolor="none", edgecolor=color,
                             label=str(visit))
                    inLegend = True
                else:
                    plt.fill(ra, dec, fill=True, alpha=alphaEdge, facecolor="none", edgecolor=color)
                plt.fill(ra, dec, fill=True, alpha=alphaEdge/2, color=color, edgecolor=color)

                # add CCD serial numbers
                if showCcds:
                    minPoint = geom.Point2D(min(ra), min(dec))
                    maxPoint = geom.Point2D(max(ra), max(dec))
                    # Use doubles in Box2D to check overlap
                    bboxDouble = geom.Box2D(minPoint, maxPoint)
                    overlaps = [not bboxDouble.overlaps(otherBbox) for otherBbox in bboxesPlotted]
                    if all(overlaps):
                        plt.text(percent(ra), percent(dec), str(ccdId), fontsize=6,
                                 horizontalalignment="center", verticalalignment="center", color="darkblue")
                    bboxesPlotted.append(bboxDouble)

    logger.info("Final list of visits (N={}) satisfying where clause and expTime > {}s: {}"
                .format(len(visitFoundList), minExpTime, visits))
    tractList = list(set(skymap.findTractIdArray(ras, decs, degrees=True)))
    tractList.sort()
    logger.info("List of tracts overlapping data:  {}".format(tractList))
    tractLimitsDict = getTractLimitsDict(skymap, tractList)
    if doTrimToTract:
        xlim, ylim = trimToTract(tractLimitsDict)
    else:
        buff = 0.1
        xlim = max(ras) + buff, min(ras) - buff
        ylim = min(decs) - buff, max(decs) + buff

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
            plt.text(tCenterRa, tCenterDec, tract, fontsize=9, alpha=alpha)
        ra, dec = bboxToRaDec(tractInfo.bbox, tractInfo.getWcs())
        plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed", alpha=alpha)
        tractArtist = mpatches.Patch(fill=False, edgecolor="k", linestyle="dashed", alpha=alpha)
        tractHandleList.append(tractArtist)
        tractStrList.append(str(tract))
        if showPatch:
            for patch in tractInfo:
                ra, dec = bboxToRaDec(patch.getInnerBBox(), tractInfo.getWcs())
                plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed",
                         alpha=alpha)
                if xlim[1] < percent(ra) < xlim[0] and ylim[0] < percent(dec) < ylim[1]:
                    plt.text(percent(ra), percent(dec),
                             ((str(patch.getIndex()))[1:-1].replace(" ", "")), fontsize=6,
                             horizontalalignment="center", verticalalignment="center", alpha=alpha)

    # add labels and save
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.axis("equal")
    if abs(xlim[1] > 99.99):
        ax.tick_params("x", labelrotation=45, pad=0, labelsize=8)
    else:
        ax.tick_params("x", labelrotation=0, pad=0, labelsize=8)
    ax.tick_params("y", labelsize=8)
    ax.set_xlabel("RA (deg)", fontsize=9)
    ax.set_ylabel("Dec (deg)", fontsize=9)

    if len(visitFoundList) > maxVisitForLegend:
        nz = matplotlib.colors.Normalize()
        colorBarScale = visitFoundList
        nz.autoscale(colorBarScale)
        cax, _ = matplotlib.colorbar.make_axes(plt.gca(), pad=0.03)
        cax.tick_params(labelsize=7)
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=nz, alpha=alphaEdge)
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
        colorBarLabel = "visit number"
        cb.set_label(colorBarLabel, rotation=-90, labelpad=13, fontsize=9)
        tractLegend = Legend(ax, tractHandleList, tractStrList, loc="upper right",
                             fancybox=True, shadow=True, fontsize=5, title_fontsize=6, title="tracts")
        ax.add_artist(tractLegend)
    else:
        ax.legend(loc="center left", bbox_to_anchor=(1.15, 0.5), fancybox=True, shadow=True,
                  fontsize=6, title="visits")
        # Create the second legend and add the artist manually.
        tractLegend = Legend(ax, tractHandleList, tractStrList, loc="center left", bbox_to_anchor=(1.0, 0.5),
                             fancybox=True, shadow=True, fontsize=6, title="tracts")
        ax.add_artist(tractLegend)

    titleStr = repo + "\n" + collections[0]
    if len(collections) > 1:
        for collection in collections[1:]:
            titleStr += "\n" + collection
    if bands is not None:
        titleStr += "\nnVisit: {}  bands: {}".format(str(len(visitFoundList)), str(bands).strip("[]\'"))
    if physicalFilters is not None:
        titleStr += "  physical filterss: {}".format(str(physicalFilters).strip("[]\'"))
    ax.set_title("{}".format(titleStr), fontsize=8)

    fig = plt.gcf()
    if saveFile is not None:
        fig.savefig(saveFile, bbox_inches="tight", dpi=150)
    else:
        fig.show()


def makeWhereInStr(parameterName, parameterList, parameterType):
    typeStr = "\'" if parameterType is str else ""
    whereInStr = parameterName + " IN (" + typeStr + str(parameterList[0])
    if len(parameterList) > 1:
        for param in parameterList[1:]:
            whereInStr += typeStr + ", " + typeStr + str(param) + typeStr
    else:
        whereInStr += typeStr
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


def trimToTract(tractLimitsDict):
    xLimMin, yLimMin = 1e12, 1e12
    xLimMax, yLimMax = -1e12, -1e12
    for tract, tractLimits in tractLimitsDict.items():
        xLimMin = min(xLimMin, min(tractLimits["ras"]))
        xLimMax = max(xLimMax, max(tractLimits["ras"]))
        yLimMin = min(yLimMin, min(tractLimits["decs"]))
        yLimMax = max(yLimMax, max(tractLimits["decs"]))
    xDelta = xLimMax - xLimMin
    yDelta = yLimMax - yLimMin
    buffFrac = 0.04
    if xDelta > yDelta:
        xLimMin -= buffFrac*yDelta
        xLimMax += buffFrac*yDelta
    else:
        yLimMin -= buffFrac*yDelta
        yLimMax += buffFrac*yDelta
    xLimMin, xLimMax, yLimMin, yLimMax = setLimitsToEqualRatio(xLimMin, xLimMax, yLimMin, yLimMax)
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
    parser.add_argument("--doTrimToTract", action="store_true", default=True,
                        help="Set plot limits based on extent of tracts plotted?")
    args = parser.parse_args()
    main(args.repo, args.collections, skymapName=args.skymapName, tracts=args.tracts, visits=args.visits,
         physicalFilters=args.physicalFilters, bands=args.bands, ccds=args.ccds, ccdKey=args.ccdKey,
         showPatch=args.showPatch, saveFile=args.saveFile, showCcds=args.showCcds,
         visitVetoFile=args.visitVetoFile, doTrimToTract=args.doTrimToTract)
