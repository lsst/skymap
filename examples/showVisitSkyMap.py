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
import matplotlib.pyplot as plt

import lsst.afw.cameraGeom as cameraGeom
import lsst.daf.butler as dafButler
import lsst.geom as geom


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


def main(repo, collections, skymapName=None, tracts=None, visits=None, bands=None,
         ccds=None, ccdKey="ccd", showPatch=False, saveFile=None, showCcds=False):
    logger.info("Making butler for collections = {} in repo {}".format(collections, repo))
    butler = dafButler.Butler(repo, collections=collections)
    instrument = butler.registry.findDataset("camera").dataId["instrument"]
    # Make a guess at the skymapName if not provided
    if skymapName is None:
        if instrument == "HSC":
            skymapName = "hsc_rings_v1"
        elif instrument == "LSSTCam-imSim":
            skymapName = "DC2"
        elif instrument == "LATISS":
            skymapName = "latiss_v1"
        elif instrument == "DECam":
            skymapName = "decam_rings_v1"
        else:
            print("Unknown skymapName for instrument: {}.  Must specify --skymapName on command line".
                  format(instrument))
    print("instrument = {} skymapName = {}".format(instrument, skymapName))
    camera = butler.get("camera", instrument=instrument)
    skymap = butler.get("skyMap", instrument=instrument, skymap=skymapName)

    # draw the CCDs
    ras, decs = [], []
    bboxesPlotted = []
    cmap = get_cmap(len(visits))
    for i_v, visit in enumerate(visits):
        print("Working on visit %d [%d of %d]" % (visit, i_v + 1, len(visits)))
        inLegend = False
        color = cmap(i_v)
        for ccd in camera:
            bbox = ccd.getBBox()
            ccdId = int(ccd.getId())

            if (ccds is None or ccdId in ccds) and ccd.getType() == cameraGeom.DetectorType.SCIENCE:
                dataId = {"visit": visit, ccdKey: ccdId}
                try:
                    wcs = butler.get("calexp.wcs", dataId)
                except LookupError as e:
                    print(e, " Skip and continuing...")
                    continue

                ra, dec = bboxToRaDec(bbox, wcs)
                ras += ra
                decs += dec
                if not inLegend:
                    plt.fill(ra, dec, fill=True, alpha=0.2, color=color, edgecolor=color,
                             label=str(visit))
                    inLegend = True
                else:
                    plt.fill(ra, dec, fill=True, alpha=0.2, color=color, edgecolor=color)

                # add CCD serial numbers
                if showCcds:
                    minPoint = geom.Point2D(min(ra), min(dec))
                    maxPoint = geom.Point2D(max(ra), max(dec))
                    # Use doubles in Box2D to check overlap
                    bboxDouble = geom.Box2D(minPoint, maxPoint)
                    overlaps = [not bboxDouble.overlaps(otherBbox) for otherBbox in bboxesPlotted]
                    if all(overlaps):
                        plt.text(percent(ra), percent(dec), str(ccdId), fontsize=6,
                                 horizontalalignment="center", verticalalignment="center", color=color)
                        plt.fill(ra, dec, fill=False, alpha=0.5, color=color, edgecolor=color)
                    bboxesPlotted.append(bboxDouble)

    buff = 0.1
    xlim = max(ras)+buff, min(ras)-buff
    ylim = min(decs)-buff, max(decs)+buff

    # draw the skymap
    if showPatch:
        alpha0 = 1.0
        for i_t, tract in enumerate(tracts):
            alpha = max(0.1, alpha0 - i_t*1.0/len(tracts))
            inLegend = False
            for tractInfo in skymap:
                if tractInfo.getId() == tract:
                    for patch in tractInfo:
                        ra, dec = bboxToRaDec(patch.getInnerBBox(), tractInfo.getWcs())
                        if not inLegend:
                            plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed",
                                        alpha=alpha, label=str(tract))
                            inLegend = True
                        else:
                            plt.fill(ra, dec, fill=False, edgecolor="k", lw=1, linestyle="dashed",
                                        alpha=alpha)
                        if xlim[1] < percent(ra) < xlim[0] and ylim[0] < percent(dec) < ylim[1]:
                            plt.text(percent(ra), percent(dec),
                                        ((str(patch.getIndex()))[1:-1].replace(" ", "")), fontsize=6,
                                        horizontalalignment="center", verticalalignment="center", alpha=alpha)

    # add labels and save
    ax = plt.gca()
    ax.set_xlabel("R.A. (deg)")
    ax.set_ylabel("Decl. (deg)")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fancybox=True, shadow=True, fontsize=6)
    collectionStr = collections[0]
    if len(collections) > 1:
        for collection in collections[1:]:
            collectionStr += "\n" + collection
    ax.set_title("{}".format(collectionStr), fontsize=9)
    fig = plt.gcf()
    if saveFile is not None:
        fig.savefig(saveFile)
    else:
        fig.show()


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
    parser.add_argument("--bands", type=int, nargs="+", default=None,
                        help="Blank-space separated list of bands to include",
                        metavar="BAND1 [BAND2 [BAND3...]")
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
    args = parser.parse_args()
    main(args.repo, args.collections, skymapName=args.skymapName, tracts=args.tracts, visits=args.visits,
         ccds=args.ccds, ccdKey=args.ccdKey, showPatch=args.showPatch, saveFile=args.saveFile,
         showCcds=args.showCcds)
