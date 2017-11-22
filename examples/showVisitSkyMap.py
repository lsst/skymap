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

from __future__ import print_function
from builtins import zip
from builtins import str

import argparse
import matplotlib.pyplot as pyplot

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.daf.persistence as dafPersist
from lsst.pipe.base.argumentParser import IdValueAction, DataIdContainer


def bboxToRaDec(bbox, wcs):
    """Get the corners of a BBox and convert them to lists of RA and Dec."""
    corners = []
    for corner in bbox.getCorners():
        p = afwGeom.Point2D(corner.getX(), corner.getY())
        coord = wcs.pixelToSky(p).toIcrs()
        corners.append([coord.getRa().asDegrees(), coord.getDec().asDegrees()])
    ra, dec = zip(*corners)
    return ra, dec


def percent(values, p=0.5):
    """Return a value a faction of the way between the min and max values in a list."""
    m = min(values)
    interval = max(values) - m
    return m + p*interval


def get_cmap(n, name='hsv'):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.
    """
    return pyplot.cm.get_cmap(name, n)


def main(rootDir, tracts, visits, ccds=None, ccdKey='ccd', showPatch=False, saveFile=None, showCcds=False):
    butler = dafPersist.Butler(rootDir)
    camera = butler.get("camera")

    # draw the CCDs
    ras, decs = [], []
    bboxesPlotted = []
    print('Number of visits = ', len(visits))
    cmap = get_cmap(len(visits))
    for i_v, visit in enumerate(visits):
        print("%r visit=%r" % (i_v + 1, visit))
        inLegend = False
        color = cmap(i_v)
        for ccd in camera:
            bbox = ccd.getBBox()
            ccdId = int(ccd.getSerial())

            if (ccds is None or ccdId in ccds) and ccd.getType() == cameraGeom.SCIENCE:
                dataId = {'visit': visit, ccdKey: ccdId}
                try:
                    md = butler.get("calexp_md", dataId)
                    wcs = afwImage.makeWcs(md)

                    ra, dec = bboxToRaDec(bbox, wcs)
                    ras += ra
                    decs += dec
                    if not inLegend:
                        pyplot.fill(ra, dec, fill=True, alpha=0.2, color=color, edgecolor=color,
                                    label=str(visit))
                        inLegend = True
                    else:
                        pyplot.fill(ra, dec, fill=True, alpha=0.2, color=color, edgecolor=color)

                    # add CCD serial numbers
                    if showCcds:
                        minPoint = afwGeom.Point2D(min(ra), min(dec))
                        maxPoint = afwGeom.Point2D(max(ra), max(dec))
                        # Use doubles in Box2D to check overlap
                        bboxDouble = afwGeom.Box2D(minPoint, maxPoint)
                        overlaps = [not bboxDouble.overlaps(otherBbox) for otherBbox in bboxesPlotted]
                        if all(overlaps):
                            pyplot.text(percent(ra), percent(dec), str(ccdId), fontsize=6,
                                        horizontalalignment='center', verticalalignment='center', color=color)
                            pyplot.fill(ra, dec, fill=False, alpha=0.5, color=color, edgecolor=color)
                        bboxesPlotted.append(bboxDouble)
                except:
                    pass

    buff = 0.1
    xlim = max(ras)+buff, min(ras)-buff
    ylim = min(decs)-buff, max(decs)+buff

    skymap = butler.get("deepCoadd_skyMap")
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
                            pyplot.fill(ra, dec, fill=False, edgecolor='k', lw=1, linestyle='dashed',
                                        alpha=alpha, label=str(tract))
                            inLegend = True
                        else:
                            pyplot.fill(ra, dec, fill=False, edgecolor='k', lw=1, linestyle='dashed',
                                        alpha=alpha)
                        if xlim[1] < percent(ra) < xlim[0] and ylim[0] < percent(dec) < ylim[1]:
                            pyplot.text(percent(ra), percent(dec),
                                        ((str(patch.getIndex()))[1:-1].replace(" ", "")), fontsize=6,
                                        horizontalalignment='center', verticalalignment='center', alpha=alpha)

    # add labels and save
    ax = pyplot.gca()
    ax.set_xlabel("R.A. (deg)")
    ax.set_ylabel("Decl. (deg)")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True, shadow=True, fontsize=6)
    fig = pyplot.gcf()
    if saveFile is not None:
        fig.savefig(saveFile)
    else:
        fig.show()


def splitId(argName):
    class SplitIdValueAction(IdValueAction):

        def __call__(self, parser, namespace, values, option_string):
            # Hack to use IdValueAction
            keyValues = [argName + "=" + str(values[0])]
            setattr(namespace, "config", "hack")
            setattr(namespace, argName, DataIdContainer())
            # Parse the data into namespace.argName.idList
            IdValueAction.__call__(self, parser, namespace, keyValues, "--"+argName)
            # Save the list into namespace.argName
            setattr(namespace, argName,
                    list({int(dataId[argName]) for dataId in getattr(namespace, argName).idList}))
    return SplitIdValueAction

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("tracts", nargs=1, action=splitId("tracts"),
                        help="Tract(s) to show", metavar="TRACT1[^TRACT2[^TRACT3...]")
    parser.add_argument("visits", nargs=1, action=splitId("visits"),
                        help="Visits to show", metavar="VISIT1[^VISIT2[^VISIT3...]")
    parser.add_argument("-c", "--ccds", nargs=1, action=splitId("ccds"), default=None,
                        help="CCDs to show", metavar="CCD1[^CCD2[^CCD3...]")
    parser.add_argument("-p", "--showPatch", action='store_true', default=False,
                        help="Show the patch boundaries")
    parser.add_argument("--saveFile", type=str, default=None,
                        help="Filename to write the plot to")
    parser.add_argument("--ccdKey", default="ccd", help="Data ID name of the CCD key")
    parser.add_argument("--showCcds", action='store_true', default=False,
                        help="Show ccd serial numbers on output image")
    args = parser.parse_args()

    main(args.root, args.tracts, visits=args.visits, ccds=args.ccds,
         ccdKey=args.ccdKey, showPatch=args.showPatch, saveFile=args.saveFile, showCcds=args.showCcds)
