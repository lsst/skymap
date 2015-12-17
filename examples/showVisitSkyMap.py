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

import argparse
import re
import matplotlib.pyplot as pyplot

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.daf.persistence as dafPersist

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

def main(rootDir, tract, visits, ccds=None, showPatch=False):
    butler = dafPersist.Butler(rootDir)
    mapper = butler.mapper
    camera = mapper.camera

    # draw the CCDs
    ras, decs = [], []
    for i_v, visit in enumerate(visits):
        print("%r visit=%r" % (i_v, visit))
        for ccd in camera:
            bbox = ccd.getBBox()
            ccdId = int(ccd.getSerial())

            if (ccds is None or ccdId in ccds) and ccdId < 104:
                dataId = {'visit': visit, 'ccd': ccdId}
                try:
                    md = butler.get("calexp_md", visit=visit, ccd=ccdId)
                    wcs = afwImage.makeWcs(md)

                    ra, dec = bboxToRaDec(bbox, wcs)
                    ras += ra
                    decs += dec
                    color = ('r', 'b', 'c', 'g', 'm')[i_v%5]
                    pyplot.fill(ra, dec, fill=True, alpha=0.2, color=color, edgecolor=color)
                except:
                    pass

    buff = 0.1
    xlim = max(ras)+buff, min(ras)-buff
    ylim = min(decs)-buff, max(decs)+buff

    # draw the skymap
    if showPatch:
        skymap = butler.get('deepCoadd_skyMap', {'tract': 0})
        for tract in skymap:
            for patch in tract:
                ra, dec = bboxToRaDec(patch.getInnerBBox(), tract.getWcs())
                pyplot.fill(ra, dec, fill=False, edgecolor='k', lw=1, linestyle='dashed')
                if xlim[1] < percent(ra) < xlim[0] and ylim[0] < percent(dec) < ylim[1]:
                    pyplot.text(percent(ra), percent(dec, 0.9), str(patch.getIndex()),
                                fontsize=6, horizontalalignment='center', verticalalignment='top')

    # add labels and save
    ax = pyplot.gca()
    ax.set_xlabel("R.A. (deg)")
    ax.set_ylabel("Decl. (deg)")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig = pyplot.gcf()
    fig.savefig("patches.png")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("tract", type=int, help="Tract to show")
    parser.add_argument("visits", help="visit to show")
    parser.add_argument("-c", "--ccds", help="specify CCDs")
    parser.add_argument("-p", "--showPatch", action='store_true', default=False,
                        help="Show the patch boundaries")
    args = parser.parse_args()

    def idSplit(id):
        if id is None:
            return id
        ids = []
        for r in id.split("^"):
            m = re.match(r"^(\d+)\.\.(\d+):?(\d+)?$", r)
            if m:
                limits = [int(v) if v else 1 for v in m.groups()]
                limits[1] += 1
                ids += range(*limits)
            else:
                ids.append(int(r))
        return ids

    main(args.root, args.tract, visits=idSplit(args.visits), ccds=idSplit(args.ccds), showPatch=args.showPatch)
