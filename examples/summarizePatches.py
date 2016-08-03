#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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

from __future__  import print_function

import argparse

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.daf.persistence as dafPersist
from lsst.pipe.base.argumentParser import IdValueAction, DataIdContainer

def wcsBBox(bbox, wcs):
    """Convert a Box2I with int pixel corners to Box2D with double WCS corners."""
    corners = []
    # get corners as WCS coordinates
    for corner in bbox.getCorners():
        p = afwGeom.Point2D(corner.getX(), corner.getY())
        coord = wcs.pixelToSky(p).toIcrs()
        corners.append([coord.getRa().asDegrees(), coord.getDec().asDegrees()])
    ra, dec = zip(*corners)

    # make Box2D with WCS coordinates
    minPoint = afwGeom.Point2D(min(ra), min(dec))
    maxPoint = afwGeom.Point2D(max(ra), max(dec))
    wcsBBox = afwGeom.Box2D(minPoint, maxPoint)
    return wcsBBox

def main(rootDir, visits=None, ccdKey='ccd', saveOverview=False, saveSummaries=None):
    """Find which ccds overlap which patch and create summaries of visits organized by patch."""
    butler = dafPersist.Butler(rootDir)
    camera = butler.get('camera')
    skymap = butler.get('deepCoadd_skyMap')

    # open files for overview and summaries if needed
    if saveOverview:
        overview = open("patchesOverview", 'w')
    if saveSummaries is not None:
        expandedFileName = "expanded_" + saveSummaries
        condensedFileName = "condensed_" + saveSummaries
        expandedSummary = open(expandedFileName, 'w')
        condensedSummary = open(condensedFileName, 'w')

    # Default to all visits if none are specified
    if visits is None:
        visits = butler.queryMetadata('calexp_md', 'visit')

    # Make BBoxes of each ccd in each visit
    skymapBBoxes = {}
    for visit in visits:
        for ccd in camera:
            ccdId = int(ccd.getSerial())
            bbox = ccd.getBBox()
            dataId = {'visit': visit, ccdKey: ccdId}
            try:
                md = butler.get('calexp_md', dataId)
                wcs = afwImage.makeWcs(md)
                ccdBBox = wcsBBox(bbox, wcs)
                skymapBBoxes[(visit, ccdId)] = ccdBBox
            except:
                pass

    # Run through patches and find overlaps with ccds
    for tract in skymap:
        for patch in tract:
            visitsToCcds = {}
            filtersToVisits = {}
            numOverlaps = 0
            patchId = patch.getIndex()
            patchBBox = wcsBBox(patch.getInnerBBox(), tract.getWcs())

            # Find which ccds the patch overlaps (if any)
            for skymapId in skymapBBoxes:
                if patchBBox.overlaps(skymapBBoxes[skymapId]):
                    visitId, ccdId = skymapId
                    if visitId not in visitsToCcds:
                        visitsToCcds[visitId] = [ccdId]
                        visitFilter = butler.queryMetadata('calexp_md', 'filter', visit=visitId)[0]
                        if visitFilter not in filtersToVisits:
                            filtersToVisits[visitFilter] = []
                        filtersToVisits[visitFilter] += [visitId]
                        numOverlaps += 1
                    else:
                        visitsToCcds[visitId] += [ccdId]
                        numOverlaps += 1

            # If there was an overlap, print and save information
            if numOverlaps != 0:
                print("patch=%r:    Found %r overlapping ccds." % (patchId, numOverlaps))

                if saveOverview:
                    overview.write("patch={}\n".format(patchId))
                    overview.write("    {}\n".format(filtersToVisits))

                if saveSummaries is not None:
                    for i_f, patchFilter in enumerate(filtersToVisits):

                        # Write expanded file - one line for each filter in each patch
                        expandedSummary.write("filter={} ".format(patchFilter))
                        expandedSummary.write("tract={} ".format(tract.getId()))
                        expandedSummary.write("patch=")
                        expandedSummary.write("{},{} ".format(patchId[0],patchId[1]))
                        visits_string = "--selectId visit=" + "^".join(str(v) for v in filtersToVisits[patchFilter])
                        expandedSummary.write(visits_string)
                        expandedSummary.write("\n")

                        # Write condensed file - one line for each patch
                        if i_f != 0:
                            condensedSummary.write("^{}".format(patchFilter))
                        else:
                            condensedSummary.write("filter={}".format(patchFilter))
                    condensedSummary.write(" tract={} ".format(tract.getId()))
                    condensedSummary.write("patch=")
                    condensedSummary.write("{},{}\n".format(patchId[0],patchId[1]))

                # Save an overview of ccd/patch overlaps
                if saveOverview:
                    for visit in visitsToCcds:
                        overview.write("        visit={}\n".format(visit))
                        for ccd in visitsToCcds[visit]:
                            overview.write("            ccd={}\n".format(ccd))
                    overview.write("\n")

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
    parser.add_argument("--visits", nargs=1, action=splitId("visits"), default=None,
                        help="Visits to show", metavar="VISIT1[^VISIT2[^VISIT3...]")
    parser.add_argument("--ccdKey", default="ccd", help="Data ID name of the CCD key")
    parser.add_argument("--saveOverview", action='store_true', default=False,
                        help="Save detalied patch summary in current directory")
    parser.add_argument("--saveSummaries", type=str, default=None,
                        help=("Save two patch summary files: "
                              "(1) condensed: one line for each patch - e.g. filter=i^z tract=0 patch=8,1; "
                              "(2) expanded: one line for each filter in each patch - e.g. "
                              "filter=i tract=0 patch=8,1 --selectId visit=232735"))
    args = parser.parse_args()

    main(args.root, visits=args.visits, ccdKey=args.ccdKey, saveOverview=args.saveOverview,
         saveSummaries=args.saveSummaries)
