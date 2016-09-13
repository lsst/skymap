from builtins import object
#
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
__all__ = ["PatchInfo"]


class PatchInfo(object):
    """Information about a patch within a tract of a sky map

    See TractInfo for more information.
    """

    def __init__(self, index, innerBBox, outerBBox):
        """Construct a PatchInfo

        @param[in] index: x,y index of patch (a pair of ints)
        @param[in] innerBBox: inner bounding box (an afwGeom.Box2I)
        @param[in] outerBBox: inner bounding box (an afwGeom.Box2I)
        """
        self._index = index
        self._innerBBox = innerBBox
        self._outerBBox = outerBBox
        if not outerBBox.contains(innerBBox):
            raise RuntimeError("outerBBox=%s does not contain innerBBox=%s" % (outerBBox, innerBBox))

    def getIndex(self):
        """Get patch index
        """
        return self._index

    def getInnerBBox(self):
        """Get inner bounding box
        """
        return self._innerBBox

    def getOuterBBox(self):
        """Get outer bounding box
        """
        return self._outerBBox

    def __eq__(self, rhs):
        """Support ==
        """
        return (self.getIndex() == rhs.getIndex()) \
            and (self.getInnerBBox() == rhs.getInnerBBox()) \
            and (self.getOuterBBox() == rhs.getOuterBBox())

    def __ne__(self, rhs):
        """Support !=
        """
        return not self.__eq__(rhs)

    def __str__(self):
        """Return a brief string representation
        """
        return "PatchInfo(index=%s)" % (self.getIndex(),)

    def __repr__(self):
        """Return a detailed string representation
        """
        return "PatchInfo(index=%s, innerBBox=%s, outerBBox=%s)" % \
            (self.getIndex(), self.getInnerBBox(), self.getOuterBBox())
