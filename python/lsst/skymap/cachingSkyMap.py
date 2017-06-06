from builtins import range
#
# LSST Data Management System
# Copyright 2008-2012 LSST Corporation.
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

from .baseSkyMap import BaseSkyMap

__all__ = ["CachingSkyMap"]

class CachingSkyMap(BaseSkyMap):
    """A SkyMap that generates its tracts on request and caches them

    A subclass should define
    * __init__ to calculate the required number of tracts (and pass it up)
    * generateTract to generate a tract

    Subclassers should also check that the arguments to the constructor are
    consistent with the below __reduce__ method.
    """

    def __init__(self, numTracts, config=None, version=0):
        super(CachingSkyMap, self).__init__(config)
        self._numTracts = numTracts
        self._tractCache = [None] * self._numTracts
        self._tractInfo = None # We shouldn't need this; we will generate tracts on demand
        self._version = version

    def __reduce__(self):
        """To support pickling

        Warning: This method assumes that the constructor should be defined:
            __init__(self, config, version=defaultVersion)
        The use of 'config' is effectively set by the registry mechanism.
        If additional optional arguments are added, this method should be
        overridden to correspond.
        """
        return (self.__class__, (self.config, self._version))

    def __iter__(self):
        """Iterator over tracts"""
        for i in range(self._numTracts):
            yield self[i]

    def __len__(self):
        """Length is number of tracts"""
        return self._numTracts

    def __getitem__(self, index):
        """Get the TractInfo for a particular index

        The tract is returned from a cache, if available, otherwise generated
        on the fly.
        """
        if index < 0 or index > self._numTracts:
            raise IndexError("Index out of range: %d vs %d" % (index, self._numTracts))
        if self._tractCache[index] is not None:
            return self._tractCache[index]
        tract = self.generateTract(index)
        self._tractCache[index] = tract
        return tract

    def generateTract(self, index):
        """Generate the TractInfo for the particular index"""
        raise NotImplementedError("Subclasses must define this method.")
