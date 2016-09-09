#!/usr/bin/env python
import lsst.utils.tests
import unittest

import lsst.afw.geom as afwGeom

from lsst.skymap.discreteSkyMap import DiscreteSkyMap
from SkyMapTestCase import SkyMapTestCase


# These are the PS1 Medium-Deep fields
coords = [(10.6750, 41.2667),  # M31
          (36.2074, -04.5833),  # XMM-LSS
          (53.1, -28.1333),  # CDFS
          (130.5917, +44.3167),  # IfA-Lynx
          (150.0, +02.2),  # COSMOS
          (161.9167, +58.0833),  # Lockman
          (185.0, +47.1167),  # NGC4258
          (213.7051, +53.0834),  # DEEP2-Field1/EGS
          (242.7875, +54.95),  # EliasN1
          (334.1875, +00.2833),  # SA22
          (352.3125, -00.4333),  # DEEP2-Field3
          (270.0, +66.56),  # North Ecliptic Pole
          ]
config = DiscreteSkyMap.ConfigClass()
config.raList = [c[0] for c in coords]
config.decList = [c[1] for c in coords]
config.radiusList = [2] * len(coords)


class DiscreteTestCase(SkyMapTestCase):

    def setUp(self):
        s_cls = SkyMapTestCase
        s_cls._NumTracts = len(coords)  # Number of tracts to expect
        s_cls._NeighborAngularSeparation = None  # Expected tract separation
        s_cls._SkyMapClass = DiscreteSkyMap  # Class of SkyMap to test
        s_cls._SkyMapConfig = config  # Configuration to use
        s_cls._SkyMapName = "discrete"  # Name of SkyMap class to test
        s_cls._numNeighbors = None  # Number of neighbours

    def testTractSeparation(self):
        self.skipTest("A particular tract separation is not important for DiscreteSkyMap")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
