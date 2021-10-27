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
import abc

import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

__all__ = ["patchBuilderRegistry", "BasePatchBuilder",
           "OldPatchBuilder", "CellPatchBuilder"]


class BasePatchBuilderConfig(pexConfig.Config):
    """Configuration that is to be shared amongst all patch builders."""
    pass


class BasePatchBuilder(pipeBase.Task, metaclass=abc.ABCMeta):
    """Base class for patch builders.

    Parameters
    ----------
    config : `lsst.pexConfig.Config`
        Input for configuring the algorithm
    """
    ConfigClass = BasePatchBuilderConfig
    _DefaultName = "patchBuilder"

    def __init__(self, config, **kwargs):
        pipeBase.Task.__init__(self, config=config, **kwargs)

    @abc.abstractmethod
    def setupPatches(self):
        """Set up the patches.

        Returns
        -------
        """
        pass


class OldPatchBuilderConfig(BasePatchBuilderConfig):
    patchInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of patches (x,y pixels)",
        dtype=int,
        length=2,
        default=(4000, 4000),
    )
    patchBorder = pexConfig.Field(
        doc="border between patch inner and outer bbox (pixels)",
        dtype=int,
        default=100,
    )


class OldPatchBuilder(BasePatchBuilder):
    ConfigClass = OldPatchBuilderConfig
    _DefaultName = "oldPatchBuilder"

    def setupPatches(self):
        pass


class CellPatchBuilderConfig(BasePatchBuilderConfig):
    cellInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of cells (x,y pixels)",
        dtype=int,
        length=2,
        default=(150, 150),
    )
    cellBorder = pexConfig.Field(
        doc="Border between cell inner and outer bbox (pixels)",
        dtype=int,
        default=50,
    )


class CellPatchBuilder(BasePatchBuilder):
    ConfigClass = CellPatchBuilderConfig
    _DefaultName = "cellPatchBuilder"

    def setupPatches(self):
        pass


patchBuilderRegistry = pexConfig.makeRegistry(
    doc="A registry of Patch Builders (subclasses of PatchBuilder)",
)

patchBuilderRegistry.register("old", OldPatchBuilder)
patchBuilderRegistry.register("cells", CellPatchBuilder)
