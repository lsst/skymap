# This file is part of skymap.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ("SkyMapDimensionPacker",)

from collections.abc import Mapping

from lsst.pex.config import Config, Field, DictField, ConfigurableField
from lsst.daf.butler import DimensionPacker, DimensionGraph, DataCoordinate
from deprecated.sphinx import deprecated


class SkyMapDimensionPackerConfig(Config):
    bands = DictField(
        "Mapping from band name to integer to use in the packed ID. "
        "The default (None) is to use a hard-coded list of common bands; "
        "pipelines that can enumerate the set of bands they are likely to see "
        "should override this.",
        keytype=str,
        itemtype=int,
        default=None,
        optional=True,
    )
    n_bands = Field(
        "Number of bands to reserve space for. "
        "If zero, bands are not included in the packed integer at all. "
        "If `None`, the size of 'bands' is used.",
        dtype=int,
        optional=True,
        default=0,
    )
    n_tracts = Field(
        "Number of tracts, or, more precisely, one greater than the maximum tract ID."
        "Default (None) obtains this value from the skymap dimension record.",
        dtype=int,
        optional=True,
        default=None,
    )
    n_patches = Field(
        "Number of patches per tract, or, more precisely, one greater than the maximum patch ID."
        "Default (None) obtains this value from the skymap dimension record.",
        dtype=int,
        optional=True,
        default=None,
    )


class SkyMapDimensionPacker(DimensionPacker):
    """A `DimensionPacker` for tract, patch and optionally band,
    given a SkyMap.

    Parameters
    ----------
    fixed : `lsst.daf.butler.DataCoordinate`
        Data ID that identifies just the ``skymap`` dimension.  Must have
        dimension records attached unless ``n_tracts`` and ``n_patches`` are
        not `None`.
    dimensions : `lsst.daf.butler.DimensionGraph`, optional
        The dimensions of data IDs packed by this instance.  Must include
        ``{skymap, tract, patch}``, and may include ``band``.  If not provided,
        this will be set to include ``band`` if ``n_bands != 0``.
    bands : `~collections.abc.Mapping` [ `str`, `int` ] or `None`, optional
        Mapping from band name to integer to use in the packed ID.  `None` uses
        a fixed set of bands defined in this class.  When calling code can
        enumerate the bands it is likely to see, providing an explicit mapping
        is preferable.
    n_bands : `int` or `None`, optional
        The number of bands to leave room for in the packed ID.  If `None`,
        this will be set to ``len(bands)``.  If ``0``, the band will not be
        included in the dimensions at all.  If ``1``, the band will be included
        in the dimensions but will not occupy any extra bits in the packed ID.
        This may be larger or smaller than ``len(bands)``, to reserve extra
        space for the future or align to byte boundaries, or support a subset
        of a larger mapping, respectively.
    n_tracts : `int` or `None`, optional
        The number of tracts to leave room for in the packed ID.  If `None`,
        this will be set via the ``skymap`` dimension record in ``fixed``.
    n_patches : `int` or `None`, optional
        The number of patches (per tract) to leave room for in the packed ID.
        If `None`, this will be set via the ``skymap`` dimension record in
        ``fixed``.

    Notes
    -----
    The standard pattern for constructing instances of this class is to use
    `make_config_field`::

        class SomeConfig(lsst.pex.config.Config):
            packer = ObservationDimensionPacker.make_config_field()

        class SomeTask(lsst.pipe.base.Task):
            ConfigClass = SomeConfig

            def run(self, ..., data_id: DataCoordinate):
                packer = self.config.packer.apply(data_id)
                packed_id = packer.pack(data_id)
                ...

    """

    SUPPORTED_FILTERS = (
        [None]
        + list("ugrizyUBGVRIZYJHK")  # split string into single chars
        + [f"N{d}" for d in (387, 515, 656, 816, 921, 1010)]  # HSC narrow-bands
        + [f"N{d}" for d in (419, 540, 708, 964)]  # DECam narrow-bands
    )
    """Sequence of supported bands used to construct a mapping from band name
    to integer when the 'bands' config option is `None` or no config is
    provided.

    This variable should no longer be modified to add new filters; pass
    ``bands`` at construction or use `from_config` instead.
    """

    ConfigClass = SkyMapDimensionPackerConfig

    @classmethod
    @deprecated(
        reason="This classmethod cannot reflect all __init__ args and will be removed after v27.",
        version="v26.0",
        category=FutureWarning,
    )
    def getIntFromFilter(cls, name):
        """Return an integer that represents the band with the given
        name.
        """
        try:
            return cls.SUPPORTED_FILTERS.index(name)
        except ValueError:
            raise NotImplementedError(f"band '{name}' not supported by this ID packer.")

    @classmethod
    @deprecated(
        reason="This classmethod cannot reflect all __init__ args and will be removed after v27.",
        version="v26.0",
        category=FutureWarning,
    )
    def getFilterNameFromInt(cls, num):
        """Return an band name from its integer representation."""
        return cls.SUPPORTED_FILTERS[num]

    @classmethod
    @deprecated(
        reason="This classmethod cannot reflect all __init__ args and will be removed after v27.",
        version="v26.0",
        category=FutureWarning,
    )
    def getMaxIntForFilters(cls):
        return len(cls.SUPPORTED_FILTERS)

    def __init__(
        self,
        fixed: DataCoordinate,
        dimensions: DimensionGraph | None = None,
        bands: Mapping[str, int] | None = None,
        n_bands: int | None = None,
        n_tracts: int | None = None,
        n_patches: int | None = None,
    ):
        if bands is None:
            bands = {b: i for i, b in enumerate(self.SUPPORTED_FILTERS)}
        if dimensions is None:
            if n_bands is None:
                n_bands = len(bands)
            dimension_names = ["tract", "patch"]
            if n_bands != 0:
                dimension_names.append("band")
            dimensions = fixed.universe.extract(dimension_names)
        else:
            if "band" not in dimensions.names:
                n_bands = 0
                if dimensions.names != {"tract", "patch", "skymap"}:
                    raise ValueError(
                        f"Invalid dimensions for skymap dimension packer with n_bands=0: {dimensions}."
                    )
            else:
                if dimensions.names != {"tract", "patch", "skymap", "band"}:
                    raise ValueError(
                        f"Invalid dimensions for skymap dimension packer with n_bands>0: {dimensions}."
                    )
                if n_bands is None:
                    n_bands = len(bands)
        if n_tracts is None:
            n_tracts = fixed.records["skymap"].tract_max
        if n_patches is None:
            n_patches = (
                fixed.records["skymap"].patch_nx_max
                * fixed.records["skymap"].patch_ny_max
            )
        super().__init__(fixed, dimensions)
        self._bands = bands
        self._n_bands = n_bands
        self._n_tracts = n_tracts
        self._n_patches = n_patches
        self._bands_list = None

    @classmethod
    def make_config_field(
        cls,
        doc: str = "How to pack tract, patch, and possibly band into an integer."
    ) -> ConfigurableField:
        """Make a config field to control how skymap data IDs are packed.

        Parameters
        ----------
        doc : `str`, optional
            Documentation for the config field.

        Returns
        -------
        field : `lsst.pex.config.ConfigurableField`
            A config field whose instance values are [wrapper proxies to]
            `SkyMapDimensionPackerConfig` instances.
        """
        return ConfigurableField(doc, target=cls.from_config, ConfigClass=cls.ConfigClass)

    @classmethod
    def from_config(
        cls, data_id: DataCoordinate, config: SkyMapDimensionPackerConfig
    ) -> SkyMapDimensionPacker:
        """Construct a dimension packer from a config object and a data ID.

        Parameters
        ----------
        data_id : `lsst.daf.butler.DataCoordinate`
            Data ID that identifies at least the ``skymap`` dimension.  Must
            have dimension records attached unless ``config.n_tracts`` and
            ``config.n_patches`` are both not `None`.
        config : `SkyMapDimensionPackerConfig`
            Configuration object.

        Returns
        -------
        packer : `SkyMapDimensionPackerConfig`
            New dimension packer.

        Notes
        -----
        This interface is provided for consistency with the `lsst.pex.config`
        "Configurable" concept, and is invoked when ``apply(data_id)`` is
        called on a config instance attribute that corresponds to a field
        created by `make_config_field`.  The constructor signature cannot play
        this role easily for backwards compatibility reasons.
        """
        return cls(
            data_id.subset(data_id.universe.extract(["skymap"])),
            n_bands=config.n_bands,
            bands=config.bands,
            n_tracts=config.n_tracts,
            n_patches=config.n_patches,
        )

    @property
    def maxBits(self) -> int:
        # Docstring inherited from DataIdPacker.maxBits
        packedMax = self._n_tracts * self._n_patches
        if self._n_bands:
            packedMax *= self._n_bands
        return (packedMax - 1).bit_length()

    def _pack(self, dataId: DataCoordinate) -> int:
        # Docstring inherited from DataIdPacker.pack
        if dataId["patch"] >= self._n_patches:
            raise ValueError(f"Patch ID {dataId['patch']} is out of bounds; expected <{self._n_patches}.")
        if dataId["tract"] >= self._n_tracts:
            raise ValueError(f"Tract ID {dataId['tract']} is out of bounds; expected <{self._n_tracts}.")
        packed = dataId["patch"] + self._n_patches * dataId["tract"]
        if self._n_bands:
            if (band_index := self._bands.get(dataId["band"])) is None:
                raise ValueError(
                    f"Band {dataId['band']!r} is not supported by SkyMapDimensionPacker "
                    f"configuration; expected one of {list(self._bands)}."
                )
            if band_index >= self._n_bands:
                raise ValueError(
                    f"Band index {band_index} for {dataId['band']!r} is out of bounds; "
                    f"expected <{self._n_bands}."
                )
            packed += self._bands[dataId["band"]] * self._n_patches * self._n_tracts
        return packed

    def unpack(self, packedId: int) -> DataCoordinate:
        # Docstring inherited from DataIdPacker.unpack
        d = {"skymap": self.fixed["skymap"]}
        if self._n_bands:
            index, packedId = divmod(packedId, (self._n_tracts * self._n_patches))
            if self._bands_list is None:
                self._bands_list = list(self._bands)
            d["band"] = self._bands_list[index]
        d["tract"], d["patch"] = divmod(packedId, self._n_patches)
        return DataCoordinate.standardize(d, graph=self.dimensions)
