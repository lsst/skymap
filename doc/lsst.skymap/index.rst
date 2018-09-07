.. py:currentmodule:: lsst.skymap

.. _lsst.skymap:

###########
lsst.skymap
###########

The skymap package provides tools for creating pixelizations of the sky to define tracts and patches for imaging data.

A sky map describes a pixelization of image data that covers most or all of the sky.
The imaging data is arranged as a sequence of overlapping rectangular "tracts".
Each tract is, in essence, a single large exposure.
However, tracts are typically too large to fit into memory, so tracts are subdivided into rectangular, possibly overlapping "patches".
The patch size is chosen to easily fit into memory.

.. toctree linking to topics related to using the module's APIs.

.. .. toctree::
..    :maxdepth: 1

.. _lsst.skymap-contributing:

Contributing
============

``lsst.skymap`` is developed at https://github.com/lsst/skymap.
You can find Jira issues for this module under the `skymap <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20skymap>`_ component.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 1

.. _lsst.skymap-pyapi:

Python API reference
====================

.. automodapi:: lsst.skymap
   :no-main-docstr:
   :no-inheritance-diagram:
.. automodapi:: lsst.skymap.detail
.. automodapi:: lsst.skymap.cachingSkyMap
.. automodapi:: lsst.skymap.healpixSkyMap
.. automodapi:: lsst.skymap.ringsSkyMap
