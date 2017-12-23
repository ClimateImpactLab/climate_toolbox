=======
History
=======

0.1.4 (current version)
-----------------------

* Support vectorized indexing for xarray >= 0.10 in :py:func:`climate_toolbox.climate_toolbox._reindex_spatial_data_to_regions` (:issue:`10`)
* Support iteratively increasing bounding box in :py:func:`~climate_toolbox.climate_toolbox._fill_holes_xr` (:issue:`11`).
* Support multiple interpolation methods (linear and cubic) in :py:func:`~climate_toolbox.climate_toolbox._fill_holes_xr` (:issue:`12`).
* Fix bug causing tests to pass no matter what

0.1.3 (2017-08-04)
------------------

* Support passing a dataset (not just a filepath) into ``load_baseline`` and ``load_bcsd`` (:issue:`4`)

0.1.2 (2017-07-25)
------------------

* merge in bug fixes

0.1.1 (2017-07-25)
-----------------------

* Various bug fixes (see :issue:`2`)


0.1.0 (2017-07-24)
------------------

* First release on PyPI.
