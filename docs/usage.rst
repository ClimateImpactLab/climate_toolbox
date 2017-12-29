=====
Usage
=====

`climate_toolbox` is a library that aids common climate data computing routines at the `Climate Impact Lab <http://impactlab.org>`_. Many of the main utilities are wrappers for `xarray <http://xarray.pydata.org>`_ functions. 

The main public functions in ``climate_toolbox`` are ``load_bcsd``, ``map_grid_to_region_segments`` and ``weighted_aggregate_grid_to_region``. 

``load_bcsd``
-------------

This function is specific to the gridded `Bias Corrected Spatially Downscaled <https://nex.nasa.gov/nex/projects/1356/>`_ (BCSD) dataset created by NASA. This method does some basic data cleaning by filling invalid values and renaming the longitude coordinates. In many, if not all, of the BCSD files, there are invalid `NaN` values and `load_bcsd` interpolates these values with a bounding box interpolation method. 

This function takes two required arguments: a file path and a varible name. The file path can be either a string formatted file path or a `xarray Dataset <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html?highlight=dataset>`_. In either case, the variable name must be in string format. 

.. code-block:: python
    
    import climate_toolbox.climate_toolbox as ctb

    path = '/data/climate/ACCESS1-0/rcp45/2018/tasmax.nc'
    ds = ctb.load_bcsd(path, 'tasmax')


You can also do the following. 

.. code-block:: python

    import climate_toolbox.climate_toolbox as ctb
    import xarray as xr

    path = '/data/climate/ACCESS1-0/rcp45/2018/tasmax.nc'
    ds = xr.open_dataset(path).load()
    ds = ctb.load_bcsd(ds, 'tasmax')

In both cases, the return value is an `xarray Dataset <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html?highlight=dataset>`_ whose `NaN` values have been removed


Mapping grids to regions
------------------------

Many times we have gridded climate data that we want to represent at an arbitrary regional level. `map_grid_to_region_segments` performs this routine with the help of `geopandas <https://geopandas.org>`_. The function takes a gridded shapefile and a regional shapefile as its two inputs. The method will return a `GeoDataFrame <http://geopandas.org/data_structures.html#geodataframe>`_ where each row represents the pixel centroid of the grid that intersects with the a given segment of a the geographic region. For example, if you have a map of the United States with 50 rows, and one row for each state, the operation would return all the grid-cell pixel centroid segments that are in the state. 

In addition to producing segment level representations, we also can compute the weighted amount of area or activity that each segment contributes within its region. So for example if in California we have have some raster file that shows the amount of cropland by pixel we can compute how much each pixel contributes to each state and then compute the cropweight for that pixel with the state. This will be used when we compute weighted climate variables. This operation is being handled by `rasterstas <http://pythonhosted.org/rasterstats/>`_ and its `zonalstats` function.

.. code-block:: python

    import climate_toolbox.climate_toolbox as ctb

    grid_shpfl = '/data/geospatial/climate_grid_025.shp'
    usa_shpfl = '/data/geospatial/usa_state_level.shp'
    raster_file = '/data/geospatial/crop_raster.tif'
    weighting = 'crop'

    usa_geodf = ctb.map_grid_to_region_segments(grid_shpfl, usa_shpfl, raster_file, weighting)

    usa_geodf.head()



Weighted Grids to Regions
-------------------------

 Weighted grid to regions performs the weighting and projection of gridded data onto regional data. This is another task common at the `Climate Impact Lab <http://impactlab.org>`_ where we want to transform gridded climate data to region specific climate data. Not only do we want to map to the region-level but we want to weight our climate data by some socioeconomic variable like population, area, or crop-planted. 

 This function wraps two main operations that are implemented in xarray. We first need to generate a common index between our grid and our regional/pixel_centroid mapping. To do this we use xarrays `advanced indexing <http://xarray.pydata.org/en/stable/indexing.html#more-advanced-indexing>`_. The second is doing some simple averaging to generate weighted-climate variables for each region.   


.. code-block:: python

    import climate_toolbox.climate_toolbox as ctb

    weights = usa_geodf


 





