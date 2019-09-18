import xarray as xr

from climate_toolbox.utils.utils import (
    rename_coords_to_lon_and_lat,
    convert_lons_split,
)

import climate_toolbox.io.sources as sources

def standardize_climate_data(
        ds, lon_name=lon_name, lat_name=lat_name, source=None):
    """
    Standardize dimensions and units for common surface climate data variables

    Standardizes:
        - dimension names: lon and lat
        - lon to -180 to 180 and
        - for known sources, drops length 0 or 1 z/alt/pres dimension
        - for known sources, standardizes variable names for tavg, tmin, tmax,
          & precip

    Parameters
    ----------
    ds:  xr.Dataset
    lon_name

    Returns
    -------
    xr.Dataset
    """

    if source is not None:
        ds = sources.standardize_source_data(ds, source=source)
    else:
        ds = rename_coords_to_lon_and_lat(ds)
        ds = convert_lons_split(ds, lon_name=lon_names)

    return ds


def load_and_standardize_climate_dataset(
        ds, lon_name=None, lat_name=None, source=None):
    """
    Read and prepare climate data

    After reading data, this method also fills NA values using linear
    interpolation, and standardizes longitude to -180:180

    Parameters
    ----------
    fp: str
        File path or dataset
    lon_name : str, optional
        Name of the longitude dimension, which will be standardized to -180 to
        180 and renamed "lon" (defualt selects from ['lon', 'lng', 'long', or
        'longitude'])
    lat_name : str, optional
        Name of the latitude dimension, which will be renamed "lat" (default)
        selects from ['lat', 'latitude']
    source: str, optional
        Source name, used to standardize variable names for data from known
        sources. See climate_toolbox.io.sources.KNOWN_SOURCES for a list of
        available sources.

    Returns
    -------
    xr.Dataset
         xarray dataset loaded into memory
    """

    if isinstance(fp, (xr.Dataset, xr.DataArray)):
        ds = fp

    else:
        with xr.open_dataset(fp) as ds:
            ds.load()

    return standardize_climate_data(
        ds, lon_name=lon_name, lat_name=lat_name, source=source)

def load_and_standardize_multiple_climate_datasets(
        file_specs,
        lon_name=None,
        lat_name=None,
        source=None):
    """
    Read and prepare multiple climate datasets

    After reading data, this method also standardizes dimension names and
    variable names for datasets from known sources.

    Parameters
    ----------
    file_specs: dict
        Dictionary of {variable name: filepath or Dataset} pairs that will be
        read in (if a path) and standardized.
    lon_name : str, optional
        Name of the longitude dimension, which will be standardized to -180 to
        180 and renamed "lon" (defualt selects from ['lon', 'lng', 'long', or
        'longitude'])
    lat_name : str, optional
        Name of the latitude dimension, which will be renamed "lat" (default)
        selects from ['lat', 'latitude']
    source: str, optional
        Source name, used to standardize variable names for data from known
        sources. See climate_toolbox.io.sources.KNOWN_SOURCES for a list of
        available sources.

    Returns
    -------
    dict
         Dictionary with variable name, Dataset pairs, with the Datasets loaded
         into memory and reformatted.

    Examples
    --------

        >>> load_and_standardize_multiple_climate_datasets(
        ...     {
        ...         'tasmin': '/shares/data/tasmin.nc',
        ...         'tasmax': '/shares/data/tasmax.nc'},
        ...     source='NASA-NEX/GDDP')  # doctest: +SKIP
        ...
        {'tasmin': <xarray.Dataset ...>, 'tasmax': <xarray.Dataset ...>}
    """

    res = {}

    for varname, path in file_specs.items():
        if isinstance(path, (xr.Dataset, xr.DataArray)):
            ds = path

        else:
            with xr.open_dataset(path) as ds:
                ds.load()

        res[varname] = standardize_climate_data(
            ds, lon_name=lon_name, lat_name=lat_name, source=source)

    return res

def load_and_standardize_climate_dataset_from_pattern(
        pattern, lon_name=None, lat_name=None, source=None, **kwargs):

    return load_and_standardize_climate_dataset(
        pattern.format(**kwargs),
        lon_name=lon_name,
        lat_name=lat_name,
        source=source)

def load_and_standardize_multiple_climate_datasets_from_patterns(
        patterns, lon_name=None, lat_name=None, source=None, **kwargs):
    """
    Read and prepare multiple climate datasets

    After reading data, this method also standardizes dimension names and
    variable names for datasets from known sources.

    Parameters
    ----------
    patterns: dict
        Dictionary of {variable name: string file pattern} pairs that will be
        populated with keyword arguments and read in, then standardized.
    lon_name : str, optional
        Name of the longitude dimension, which will be standardized to -180 to
        180 and renamed "lon" (defualt selects from ['lon', 'lng', 'long', or
        'longitude'])
    lat_name : str, optional
        Name of the latitude dimension, which will be renamed "lat" (default)
        selects from ['lat', 'latitude']
    source: str, optional
        Source name, used to standardize variable names for data from known
        sources. See climate_toolbox.io.sources.KNOWN_SOURCES for a list of
        available sources.

    **kwargs used to populate file patterns

    Returns
    -------
    dict
         Dictionary with variable name, Dataset pairs, with the Datasets loaded
         into memory and reformatted.

    Examples
    --------

        >>> load_and_standardize_multiple_climate_datasets_from_patterns(
        ...     {
        ...         'tasmin': '/shares/data/tasmin/{rcp}/{model}/{year}.nc',
        ...         'tasmax': '/shares/data/tasmax/{rcp}/{model}/{year}.nc'},
        ...     source='NASA-NEX/GDDP',
        ...     rcp='rcp85',
        ...     model='CCSM4',
        ...     year=2005)  # doctest: +SKIP
        ...
        {'tasmin': <xarray.Dataset ...>, 'tasmax': <xarray.Dataset ...>}
    """

    res = {
        k: load_and_standardize_climate_dataset(
            pattern.format(**kwargs),
            lon_name=lon_name,
            lat_name=lat_name,
            source=source)
        for k, pattern in patterns.items()}

    return res
