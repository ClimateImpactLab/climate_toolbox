import xarray as xr

from climate_toolbox.utils.utils import *


def standardize_climate_data(ds):
    """
    Read climate data and standardize units to:
        - lon and lat,
        - lon to -180 to 180 and

    Parameters
    ----------
    ds:  xr.Dataset

    Returns
    -------
    xr.Dataset
    """

    ds = rename_coords_to_lon_and_lat(ds)
    ds = convert_lons_split(ds, lon_name='lon')

    return ds


def load_bcsd(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    """
    Read and prepare climate data

    After reading data, this method also fills NA values using linear
    interpolation, and standardizes longitude to -180:180

    Parameters
    ----------
    fp: str
        File path or dataset

    varname: str
        Variable name to be read

    lon_name : str, optional
        Name of the longitude dimension (defualt selects from ['lon' or
        'longitude'])

    Returns
    -------
    xr.Dataset
         xarray dataset loaded into memory
    """

    if lon_name is not None:
        lon_names = [lon_name]

    if hasattr(fp, 'sel_points'):
        ds = fp

    else:
        with xr.open_dataset(fp) as ds:
            ds.load()

    return standardize_climate_data(ds)


def load_gmfd(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    pass


def load_best(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    pass

