import xarray as xr

from climate_toolbox.utils.utils import *


def load_climate_data(data_type, file_path):
    """ load_climate_data(data_type, file_path)
    Read and prepare climate data

    :param data_type: str
        datatype to be read, supported types:
        bcsd, gmfd, best, era5

    :param file_path: str
        File path

    :return: ds: xr.Dataset
         xarray dataset loaded in memory
    """

    return _load_climate_data(
        _find_loader(data_type),
        file_path
    )


def _load_climate_data(loader, file_path):
    with xr.open_dataset(file_path) as ds:
        ds.load()

    return loader(ds)


def load_min_max_temperatures(data_type, file_path_tmin, file_path_tmax):
    """ load_min_max_temperatures(data_type, file_path_tmin, file_path_tmax)

    :param data_type: str
        datatype to be read, supported types:
        bcsd, gmfd, best, era5

    :param file_path_tmin: path for min temperature
    :param file_path_tmax: path for max temperature
    :return:
        ds_tasmax: xr.Dataset, ds_tasmin: xr.Dataset
    """

    return _load_min_max_temperatures(
        _find_loader(data_type),
        file_path_tmin,
        file_path_tmax
    )


def _load_min_max_temperatures(loader, file_path_tmin, file_path_tmax):
    with xr.open_dataset(file_path_tmin) as ds_tasmin:
        ds_tasmin.load()
    with xr.open_dataset(file_path_tmax) as ds_tasmax:
        ds_tasmax.load()

    return loader(ds_tasmin), loader(ds_tasmax)


def _find_loader(data_type):
    """ Helper function to find climate data loader """

    data_type = data_type.lower()

    if 'bcsd' in data_type:
        loader = load_bcsd
    elif 'gmfd' in data_type:
        loader = load_gmfd
    elif 'best' in data_type:
        loader = load_best
    elif 'era' in data_type:
        loader = load_era5
    else:
        raise TypeError("'" + data_type + "' not supported. Supported data "
                                          "types are: NASA BCSD, GMFD, BEST, ERA5.")
    return loader


def standardize_climate_data(ds):
    """ standardize_climate_data(ds)
    Standardize climate data units to:
        - lon and lat,
        - lon to -180 to 180 and

    :param ds: xr.Dataset
    :return: ds: xr.Dataset
    """

    ds = rename_coords_to_lon_and_lat(ds)
    ds = convert_lons_split(ds, lon_name='lon')

    return ds


def load_bcsd(ds):
    return ds


def load_gmfd(ds):
    if 'tmin' in ds.data_vars or 'tmax' in ds.data_vars:
        return standardize_climate_data(ds)
    if 'lat' not in ds.coords or 'lon' not in ds.coords:
        ds = rename_coords_to_lon_and_lat(ds)
    return convert_lons_split(ds, lon_name="lon")


def load_best():
    pass


def load_era5():
    pass
