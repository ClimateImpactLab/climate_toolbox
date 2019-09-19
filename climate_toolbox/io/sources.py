
from __future__ import absolute_import

from climate_toolbox.utils.utils import (
    rename_coords_to_lon_and_lat,
    convert_lons_split,
)

__all__ = ['KNOWN_SOURCES', 'standardize_source_data']

def standardize_berkeley_earth_data(ds):
    raise NotImplementedError(
        "needs work to combine climatology and temperature variables")

def standardize_kelvin_temperature_data(ds):
    '''
    Convert kelvin to degrees C and assert reasonable surface temperature bounds
    '''

    if 'tas' in ds.data_vars:
        all_values_null_or_valid = (
            ((ds['tas'] > -55 + 273.15) & (ds['tas'] < 50 + 273.15))
            | (ds['tas'].isnull())).all()

        assert all_values_null_or_valid, "values outside range [-55, 55] C"

        # convert to Celsius
        ds['tas'] = ds['tas'] - 273.15

    if 'tasmin' in ds.data_vars:
        all_values_null_or_valid = (
            ((ds['tasmin'] > -75 + 273.15) & (ds['tasmin'] < 40 + 273.15))
            | (ds['tasmin'].isnull())).all()

        assert all_values_null_or_valid, "values outside range [-55, 55] C"

        # convert to Celsius
        ds['tasmin'] = ds['tasmin'] - 273.15

    if 'tasmax' in ds.data_vars:
        all_values_null_or_valid = (
            ((ds['tasmax'] > -45 + 273.15) & (ds['tasmax'] < 55 + 273.15))
            | (ds['tasmax'].isnull())).all()

        assert all_values_null_or_valid, "values outside range [-55, 55] C"

        # convert to Celsius
        ds['tasmax'] = ds['tasmax'] - 273.15

def standardize_global_meterological_forcing_dataset_v1(ds):

    # standardize dimension names
    ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})

    # remove length-1 vertical dimension
    assert len(ds.z) <= 1, 'expected length 1 z dimention'
    ds = ds.isel(z=0, drop=True)

    # standardize variable names
    ds = ds.rename({
        'tmax': 'tasmax',
        'tmin': 'tasmin',
        'tavg': 'tas',
        'prcp': 'precip'})

    standardize_kelvin_temperature_data(ds)

    return ds

def standardize_global_meterological_forcing_dataset_v3(ds):

    assert len(ds.dims) == 3, 'expected 3 dimensions: (lon, lat, time)'

    for dim in ['lon', 'lat', 'time']:
        assert dim in ds.dims, "dimension {} not found".format(dim)

    # standardize variable names
    ds = ds.rename({
        'tmax': 'tasmax',
        'tmin': 'tasmin',
        'tavg': 'tas',
        'prcp': 'precip'})

    standardize_kelvin_temperature_data(ds)

    return ds

def standardize_nasa_nex_gddp_dataset(ds):
    raise NotImplementedError(
        "Mike couldn't find these on sacagawea and got lazy")

def standardize_climate_impact_lab_smme_pattern_scaled_nasa_nex_dataset(ds):
    raise NotImplementedError(
        "Mike couldn't find these on sacagawea and got lazy")

_source_standardizer_lookup = {
    'BerkeleyEarth': standardize_berkeley_earth_data,
    'GMFDv1': standardize_global_meterological_forcing_dataset_v1,
    'GMFDv3': standardize_global_meterological_forcing_dataset_v3,
    'NASA-NEX/GDDP': standardize_nasa_nex_gddp_dataset,
    'CIL-SMME-NASA/NEX-GDDP': (
        standardize_climate_impact_lab_smme_pattern_scaled_nasa_nex_dataset),
}

KNOWN_SOURCES = list(_source_standardizer_lookup.keys())

def standardize_source_data(ds, source):
    '''
    Calls source-specific functions which reformat climate data to a standard format

    Parameters
    ----------
    ds : Dataset
        :py:class:`xarray.Dataset` to standardize
    source: str
        Known source with an associated standardization function. See
        climate_toolbox.io.sources.KNOWN_SOURCES for a list of available
        sources.

    Returns
    -------
    ds : Dataset
        Standardized dataset
    '''

    ds = _source_standardizer_lookup[source](ds)
    return ds
