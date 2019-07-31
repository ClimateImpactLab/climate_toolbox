#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `climate_toolbox` package."""

import pytest

from climate_toolbox.utils.utils import *
from climate_toolbox.aggregations.aggregations import \
    _reindex_spatial_data_to_regions, _aggregate_reindexed_data_to_regions
from climate_toolbox.transformations.transformations import snyder_edd, snyder_gdd
from climate_toolbox.io import *

import numpy as np
import pandas as pd
import xarray as xr


# python utils

# create pytest resource to be used across tests

@pytest.fixture(scope='session')
def lat():

    return np.arange(-89.875, 90, 2)


@pytest.fixture(scope='session')
def lon():

    return np.arange(0.125, 360.0, 2)


@pytest.fixture(scope='session')
def time():

    return pd.date_range(start=pd.datetime(2000, 1, 1), periods=10, freq='D')


@pytest.fixture
def clim_data(lat, lon, time):
    """
    Generate fake climate data to do tests on
    """
    np.random.seed(42)
    temp = np.random.rand(len(lat), len(lon), len(time))*100

    ds = xr.Dataset({'temperature': (['lat', 'lon', 'time'],  temp)},
                    coords={'lon': lon,
                            'lat': lat,
                            'time': time})

    return ds


def test_clim_data(clim_data):

    assert not clim_data.temperature.isnull().any()


@pytest.fixture
def make_holes(clim_data, lat):

    N = len(lat) - 1

    tmp = clim_data['temperature'].values
    array = np.random.randint(0, N, 1)
    tmp[array] = np.nan

    clim_data['temperature'].values = tmp

    return clim_data


@pytest.fixture
def weights(lat, lon):

    df = pd.DataFrame()
    lats = np.random.choice(lat, 100)
    lons = np.random.choice(lon, 100)
    df['lat'] = lats
    df['lon'] = lons

    df['areawt'] = np.random.random(len(df['lon']))
    tmp = np.random.random(len(df['lon']))
    tmp[::5] = np.nan
    df['popwt'] = tmp
    df['hierid'] = np.random.choice(np.arange(1, 25), len(lats))

    mapping = {
        h: np.random.choice(np.arange(1, 5)) for h in df['hierid'].values}

    df['ISO'] = [mapping[i] for i in df['hierid']]
    df.index.names = ['reshape_index']

    return df


def test_reindex_spatial_weights(clim_data, weights):

    assert not clim_data.temperature.isnull().any()

    ds = _reindex_spatial_data_to_regions(clim_data, weights)

    assert ds.temperature.shape == (len(ds['lon']), len(ds['time']))
    assert 'reshape_index' in ds.dims


def test_weighting(clim_data, weights):

    assert np.isnan(weights['popwt'].values).any()
    ds = _reindex_spatial_data_to_regions(clim_data, weights)
    assert not ds.temperature.isnull().any()

    wtd = _aggregate_reindexed_data_to_regions(
        ds, 'temperature', 'popwt', 'ISO', weights)

    assert not wtd.temperature.isnull().any()

    wtd = _aggregate_reindexed_data_to_regions(
        ds, 'temperature', 'areawt', 'ISO', weights)

    assert not wtd.temperature.isnull().any()


def test_rename_coords_to_lon_and_lat():
    ds = xr.Dataset(
        coords={'z': [1.20, 2.58], 'long': [156.6, 38.48]})

    ds = rename_coords_to_lon_and_lat(ds)
    coords = ds.coords

    assert 'z' not in coords
    assert coords.z is None
    assert 'lon' in coords and 'long' not in coords


def test_rename_coords_to_lon_and_lat():
    ds = xr.Dataset(
        coords={'latitude': [71.32, 72.58], 'longitude': [156.6, 38.48]})

    ds = rename_coords_to_lon_and_lat(ds)
    coords = ds.coords

    assert 'lat' in coords and 'latitude' not in coords
    assert 'lon' in coords and 'longitude' not in coords


def test_rename_coords_to_longitude_and_latitude():
    ds = xr.Dataset(
        coords={'lat': [71.32, 72.58], 'lon': [156.6, 38.48]})
    ds = rename_coords_to_longitude_and_latitude(ds)
    coords = ds.coords

    assert 'latitude' in coords and 'lat' not in coords
    assert 'longitude' in coords and 'lon' not in coords


def test_rename_coords_to_longitude_and_latitude_with_clim_data(clim_data):
    ds = rename_coords_to_longitude_and_latitude(clim_data)
    coords = ds.coords

    assert 'latitude' in coords and 'lat' not in coords
    assert 'longitude' in coords and 'lon' not in coords


def test_convert_lons_mono():
    ds = xr.Dataset(coords={'lon': [-156.6, -38.48]})
    expected = np.array([203.4, 321.52])

    ds = convert_lons_mono(ds, lon_name='lon')

    np.testing.assert_array_equal(ds.lon.values, expected)


def test_convert_lons_split():
    ds = xr.Dataset(coords={'longitude': [300, 320]})
    expected = np.array([-60, -40])

    ds = convert_lons_split(ds)

    np.testing.assert_array_equal(ds.longitude.values, expected)


def test_remove_leap_days():
    da = xr.DataArray(
        np.random.rand(4, 3),
        [('time', pd.date_range('2000-02-27', periods=4)),
         ('space', ['IA', 'IL', 'IN'])])
    leap_day = np.datetime64('2000-02-29')

    da = remove_leap_days(da)

    assert leap_day not in da.coords['time'].values


def test_remove_leap_days_with_clim_data(clim_data):
    leap_day = np.datetime64('2000-02-29')

    da = remove_leap_days(clim_data)

    assert leap_day not in da.coords['time'].values


def test_convert_kelvin_to_celsius(clim_data):
    ds = convert_kelvin_to_celsius(clim_data, 'temperature')

    assert 'C' in ds.data_vars['temperature'].units


def test_standardize_climate_data(clim_data):
    ds = standardize_climate_data(clim_data)

    coordinates = ds.coords

    assert 'lat' in coordinates and 'latitude' not in coordinates
    assert 'lon' in coordinates and 'longitude' not in coordinates
    
    
def test_snyder_edd():
    ds_tasmax = xr.Dataset(
        data_vars={'tmax': (
            '(latitude, longitude)', [280.4963, 280.7887], {'units': 'K'})},
        coords={'latitude': [-33.625, -33.375], 'longitude': [286.125, 286.375]})

    ds_tasmin = xr.Dataset(
        data_vars={'tmin': (
            '(latitude, longitude)', [278.902, 278.23163], {'units': 'K'})},
        coords={'latitude': [-33.625, -33.375], 'longitude': [286.125, 286.375]})

    threshold = 8

    res = snyder_edd(
                ds_tasmin.tmin,
                ds_tasmax.tmax,
                threshold=273.15 + threshold)
    
    assert res.units == 'degreedays_281.15K'
    assert res.sum().item(0) == 0.0

    
def test_snyder_gdd():
    ds_tasmax = xr.Dataset(
        data_vars={'tmax': ('(latitude, longitude)', [280.4963, 280.7887], {'units': 'K'})},
        coords={'latitude': [-33.625, -33.375], 'longitude': [286.125, 286.375]})

    ds_tasmin = xr.Dataset(
        data_vars={'tmin': ('(latitude, longitude)', [278.902, 278.23163], {'units': 'K'})},
        coords={'latitude': [-33.625, -33.375], 'longitude': [286.125, 286.375]})
    
    res = snyder_gdd(
            ds_tasmin.tmin,
            ds_tasmax.tmax,
            threshold_low=273.15 + 1,
            threshold_high=273.15 + 8)
    
    assert not res.units == 'degreedays_281.15K'
    assert res.units == 'degreedays_274.15-281.15K'
    assert res.sum().item(0) == pytest.approx(11, 0.1)
