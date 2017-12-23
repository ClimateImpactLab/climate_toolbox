#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `climate_toolbox` package."""

import pytest

from climate_toolbox import climate_toolbox as ctb

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
    '''
    Generate fake climate data to do tests on
    '''
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


# create pytest resource to be used acress tessts
def test_fill_holes(make_holes):

    assert make_holes.temperature.isnull().any()
    ctb._fill_holes_xr(make_holes, 'temperature', minlat=-90, maxlat=90)

    assert not make_holes.temperature.isnull().any()


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

    ds = ctb._reindex_spatial_data_to_regions(clim_data, weights)

    assert ds.temperature.shape == (len(ds['lon']), len(ds['time']))
    assert 'reshape_index' in ds.dims


def test_weighting(clim_data, weights):

    assert np.isnan(weights['popwt'].values).any()
    ds = ctb._reindex_spatial_data_to_regions(clim_data, weights)
    assert not ds.temperature.isnull().any()

    wtd = ctb._aggregate_reindexed_data_to_regions(
        ds, 'temperature', 'popwt', 'ISO', weights)

    assert not wtd.temperature.isnull().any()

    wtd = ctb._aggregate_reindexed_data_to_regions(
        ds, 'temperature', 'areawt', 'ISO', weights)

    assert not wtd.temperature.isnull().any()
