#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `climate_toolbox` package."""

import pytest

from click.testing import CliRunner

from climate_toolbox import climate_toolbox as ctb
from climate_toolbox import cli

import numpy as np
import pandas as pd
import pytest
import xarray as xr
import itertools


#python utils

#create pytest resource to be used across tests
@pytest.fixture
def clim_data():
    '''
    Generate fake climate data to do tests on
    '''
    np.random.seed(42)
    temp = np.random.rand(720, 1440,365)*100
    lat = np.arange(-89.875,90,0.25)
    lon = np.arange(0.125, 360., .25)
    time = pd.date_range(start=pd.datetime(2000, 1, 1), periods=365)
    ds = xr.Dataset({'temperature': (['lat', 'lon', 'time'],  temp)}, 
                    coords={'lon':  lon,
                            'lat':  lat,
                            'time':time})

    return ds


def test_clim_data(clim_data):

    assert clim_data.isnull().any() == False


@pytest.fixture
def make_holes(clim_data):

    tmp = clim_data['temperature'].values
    array = np.random.randint(0,719, 1)
    tmp[array] = np.nan

    clim_data['temperature'].values = tmp

    return clim_data


#create pytest resource to be used acress tessts
def test_fill_holes(make_holes):

    assert make_holes.isnull().any() == True
    ctb._fill_holes_xr(make_holes, 'temperature')

    assert make_holes.isnull().any() == False

@pytest.fixture
def weights():

    df = pd.DataFrame()
    lat = np.random.choice(np.arange(-89.875,90,0.25), 1000)
    lon = np.random.choice(np.arange(0.125, 360., .25), 1000)
    df['lat'] = lat
    df['lon'] = lon
    
    df['areawt'] = np.random.random(len(df['lon']))
    tmp = np.random.random(len(df['lon']))
    tmp[::5] = np.nan
    df['popwt'] = tmp
    df['hierid'] = np.random.choice(np.arange(1,250), len(lat))
    mapping = {h: np.random.choice(np.arange(1,20)) for h in df['hierid'].values}
    df['ISO'] = [mapping[i] for i in df['hierid']]
    df.index.names = ['reshape_index']

    return df


def test_reindex_spatial_weights(clim_data, weights):

    assert clim_data.isnull().any() == False

    ds = ctb._reindex_spatial_data_to_regions(clim_data, weights)

    assert ds.temperature.shape == (len(ds['lon']), len(ds['time']))
    assert 'reshape_index' in ds.dims




def test_weighting(clim_data, weights):

    assert np.isnan(weights['popwt'].values).any() == True
    ds = ctb._reindex_spatial_data_to_regions(clim_data, weights)
    assert ds.temperature.isnull().any() == False

    wtd = ctb._aggregate_reindexed_data_to_regions(ds, 'temperature', 'popwt', 'ISO', weights)
    assert wtd.temperature.isnull().any() == False

    wtd = ctb._aggregate_reindexed_data_to_regions(ds, 'temperature', 'areawt', 'ISO', weights)
    assert wtd.temperature.isnull().any() == False







@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'climate_toolbox.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output



