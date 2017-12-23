'''
This file describes the process for computing weighted climate data
'''

import xarray as xr
import numpy as np
import pandas as pd
from scipy.ndimage import label
from scipy.interpolate import griddata
from six import string_types
import itertools
import toolz
import warnings

from distutils.version import LooseVersion

WEIGHTS_FILE = (
    'GCP/spatial/world-combo-new/segment_weights/' +
    'agglomerated-world-new_BCSD_grid_segment_weights_area_pop.csv')
''' filepath to default reshape weights file '''


'''
=================
Private Functions
=================
'''


def _fill_holes_xr(
        ds,
        varname,
        broadcast_dims=('time',),
        lon_name='lon',
        lat_name='lat',
        gridsize=0.25,
        minlat=-85,
        maxlat=85,
        method='linear'):
    '''
    Fill NA values inplace in a gridded dataset

    Parameters
    ----------

    ds : xarray.Dataset
        name of the dataset with variable to be modified

    varname : str
        name of the variable to be interpolated

    broadcast_dims : tuple of strings, optional
        tuple of dimension names to broadcast the interpolation step over
        (default 'time')

    lon_name : str, optional
        name of the longitude dimension (default 'lon')

    lat_name : str, optional
        name of the latitude dimension (default 'lat')

    gridsize : float, optional
        size of the lat/lon grid. Important for creating a bounding box around
        NaN regions (default 0.25)

    minlat : float, optional
        latitude below which no values will be interpolated (default -85)

    minlon : float, optional
        latitude above which no values will be interpolated (default 85)

    '''
    if isinstance(broadcast_dims, string_types):
        broadcast_dims = (broadcast_dims, )

    ravel_lons, ravel_lats = (
        np.meshgrid(ds.coords[lon_name].values, ds.coords[lat_name].values))

    # remove infinite values
    ds[varname] = (
        ds[varname]
        .where(~np.isinf(ds[varname]))
        .where(ds[varname] < 1e10))

    for indexers in itertools.product(*tuple(
            [range(len(ds.coords[c])) for c in broadcast_dims])):

        slicer_dict = dict(zip(broadcast_dims, indexers))

        slicer = tuple([
                slicer_dict[c]
                if c in broadcast_dims
                else slice(None, None, None)
                for c in ds[varname].dims])

        sliced = ds[varname].values.__getitem__(slicer)

        if not np.isnan(sliced).any():
            continue

        # filled = _fill_holes(
        #     var=np.ma.masked_invalid(sliced),
        #     lat=ravel_lats,
        #     lon=ravel_lons,
        #     gridsize=0.25,
        #     minlat=-85,
        #     maxlat=85)

        iterative_fill_holes(
            da=ds[varname][slicer_dict],
            method=method)


def iterative_fill_holes(da, lat='lat', lon='lon', method='linear'):
    '''
    Interpolates missing gridded data using a progressively widening bounding box


    Parameters
    ----------
    da: DataArray with dims lat and lon

    lat: str

    lon: str

    method: str
        options include 'cubic' 1D and 2D, 'linear', 'nearest'

    Returns
        DataArray
    '''

    attempts = 0
    max_attempts = 10

    while da.isnull().any():
        
        attempts += 1
        
        if attempts > max_attempts:
            warnings.warn('Maximum allowed attempts exceeded in iterative_fill_holes')
            break

        var = np.ma.masked_invalid(da.values)

        missing = np.where((var.mask))
        mp = np.zeros(var.shape)
        mp[missing] = 1

        ptch, n_ptch = label(mp)

        if n_ptch * 2 > max_attempts:
            max_attempts = n_ptch * 2

        ipatch = np.random.choice(range(n_ptch))
            
        ilat = np.where((ptch == ipatch).any(axis=1))[0]
        imin_lat = min(ilat)
        imax_lat = max(ilat)
        
        ilon = np.where((ptch == ipatch).any(axis=0))[0]
        imin_lon = min(ilon)
        imax_lon = max(ilon)

        while (da.isnull()).isel(lat=imin_lat, lon=slice(imin_lon, imax_lon)).any():
            imin_lat = max(imin_lat - 1, 0)
            if imin_lat == 0:
                break

        while (da.isnull()).isel(lat=imax_lat, lon=slice(imin_lon, imax_lon)).any():
            imax_lat = min(imax_lat + 1, len(da.lat) - 1)
            if imax_lat == len(da.lat) - 1:
                break

        while (da.isnull()).isel(lat=slice(imin_lat, imax_lat), lon=imin_lon).any():
            imin_lon = max(imin_lon - 1, 0)
            if imin_lon == 0:
                break

        while (da.isnull()).isel(lat=slice(imin_lat, imax_lat), lon=imax_lon).any():
            imax_lon = min(imax_lon + 1, len(da.lon) - 1)
            if imax_lon == len(da.lon) - 1:
                break

        ravel_lats, ravel_lons = (
            np.meshgrid(da.coords[lat].values, da.coords[lon].values))

        inds_lon, inds_lat = np.where(
                (ravel_lats >= ravel_lats[imin_lat]) &
                (ravel_lats <= ravel_lats[imax_lat]) &
                (ravel_lons >= ravel_lons[imin_lon]) &
                (ravel_lons <= ravel_lons[imax_lon]))

        var_box = var[inds_lat, inds_lon]
        lat_box = ravel_lats[inds_lon, inds_lat]
        lon_box = ravel_lons[inds_lon, inds_lat]
        not_missing = np.where(~var_box.mask)
        points = np.column_stack([lat_box[not_missing].T, lon_box[not_missing].T])
        values = var_box[~var_box.mask]

        da.values[inds_lat, inds_lon] = griddata(
                points,
                values,
                (lat_box, lon_box),
                method=method)


def _fill_holes(var, lat, lon, gridsize=0.25, minlat=-85, maxlat=85):
    '''
    Interpolates the missing values between points on grid

    Parameters
    ----------
    var: masked np.array
        array of climate values

    lat: masked np.array
        array of latitude values

    lon: masked np. array
        array of longitude values

    gridsize: int
        corresponds to degrees on the grid for climate data

    minlat: int
        corresponds to min latitude values to include. Used to remove poles

    maxlat: int
        corresponds to max lat values to include. Used to remove poles


    '''
    # fill the missing value regions by linear interpolation
    # pass if no missing values
    if not np.ma.is_masked(var):
        return var
    # or the missing values are only in polar regions
    if not var.mask[20: -20, :].any():
        return var

    # fill the holes
    var_filled = var[:]
    missing = np.where((var.mask) & (lat > minlat) & (lat < maxlat))
    mp = np.zeros(var.shape)
    mp[missing] = 1
    ptch, n_ptch = label(mp)

    for p in range(1, n_ptch+1):

        ind_ptch = np.where(ptch == p)
        lat_ptch = lat[ind_ptch]
        lon_ptch = lon[ind_ptch]

        ind_box = np.where(
                (lat <= np.max(lat_ptch)+gridsize) &
                (lat >= np.min(lat_ptch)-gridsize) &
                (lon <= np.max(lon_ptch)+gridsize) &
                (lon >= np.min(lon_ptch)-gridsize))

        var_box = var[ind_box]
        lat_box = lat[ind_box]
        lon_box = lon[ind_box]
        not_missing = np.where(~var_box.mask)
        points = np.column_stack([lon_box[not_missing], lat_box[not_missing]])
        values = var_box[~var_box.mask]
        var_filled[ind_box] = griddata(
                points,
                values,
                (lon_box, lat_box),
                method='linear')

    return var_filled


def _standardize_longitude_dimension(ds, lon_names=['lon', 'longitude']):
    '''
    Rescales the lat and lon coordinates to ensure lat is within (-90,90)
    and lon is within (-180, 180). Renames coordinates
    from lon to longitude and from lat to latitude. Sorts any new
    rescaled coordinated.

    Parameters
    ----------
    ds: xarray.DataSet

    Returns
    -------
    ds: xarray.DataSet

    .. note:: this will be unnecessary if we standardize inputs. We can
    scale the longitude dim to between (-180, 180)

    '''

    coords = np.array(ds.coords.keys())

    assert len(coords[np.in1d(coords, lon_names)]) == 1
    _lon_coord = coords[np.in1d(coords, ['longitude', 'lon'])][0]

    ds = ds.rename({_lon_coord: '_longitude'})

    # Adjust lat and lon to make sure they are within (-90, 90) and (-180, 180)
    ds['_longitude_adjusted'] = (
        (ds._longitude - 360)
        .where(ds._longitude > 180)
        .fillna(ds._longitude))

    # reassign the new coords to as the main lon coords
    ds = (
        ds
        .swap_dims({'_longitude': '_longitude_adjusted'})
        .reindex({'_longitude_adjusted': sorted(ds._longitude_adjusted)}))

    if '_longitude' in ds.coords:
        ds = ds.drop('_longitude')

    ds = ds.rename({'_longitude_adjusted': _lon_coord})

    return ds


@toolz.memoize
def _prepare_spatial_weights_data(weights_file=None):
    '''
    Rescales the pix_cent_x colum values

    Requires the :py:mod:`datafs` package.

    Parameters
    ----------
    weights_file: str
        location of file used for weighting


    .. note:: unnecessary if we can standardize our input
    '''

    import datafs

    if weights_file is None:
        weights_file = WEIGHTS_FILE

        api = datafs.get_api()
        archive = api.get_archive(weights_file)

        with archive.open('r') as f:
            df = pd.read_csv(f)
    else:
        df = pd.read_csv(weights_file)

    # Re-label out-of-bounds pixel centers
    df.set_value((df['pix_cent_x'] == 180.125), 'pix_cent_x', -179.875)

    # probably totally unnecessary
    df.drop_duplicates()
    df.index.names = ['reshape_index']

    df.rename(
        columns={'pix_cent_x': 'lon', 'pix_cent_y': 'lat'},
        inplace=True)

    return df


def _reindex_spatial_data_to_regions(ds, df):
    '''
    Reindexes spatial and segment weight data to regions
    Enables region index-based math operations
    Parameters
    ----------
    ds: xarray Dataset
    df: pandas DataFrame
    Returns
    -------
    Xarray DataArray
    '''

    # use vectorized indexing in xarray >= 0.10
    if LooseVersion(xr.__version__) > LooseVersion('0.9.999'):

        lon_indexer = xr.DataArray(df.lon.values, dims=('reshape_index', ))
        lat_indexer = xr.DataArray(df.lat.values, dims=('reshape_index', ))

        return ds.sel(lon=lon_indexer, lat=lat_indexer)

    else:
        res = ds.sel_points(
            'reshape_index',
            lat=df.lat.values,
            lon=df.lon.values)

        return res


def _aggregate_reindexed_data_to_regions(
        ds,
        variable,
        aggwt,
        agglev,
        weights,
        backup_aggwt='areawt'):
    '''
    Performs weighted avg for climate variable by region

    Parameters
    ----------

    ds: xarray.DataArray

    variable: str
        name of the data variable

    aggwt: str
        variable to weight by (i.e popwt, areawt, cropwt)

    agglev: str
        indicates which regional id scheme to select in the dataframe

    weight: pd.DataFrame
        pandas DataFrame of weights

    backup_aggwt: str, optional
        aggregation weight to use in regions with no aggwt data (default
        'areawt')

    '''

    ds.coords[agglev] = xr.DataArray(
                weights[agglev].values,
                dims={'reshape_index': weights.index.values})

    # format weights
    ds[aggwt] = xr.DataArray(
                weights[aggwt].values,
                dims={'reshape_index': weights.index.values})

    ds[aggwt] = (
        ds[aggwt]
        .where(ds[aggwt] > 0)
        .fillna(weights[backup_aggwt].values))

    weighted = xr.Dataset({
        variable: (
            (
                (ds[variable]*ds[aggwt])
                .groupby(agglev)
                .sum(dim='reshape_index')) /
            (
                ds[aggwt]
                .groupby(agglev)
                .sum(dim='reshape_index')))})

    return weighted


'''
================
Public Functions
================
'''


def load_bcsd(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    '''
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
    '''

    if lon_name is not None:
        lon_names = [lon_name]

    if hasattr(fp, 'sel_points'):
        ds = fp

    else:
        with xr.open_dataset(fp) as ds:
            ds.load()

    _fill_holes_xr(ds, varname, broadcast_dims=broadcast_dims)
    return _standardize_longitude_dimension(ds, lon_names=lon_names)


def load_baseline(fp, varname, lon_name='lon', broadcast_dims=None):
    '''
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
    '''

    if lon_name is not None:
        lon_names = [lon_name]

    if broadcast_dims is None:
        broadcast_dims = tuple([])

    if hasattr(fp, 'sel_points'):
        ds = fp

    else:
        with xr.open_dataset(fp) as ds:
            ds.load()

    if 'lat' in ds.data_vars:
        ds = ds.set_coords('lat')
        ds = ds.swap_dims({'nlat': 'lat'})

    if 'lon' in ds.data_vars:
        ds = ds.set_coords('lon')
        ds = ds.swap_dims({'nlon': 'lon'})

    _fill_holes_xr(ds, varname, broadcast_dims=broadcast_dims)
    return _standardize_longitude_dimension(ds, lon_names=lon_names)


def weighted_aggregate_grid_to_regions(
        ds,
        variable,
        aggwt,
        agglev,
        weights=None):
    '''
    Computes the weighted reshape of gridded data

    Parameters
    ----------
    ds : xr.Dataset
        xarray Dataset to be aggregated. Must have 'lat' and 'lon' in the
        coordinates.

    variable : str
        name of the variable to be aggregated

    aggwt : str
        Weighting variable (e.g. 'popwt', 'areawt'). This must be a column name
        in the weights file.

    agglev : str
        Target regional aggregation level (e.g. 'ISO', 'hierid'). This must be
        a column name in the weights file.

    weights : str, optional
        Regional aggregation weights (default agglomerated-world-new BCSD
        segment weights)

    Returns
    -------
    ds: xr.Dataset
        weighted and averaged dataset based on agglev
    '''

    if weights is None:
        weights = _prepare_spatial_weights_data()

    ds = _reindex_spatial_data_to_regions(ds, weights)
    ds = _aggregate_reindexed_data_to_regions(
        ds,
        variable,
        aggwt,
        agglev,
        weights)

    return ds
