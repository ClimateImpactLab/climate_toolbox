'''
This file describes the process for computing weighted climate data
'''
import toolz
import itertools
import xarray as xr
import numpy as np
import pandas as pd
<<<<<<< HEAD
import geopandas as gpd
import datafs
import time
=======
>>>>>>> master
from scipy.ndimage import label
from scipy.interpolate import griddata
import geopandas as gpd
from rasterstats import zonal_stats
from six import string_types
<<<<<<< HEAD
from itertools import ifilterfalse

=======
import itertools
import toolz
import warnings

from distutils.version import LooseVersion
>>>>>>> master

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

    method :

    '''
    if isinstance(broadcast_dims, string_types):
        broadcast_dims = (broadcast_dims, )

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

        iterative_fill_holes(
            da=ds[varname][slicer_dict],
            lat_name=lat_name,
            lon_name=lon_name,
            method=method)


def iterative_fill_holes(da, lat_name='lat', lon_name='lon', method='linear'):
    '''
    Interpolates missing data within a progressively widening bounding box


    Parameters
    ----------
    da: DataArray with dims lat and lon

    lat_name: str

    lon_name: str

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
            warnings.warn(
                'Maximum allowed attempts exceeded in iterative_fill_holes')

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

        while (da.isnull()).isel(**{
                lat_name: imin_lat,
                lon_name: slice(imin_lon, imax_lon)}).any():

            imin_lat = max(imin_lat - 1, 0)
            if imin_lat == 0:
                break

        while (da.isnull()).isel(**{
                lat_name: imax_lat,
                lon_name: slice(imin_lon, imax_lon)}).any():

            imax_lat = min(imax_lat + 1, len(da.coords[lat_name]) - 1)
            if imax_lat == len(da.coords[lat_name]) - 1:
                break

        while (da.isnull()).isel(**{
                lat_name: slice(imin_lat, imax_lat),
                lon_name: imin_lon}).any():

            imin_lon = max(imin_lon - 1, 0)
            if imin_lon == 0:
                break

        while (da.isnull()).isel(**{
                lat_name: slice(imin_lat, imax_lat),
                lon_name: imax_lon}).any():

            imax_lon = min(imax_lon + 1, len(da.lon) - 1)
            if imax_lon == len(da.lon) - 1:
                break

        ravel_lats, ravel_lons = (
            np.meshgrid(
                da.coords[lat_name].values, da.coords[lon_name].values))

        inds_lon, inds_lat = np.where(
                (ravel_lats >= ravel_lats[imin_lat]) &
                (ravel_lats <= ravel_lats[imax_lat]) &
                (ravel_lons >= ravel_lons[imin_lon]) &
                (ravel_lons <= ravel_lons[imax_lon]))

        var_box = var[inds_lat, inds_lon]
        lat_box = ravel_lats[inds_lon, inds_lat]
        lon_box = ravel_lons[inds_lon, inds_lat]
        not_missing = np.where(~var_box.mask)

        points = np.column_stack(
            [lat_box[not_missing].T, lon_box[not_missing].T])

        values = var_box[~var_box.mask]

        da.values[inds_lat, inds_lon] = griddata(
                points,
                values,
                (lat_box, lon_box),
                method=method)


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

    dims = np.array(ds.dims)

    assert len(dims[np.in1d(dims, lon_names)]) == 1
    _lon_coord = dims[np.in1d(dims, ['longitude', 'lon'])][0]

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

    if '_longitude' in ds.dims:
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

 
def _intersect_grid_admin(grid_shp, admin_shp):
    '''
    Generate a GeoDataFrame of intersected segments of two shapefiles 
    
    Parameters
    ----------

    grid_shp: geopandas dataframe
    admin_shp: geopandas dataframe
    n_jobs: num_jobs to run in parallel


    Returns
    -------
    geopandas DataFrame
    geometries and attributes for each pixelated segment of an admin region

    '''

    import geopandas as gpd
    
    #generate an rtree index of gridded climate data
    grid_idx = grid_shp.sindex
    #iterate over admin geometry rows and 
    segments = []
    for feature in admin_shp.itertuples():

        segments.append(_gen_segments(grid_idx=grid_idx, grid_shp=grid_shp, 
                          feature=feature._asdict()))
    flatten = [item for sublist in segments for item in sublist]
    return gpd.GeoDataFrame(flatten, columns=flatten[0].keys())


def _gen_segments(grid_idx, grid_shp, feature):
    '''
    Function applied to each administrative section to determine the gridded segments whose
    boundaries are within or intersect the admin polygon
    
    Parameters
    ----------
    grid_idx: spatial index
        Rtree index of gridded dataset
    grid_shp: shapefile
        gridded shapefile
    feature: OrderedDict
        dictionary of attributes for a admin feature
        
    Returns
    -------
        list
            pixel-level centroids for all segments in an administrative feature
    '''
    t1 = time.time()

    segments = []

    #if geometry is a string, turn it into a shapely geometry object
    if not feature['geometry'].is_valid:
        feature_geom = feature['geometry'].buffer(0)

    else: feature_geom = feature['geometry']

    #get the indexes of climate grid pixe;s that may intersect admin grids
    near_grid_ids = list(grid_idx.intersection(feature_geom.bounds))
    #get the actual grid cells
    near_grid_features = grid_shp.iloc[near_grid_ids]

    #get the list of intersecting pixels
    intersect_list = filter(feature_geom.intersects, near_grid_features.geometry)

    #get list of pixels that intersect partially
    boundary_list = ifilterfalse(feature_geom.contains, intersect_list)

    #get list of pixels that are completely interior
    interior_list = filter(feature_geom.contains, intersect_list)

    #construct a geometry segment for each interior pixel
    for geom in interior_list:
        segments.append(_gen_segment(geom, feature, method='interior'))

    #construct a geometry segment for each intersecting pixel
    for geom in boundary_list:
        segments.append(_gen_segment(geom, feature, method='segment'))

    t2 = time.time()
    print('generate_segments: ', t2 -t1)
    print('number_of_intersecting_segments:', len(near_grid_features))

    return segments

def _gen_segment(grid_geom, admin_feat, method=None):
    '''
    Computes the geometry and set of attributes associated with
    a segment of a shapefile. 

    Parameters
    ----------
    grid: shapely geometry

    admin: shapely geometry

    method: str 
        'segment', 'overlap', 'interior', 'centroid'

    Returns
    -------
        dict
        dictionary for one segment representing a shapely geometry 
        and segment features for pixel centroid 
    '''
    
    if method == 'interior':
        seg_geom = grid_geom
    #finds the geometry of the intersecting segment
    elif method == 'segment':
            seg_geom = admin_feat['geometry'].intersection(grid_geom)
    #set the boundaries of the grid_geom centroid as pix_cent_x and pix_cent_y
        
    properties = {'lon': grid_geom.centroid.bounds[0],
                 'lat': grid_geom.centroid.bounds[1]} 
        
    #need to combine features from admin and geometry of pix
    properties.update({k: v for k, v in admin_feat.items() if k not in ['geometry', 'Index']})
    properties['geometry'] = seg_geom
    
    return properties

def _generate_weights(intersected, raster_path, weighting, header_ids=None):
    '''
    calculates the weighted values for each segment
    Will calculate a area weighted average by default

    Parameters
    ----------

    intersected: GeoDataFramae
        geopandas GeoDataFrame of segmented points of the administrative region

    raster_path: str
        file path to raster data

    weighting: str
       the raster dataset variable you are measuring: 'crop', 'pop', etc
    
    header_ids: list
        list of column headers to set for new dataframe/csv

    Returns
    -------
    pandas DataFrame of segment by weighting 
    '''

    #main utility for computing statistics for each pixel
    stats = zonal_stats(intersected, raster_path)
    
    segs = []

    weighting_list = [weighting, 'area']

    #loop through our values
    for i, seg in enumerate(intersected.itertuples()):
        props = {}
        #this simply generates the ara in the square pixel
        #whatever the grid resolution is the area/pixel = resolution**2
        props['area'] = seg.geometry.area
        if stats[i]['mean'] is not None:
            #compute a total for that segment 
            props[weighting] = stats[i]['count']*stats[i]['mean']
        else: props[weighting] = 0
        
        #these are named tuples so we read them asdict
        props.update({k:v for k, v in seg._asdict().items()})
        segs.append(props)

    #restructure our data a bit
    headers = header_ids + weighting_list + ['lat' ,'lon']
    seg_df = pd.DataFrame(segs, columns=headers, index=range(len(segs)))

    #this provides totals across the regional designations we care about
    #returns the total of the weights
    totals_df = seg_df.groupby(header_ids).transform('sum')[weighting_list]
    #some column relabeling
    totals_df.columns = [wt + 'total' for wt in weighting_list]
    #construct our new dataframe
    df = pd.concat([seg_df, totals_df], axis=1)
    
    #compute the wts per segment
    for wt in weighting_list:
        df[wt+'wt'] = df[wt] / df[wt+'total']

    return df

'''
================
Public Functions
================
'''

def grid_to_region(grid_gdf, admin_gdf, raster_fp, weighting, header_ids=None, intersected_fp=None):
    '''
    Constructs dataframe of weights, feature values, and pixel centroids
    used in computing the weighted climate for a given grid segment

    Parameters
    ----------

    grid_gdf: GeoDataFrame
         climate grid DataFrame

    admin_gdf: GeoDataFrame
        administrative GeoDataFrame

    raster_fp: str
        file path of raster file
    
    weighting: str
        weighting that raster file is measuring i.e. crop, irrigation, population, etc

    header_ids: list
        list of string-formatted column header values for new dataframe

    intersected_fp: str

        file path to read/write intersected GeoDataFrame

    Returns
    -------
        pd.DataFrame of weights by regional feature set
    '''

    #produce the intersected grid_file
    print('intersecting')
    intersected = _intersect_grid_admin(grid_gdf, admin_gdf)

    print('writing intersection file')
    if intersected_fp:
        intersected.to_file(intersected_fp)

    #write to disk, so we can read it in for our weight generation

    #create weights dataframe
    print('generating weights')
    df = _generate_weights(intersected, raster_fp, weighting, header_ids=header_ids)

    return df

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
