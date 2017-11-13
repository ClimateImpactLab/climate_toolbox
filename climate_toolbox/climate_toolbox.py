'''
This file describes the process for computing weighted climate data
'''
import toolz
import itertools
import xarray as xr
import numpy as np
import pandas as pd
import datafs
from scipy.ndimage import label
from scipy.interpolate import griddata
from shapely.geometry import shape, mapping
from shapely.prepared import prep
from shapely import speedups
import geopandas as gpd
from rtree.index import Index
from rasterstats import zonal_stats
from six import string_types

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
        maxlat=85):
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

        filled = _fill_holes(
            var=np.ma.masked_invalid(sliced),
            lat=ravel_lats,
            lon=ravel_lons,
            gridsize=0.25,
            minlat=-85,
            maxlat=85)

        ds[varname][slicer_dict] = filled


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

    Parameters
    ----------
    weights_file: str
        location of file used for weighting


    .. note:: unnecessary if we can standardize our input
    '''

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


def intersect_grid_admin(grid_shp, admin_shp, method):
    '''
    Finds intersection point between two shapefiles, 
    constructing a new segment-based shapefile with 
    polygons and centroids for each new segment
    
    Parameters
    ----------

    shape1: geopandas dataframe
    shape2: geopandas dataframe
    method: str


    Returns
    -------
    geopandas dataframe
    '''

    #iterate over each segment in shape2 and find the intersection
    #points between it and shape1
    #Evaluate if the segment is an intersection or boundary

    from joblib import Parallel, delayed
    import multiprocessing


    grid_idx = grid_shp.sindex

    if n_jobs == None:
        n_jobs = multiprocessing.cpu_count()

    with Parallel(n_jobs=n_jobs) as paralllelize:
        segment_features = []
        for i, feature in enumerate(admin_shp.itertuples(), 1)  :
            if not feature.geometry.is_valid:
                feature_geom = feature.geometry.buffer(0)

            feature_geom = feature.geometry

            near_grid_ids = list(grid_idx.intersection(feature_geom.bounds))

            #get intersecting pixels
            intersect_list = filter(admin_prepped.intersects, near_grid_gids)
            # then, get list of pixels that are completely interior
            interior_list = filter(admin_prepped.contains, intersect_list)


            new_segments = []
            admin_segment = paralllelize(delayed(generate_segment)(grid_geom, feature, method=interior) for grid_geom in interior_list)



  
def intersect_grid_admin(grid_shp, admin_shp, method=None):

    '''
    Finds intersection point between two shapefiles, 
    constructing a new segment-based shapefile with 
    polygons and centroids for each new segment
    
    Parameters
    ----------

    shape1: geopandas dataframe
    shape2: geopandas dataframe
    method: str


    Returns
    -------
    geopandas dataframe
    '''

    #iterate over each segment in shape2 and find the intersection
    #points between it and shape1
    #Evaluate if the segment is an intersection or boundary
    
    
    #create an rtree index 
    grid_idx = grid_shp.sindex
    segments = []

    #since we are dealing with a pandas dataframe, we can use itertuples for fast iteration
    for i, feature in enumerate(admin_shp.itertuples(), 1):
        #print(feature)
        #Passing a distance of 0 buffer cleans self-touching or self-crossing polygons
        if not feature.geometry.is_valid:
            feature_geom = feature.geometry.buffer(0)

        else: feature_geom = feature.geometry
        #get the indexes of pixes that may intersect
        near_grid_ids = list(grid_idx.intersection(feature_geom.bounds))
        #print(near_grid_ids)
        #get the actual grid cells
        near_grid_geoms = grid_shp.iloc[near_grid_ids]
        #print(type(near_grid_geoms), len(near_grid_geoms))
        #print(type(feature_geom))
        intersect_list = filter(feature_geom.intersects, near_grid_geoms.geometry)
        boundary_list = ifilterfalse(feature_geom.contains, intersect_list)

        interior_list = filter(feature_geom.contains, intersect_list)
        
        
        new_segments = [] 
        for geom in interior_list:
            print('gen_segment: interior')
            new_segments.append(gen_segment(geom, feature, method='interior'))
            
        for geom in boundary_list:
            print('gen_segment: boundary')

            new_segments.append(gen_segment(geom, feature, method=method))
            
        
        segments.append(new_segments)
        
        
    flattened = [item for sublist in segments for item in sublist]
    cols = flattened[0].keys()
    segment_df = pd.GeoDataFrame(flattened, columns=cols, index=range(len(flattened)))

    return segment_df


def generate_segment(grid_geom, admin_geom, method=None):
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
        dictionary for one segment representing a shapely geometry and segment features
        such as key,values for pixel centroid 
    '''
    

    if method == 'interior' || 'overlap':
        seg_geom = grid_geom
    elif method == 'segment':
            seg_geom = admin_geom.intersection(grid_geom)
    elif method == 'centroid':
        if admin_geom.contains(grid_geom.centroid):
            seg_geom = grid_geom
        else:
            return None
    if seg_geom.geom_type in ['Polygon', 'MultiPolygon']:
    
        props = {'pix_cent_x': grid_geom.centroid.bounds[0],
                 'pix_cent_y': grid_geom.centroid.bounds[1]}

    props.update({k: v for k, v in admin_feat._asdict().items() if k != 'geometry'})

    props['geometry'] = seg_geom

    return props



def generate_weights(shapefile_path, raster_path, weighting):
    '''
    calculates the weighted values for each segment

    Parameters
    ----------

    shapefile_path: str
        file path to shapefile data

    raster_path: str
        file path to raster data

    weighting

    Returns
    -------
    pandas DataFrame of segment by weighting 
    '''
    import geopandas as gpd
    from rasterstats import zonal_stats

    df = pd.DataFrame
    stats = zonal_stats(shapefile_path, raster_path)

    segment_layer = gpd.read_file(shapefile_path)
    
    for i, layer in enumerate(segment_layer):
        headers = {}
        headers['area'] = layer['geometry'].area
        if stats[i]['mean'] is not None:
            headers[weighting] = stats[i]['count'] * stats[i]['mean']

        else:
            headers[weighting] = 0























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
