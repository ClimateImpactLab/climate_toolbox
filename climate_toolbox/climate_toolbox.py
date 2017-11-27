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

 
def _intersect_grid_admin(grid_shp, admin_shp, n_jobs=1):
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
    
    from multiprocessing import cpu_count
    from joblib import delayed, Parallel
    import geopandas as gpd
    
    if not n_jobs:
        p  = Parallel(n_jobs=cpu_count(), verbose=2)
        
    else: p = Parallel(n_jobs=n_jobs, verbose=2)
    
    #generate an rtree index of gridded climate data
    grid_idx = grid_shp.sindex
    #iterate over admin geometry rows and 
    segments = []
    segments.append(
                    p(delayed(gen_segments)
                         (grid_idx=grid_idx, grid_shp=grid_shp, 
                          feature=feature._asdict()) for feature in admin_shp.itertuples()))
    
    flattened = [item for sublist in segments for item in sublist]
    
    return gpd.GeoDataFrame(flattened, columns=flattened[0].keys())


def _gen_segments(grid_idx, feature):
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
    from itertools import ifilterfalse


    segments = []

    #if geometry is a string, turn it into a shapely geometry object
    if not feature.geometry.is_valid:
        feature_geom = feature.geometry.buffer(0)

    else: feature_geom = feature.geometry

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
        segments.append(gen_segment(geom, feature._asdict(), method='interior'))

    #construct a geometry segment for each intersecting pixel
    for geom in boundary_list:
        segments.append(gen_segment(geom, feature._asdict(), method='segment'))

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
            seg_geom = admin_feat.geometry.intersection(grid_geom)
    #set the boundaries of the grid_geom centroid as pix_cent_x and pix_cent_y
        
    properties = {'pix_cent_x': grid_geom.centroid.bounds[0],
                 'pix_cent_y': grid_geom.centroid.bounds[1]} 
        
    #need to combine features from admin and geometry of pix
    properties.update({k: v for k, v in admin_feat.items() if k not in ['geometry', 'Index']})
    properties['geometry'] = seg_geom
    
    return properties


def _generate_weights(shapefile_path, raster_path, weighting, header_ids=None):
    '''
    calculates the weighted values for each segment
    Will calculate a area weighted average by default

    Parameters
    ----------

    segment_shapefile_path: str
        file path to segmented_shapefile data

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
    import geopandas as gpd
    import pandas as pd
    from rasterstats import zonal_stats
    

    #main utility for computing statistics for each pixel
    stats = zonal_stats(shapefile_path, raster_path)

    segment_layer = gpd.read_file(shapefile_path)
    
    segs = []

    weighting_list = [weighting, 'area']

    #loop through our values
    for i, seg in enumerate(segment_layer.itertuples()):
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
    headers = header_ids + weighting_list + ['pix_cent_x' ,'pix_cent_y']
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

def create_segment_weights_df(grid_shp, admin_shp, raster_fp, weighting, intersected_fp=None, header_ids=None, n_jobs=None):
    '''
    Constructs dataframe of weights, feature values, and pixel centroids 
    used in computing the weighted climate for a given grid segment

    Parameters
    ----------

    grid_shp: GeoDataFrame
         climate grid DataFrame

    admin_shp: GeoDataFrame
        administrative GeoDataFrame

    raster_fp: str
        file path of raster file

    intersected_fp: str
        file path to read/write intersected GeoDataFrame

    weighting: str
        weighting that raster file is measuring i.e. crop, irrigation, population, etc

    intersection_method: str
        how to evaluate intersections of polygons i.e centroid, segment etc

    header_ids: list
        column header values for new dataframe

    Returns
    -------
        pd.DataFrame of weights by regional feature set
    '''

    #produce the intersected grid_file
    intersected = _intersect_grid_admin(grid_shp, admin_shp, n_jobs=n_jobs)

    if intersected_fp:
        intersected.to_file(intersected_fp)

    #write to disk, so we can read it in for our weight generation

    #create weights dataframe
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
