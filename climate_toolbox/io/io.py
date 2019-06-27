import xarray as xr

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


def load_gmfd(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    pass


def load_best(fp, varname, lon_name='lon', broadcast_dims=('time',)):
    pass
