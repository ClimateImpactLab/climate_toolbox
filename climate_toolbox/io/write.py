
import os

def netcdf_writer(ds, path, **unused_kwargs):
    ''' write to a NetCDF dataset '''

    os.makedirs(os.path.dirname(write_file), exist_ok=True)

    # Coerce attributes to string (required by NetCDF spec)
    ds = ds.copy(deep=False)
    ds.attrs.update({k: str(v) for k, v in ds.attrs.items()})
    for k, v in ds.data_vars.items():
        v.attrs.update({k: str(v) for k, v in v.attrs.items()})

    ds.to_netcdf(path)


def netcdf_write_to_pattern(ds, pattern, **pattern_kwargs):
    netcdf_writer(ds, pattern.format(**pattern_kwargs))


def atomic_netcdf_writer(ds, path, metadata, temp_path=None, **unused_kwargs):
    '''
    atomic write to a NetCDF dataset

    Safe version of climate_toolbox.io.write.netcdf_writer, in which the
    Dataset will be written to a temporary location and then moved into the
    final destination on write completion. This ensures a complete write
    when working in unstable or preemptable computing environments.

    Note that the user is responsible for checking to ensure all writes were
    successful and for deleting any incomplete writes at the temp_path
    location.

    Parameters
    ----------
    ds : Dataset
        :py:class:`xarray.Dataset` object to write
    path : str
        Filepath of write destination
    temp_path : str, optional
        Filepath of temporary write destination. Default uses `path + '~'`.

    '''

    if temp_path is None:
        temp_path = path + '~'

    # Coerce attributes to string (required by NetCDF spec)
    ds = ds.copy(deep=False)
    ds.attrs.update({k: str(v) for k, v in ds.attrs.items()})
    for k, v in ds.data_vars.items():
        v.attrs.update({k: str(v) for k, v in v.attrs.items()})

    ds.to_netcdf(temp_path)
    os.rename(temp_path, path)


def atomic_netcdf_write_to_pattern(ds, pattern, **pattern_kwargs):
    atomic_netcdf_writer(ds, pattern.format(**pattern_kwargs))


def fgh_header_writer(ds, path, **unused_kwargs):

    metacsv.to_header(
        path,
        attrs=dict(ds.attrs),
        variables={v: dict(ds[v].attrs) for v in ds.data_vars.keys()})


def fgh_header_write_to_pattern(ds, pattern, **pattern_kwargs):
    fgh_header_writer(ds, pattern.format(**pattern_kwargs))


def csv_writer(ds, path, metadata, **unused_kwargs):
    raise NotImplementedError
