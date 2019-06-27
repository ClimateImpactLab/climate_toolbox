import xarray as xr
import numpy as np

from climate_toolbox.utils.utils import remove_leap_days

def edd_ag(ds_tasmax, ds_tasmin, threshold):
    '''
    Note: there are implicitly three cases:

        1. tmax > threshold & tmin < threshold
        2. tmax > threshold & tmin > threshold
        3. tmax <= threshold
    
    Case (1) is the first part of the np.where() statement.
    Case (2) is also the first part of this statement, which returns 'NA' in
    this case and so is not included in the final summation of HDD. (I.e., in
    case (2), the HDD should be zero. Instead we get 'NA', which functions as
    the equivalent of a zero value in the subsequent code.)
    Case (3) is, of course, the second part of the np.where() statement.
    
    Parameters
    ----------
    ds : Dataset
        xarray.Dataset with two variables: tasmin and tasmax. 
        tasmin and tasmax are in Kelvin and are indexed by
        impact region (``hierid``) and day (``time``) in a
        365-day calendar.
        
    Returns
    -------
    ds : Dataset
        xarray.Dataset with dimensions ``(hierid, threshold)``
    '''

    # convert from K to C
    tmax = (ds_tasmax.tasmax - 273.15)
    tmin = (ds_tasmin.tasmin - 273.15)

    #get the 
    snyder_m = (tmax + tmin)/2
    snyder_w = (tmax - tmin)/2
    snyder_theta = np.arcsin( (threshold - snyder_m)/snyder_w )

    transdata = np.where(
        tmin.values < threshold,
            np.where(
                tmax.values > threshold,
            ((snyder_m.values - threshold) * (np.pi/2 - snyder_theta.values) +
                         snyder_w.values * np.cos(snyder_theta.values) ) / np.pi, 0),
        snyder_m.values - threshold)

    
    return xr.DataArray(transdata, dims=tmax.dims, coords=tmax.coords)


def tas_poly(ds, power):
    '''
    Daily average temperature (degrees C), raised to a power

    Leap years are removed before counting days (uses a 365 day
    calendar).
    '''

    powername = ordinal(power)

    description = format_docstr(('''
            Daily average temperature (degrees C){raised}

            Leap years are removed before counting days (uses a 365 day
            calendar).
            '''.format(
                raised='' if power == 1 else (
                    ' raised to the {powername} power'
                    .format(powername=powername)))).strip())

    ds1 = xr.Dataset()

    # remove leap years
    ds = remove_leap_days(ds)

    # do transformation
    ds1[varname] = (ds.tas - 273.15)**power

    # Replace datetime64[ns] 'time' with YYYYDDD int 'day'
    if ds.dims['time'] > 365:
        raise ValueError

    ds1.coords['day'] = ds['time.year']*1000 + np.arange(1, len(ds.time)+1)
    ds1 = ds1.swap_dims({'time': 'day'})
    ds1 = ds1.drop('time')

    ds1 = ds1.rename({'day': 'time'})

    # document variable
    ds1[varname].attrs['units'] = (
        'C^{}'.format(power) if power > 1 else 'C')

    ds1[varname].attrs['long_title'] = description.splitlines()[0]
    ds1[varname].attrs['description'] = description
    ds1[varname].attrs['variable'] = varname

    return ds1