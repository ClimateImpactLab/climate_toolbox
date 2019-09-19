import xarray as xr
import numpy as np

from climate_toolbox.utils.utils import \
    remove_leap_days, convert_kelvin_to_celsius


def snyder_edd(tasmin_tasmax, threshold, **unused_kwargs):
    r"""
    Snyder exceedance degree days/cooling degree days

    Similarly to Snyder HDDs, Snyder exceedance degree days for any given day
    are given by the integral between the sinosiod-interpolated temperature and
    the threshold.

    The closed form solution is given by:

    .. math::

        EDD_{P} = \sum_{d \in P} EDD_d

    where

    .. math::

        EED_d =
            \begin{cases}
                ( (M - e)(\pi /2 - \theta) + w \cos(\theta) ) / \pi, & \text{if } tmin_d < e < tmax_d \\
                0 , & \text{if } tmax_d < e \\
                M - e, & \text{otherwise}
            \end{cases}

    and

    .. math::

        \begin{array}{rll}
            M & = & (tmax_d + tmin_d)/2 \\
            w & = & (tmax_d-tmin_d)/2 \\
            \theta & = & \arcsin( (e-M)/w ) \\
        \end{array}

    Parameters
    ----------

    tasmin_tasmax : tuple of xarray.Datasets
        tuple containing two :py:class:`xarray.DataArray` objects, in the order
        ``(tasmin, tasmax)``. tasmin should contain a variable ``tasmin`` with
        daily minimum surface air temperature and tasmax should contain a
        variable ``tasmax`` with daily minimum surface temperature. The units
        should be the same as the threshold.

    threshold : int, float, xarray.DataArray
        Threshold (degrees C)

    Returns
    -------

    edd : xarray.DataArray
        Snyder exceedance degree days (degreedays)

    """

    tasmin = tasmin_tasmax[0].tasmin
    tasmax = tasmin_tasmax[1].tasmax

    # Check for unit agreement
    assert tasmin.units == tasmax.units

    # check to make sure tasmax > tasmin everywhere
    assert not (tasmax < tasmin).any(), "values encountered where tasmin > tasmax"

    # compute useful quantities for use in the transformation
    snyder_mean = ((tasmax + tasmin)/2)
    snyder_width = ((tasmax - tasmin)/2)
    snyder_theta = np.arcsin((threshold - snyder_mean)/snyder_width)

    # the trasnformation is computed using numpy arrays, taking advantage of
    # numpy's second where clause. Note that in the current dev build of
    # xarray, xr.where allows this functionality. As soon as this goes live,
    # this block can be replaced with xarray
    res = xr.where(
        tasmin < threshold,
        xr.where(
            tasmax > threshold,
            ((snyder_mean - threshold) * (np.pi/2 - snyder_theta)
                + (snyder_width * np.cos(snyder_theta))) / np.pi,
            0),
        snyder_mean - threshold)

    res.attrs['units'] = (
        'degreedays_{}{}'.format(threshold, tasmax.attrs.get('units', '')))

    out_ds = xr.Dataset({'edd': res})

    out_ds['edd'].attrs['long_title'] = 'Exceedance degree days'
    out_ds['edd'].attrs['description'] = (
        'Exceedance degree days computed using the synder diurnal cycle '
        'interpolation method, computed above a threshold of '
        '{threshold}{units}.'
        .format(
            threshold=threshold,
            units=tasmax.attrs.get('units', '')))

    return out_ds


def snyder_gdd(tasmin_tasmax, threshold_low, threshold_high, **unused_kwargs):
    r"""
    Snyder growing degree days

    Growing degree days are the difference between EDD measures at two
    thresholds.

    .. math::

        {GDD}_{T_{low}, T_{high}, y, i} = {EDD}_{T_{low}, y, i} - {EDD}_{T_{high}, y, i}

    Note that where :math:`tas_{d,i}>{T_{high}}`, GDD will be a constant value
    :math:`T_{high}-T_{low}`. Thus, this measure is only useful when another
    measure, e.g. :math:`{EDD}_{T_{high}}`, sometimes referred to as
    *killing degree days*, is used as an additional predictor.


    Parameters
    ----------

    tasmin_tasmax : tuple of xarray.DataArrays
        tuple containing two :py:class:`xarray.DataArray` objects, in the order
        ``(tasmin, tasmax)``. tasmin should contain daily minimum surface air
        temperature and tasmax should contain daily minimum surface air
        temperature. The units should be the same as the thresholds.

    threshold_low : int, float, xarray.DataArray
        Lower threshold (same units as datasets in tasmin_tasmax)

    threshold_high : int, float, xarray.DataArray
        Upper threshold (same units as datasets in tasmin_tasmax)

    Returns
    -------

    gdd : xarray.DataArray
        Snyder growing degree days (degreedays)

    """

    tasmin, tasmax = tasmin_tasmax

    # Check for unit agreement
    assert tasmin.tasmin.units == tasmax.tasmax.units

    res = (
        snyder_edd((tasmin, tasmax), threshold_low).edd
        - snyder_edd((tasmin, tasmax), threshold_high).edd)

    res.attrs['units'] = (
        'degreedays_{}-{}{}'
        .format(threshold_low, threshold_high, tasmax.tasmax.attrs.get('units', '')))

    out_ds = xr.Dataset({'gdd': res})

    out_ds['gdd'].attrs['long_title'] = 'Growing degree days'
    out_ds['gdd'].attrs['description'] = (
        'Growing degree days computed using the synder diurnal cycle '
        'interpolation method, computed between thresholds of '
        '{threshold_low}{units} and {threshold_high}{units}.'
        .format(
            threshold_high=threshold_high,
            threshold_low=threshold_low,
            units=tasmax.tasmax.attrs.get('units', '')))

    return out_ds


def validate_snyder_edd(ds, thresholds, assert_no_nans=False, **unused_kwargs):
    '''
    TODO:
        This is a CIL-specific validation function. we need a way of
        standardizing this.
    '''

    msg_null = 'hierid dims do not match 24378'

    assert ds.hierid.shape == (24378,), msg_null

    for threshold in thresholds:
        assert threshold in list(ds.refTemp)

    if assert_no_nans:
        for v in ds.data_vars.keys():
            assert ds[v].notnull().all(), "NaNs encountered in {}".format(v)

    return


def tas_poly(ds, power, variable, **unused_kwargs):
    """
    Daily average temperature (degrees C), raised to a power

    Leap years are removed before counting days (uses a 365 day
    calendar).
    """

    description = ('''
            Daily average temperature (degrees C){raised}

            Leap years are removed before counting days (uses a 365 day
            calendar).
            '''.format(
                raised='' if power == 1 else (
                    ' raised to the {powername} power'
                    .format(powername=powername)))).strip()

    out_ds = xr.Dataset()

    # remove leap years
    ds = remove_leap_days(ds)

    # do transformation
    out_ds[variable] = (ds.tas)**power

    # ======================== MOVE THE FOLLOWING TO CIL WRAPPER ==============
    # This is the worst and is impactlab specific until brewster can update the
    # projection system to handle datetime objects. we should put this in
    # climate_transform_specs as a custom writer.

    # Replace datetime64[ns] 'time' with YYYYDDD int 'day'
    if ds.dims['time'] > 365:
        raise ValueError
    out_ds.coords['day'] = ds['time.year']*1000 + np.arange(1, len(ds.time)+1)
    out_ds = out_ds.swap_dims({'time': 'day'})
    out_ds = out_ds.drop('time')
    out_ds = out_ds.rename({'day': 'time'})
    # ======================== MOVE THE ABOVE TO CIL WRAPPER ==================

    # document variable
    out_ds[variable].attrs['units'] = (
        'C^{}'.format(power) if power > 1 else 'C')

    out_ds[variable].attrs['long_title'] = description.splitlines()[0]
    out_ds[variable].attrs['description'] = description
    out_ds[variable].attrs['variable'] = variable

    return out_ds


def validate_tas_poly(
        ds, power, valid_tas_range=(-20, 55), assert_no_nans=True):

    mint = ds.tas.sel(lat=slice(-60, 70)).min()
    maxt = ds.tas.sel(lat=slice(-60, 70)).max()

    min_msg = "min value {} outside allowed range".format(mint)
    assert mint >= min(0, valid_tas_range[0]**power), min_msg

    max_msg = "max value {} outside allowed range".format(mint)
    assert maxt <= (valid_tas_range[1]**power), max_msg

    if assert_no_nans:
        assert ds.tas.notnull().all(), "NaNs found in the data"


def ordinal(n):
    """ Converts numbers into ordinal strings """

    return (
        "%d%s" %
        (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10::4]))
