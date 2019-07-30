import xarray as xr
import numpy as np

from climate_toolbox.utils.utils import \
    remove_leap_days, convert_kelvin_to_celsius


def snyder_edd(tasmin, tasmax, threshold):
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

    tasmin : xarray.DataArray
        Daily minimum temperature (degrees C)

    tasmax : xarray.DataArray
        Daily maximum temperature (degrees C)

    threshold : int, float, xarray.DataArray
        Threshold (degrees C)

    Returns
    -------

    edd : xarray.DataArray
        Snyder exceedance degree days (degreedays)

    """

    # Check for unit agreement
    assert tasmin.units == tasmax.units

    # check to make sure tasmax > tasmin everywhere
    assert not (tasmax < tasmin).any(), "values encountered where tasmin > tasmax"

    # compute useful quantities for use in the transformation
    snyder_mean = ((tasmax + tasmin)/2)
    snyder_width = ((tasmax - tasmin)/2)
    snyder_theta = xr.ufuncs.arcsin((threshold - snyder_mean)/snyder_width)

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
        'degreedays_{}{}'.format(threshold, tasmax.attrs['units']))

    return res


def snyder_gdd(tasmin, tasmax, threshold_low, threshold_high):
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

    tasmin : xarray.DataArray
        Daily minimum temperature (degrees C)

    tasmax : xarray.DataArray
        Daily maximum temperature (degrees C)

    threshold_low : int, float, xarray.DataArray
        Lower threshold (degrees C)

    threshold_high : int, float, xarray.DataArray
        Upper threshold (degrees C)

    Returns
    -------

    gdd : xarray.DataArray
        Snyder growing degree days (degreedays)

    """

    # Check for unit agreement
    assert tasmin.units == tasmax.units

    res = (
        snyder_edd(tasmin, tasmax, threshold_low)
        - snyder_edd(tasmin, tasmax, threshold_high))

    res.attrs['units'] = (
        'degreedays_{}{}'.format(threshold, tasmax.attrs['units']))

    return res


def validate_edd_snyder_agriculture(ds, thresholds):
    msg_null = 'hierid dims do not match 24378'

    assert ds.hierid.shape == (24378,), msg_null

    for threshold in thresholds:
        assert threshold in list(ds.refTemp)
    return


def tas_poly(ds, power, varname):
    """
    Daily average temperature (degrees C), raised to a power

    Leap years are removed before counting days (uses a 365 day
    calendar).
    """

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


def ordinal(n):
    """ Converts numbers into ordinal strings """

    return (
        "%d%s" %
        (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10::4]))
