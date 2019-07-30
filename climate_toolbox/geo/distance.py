import numpy as np
#             model             major (km)   minor (km)     flattening
ELLIPSOIDS = {'WGS-84':        (6378.137,    6356.7523142,  1 / 298.257223563),
              'GRS-80':        (6378.137,    6356.7523141,  1 / 298.257222101),
              'Airy (1830)':   (6377.563396, 6356.256909,   1 / 299.3249646),
              'Intl 1924':     (6378.388,    6356.911946,   1 / 297.0),
              'Clarke (1880)': (6378.249145, 6356.51486955, 1 / 293.465),
              'GRS-67':        (6378.1600,   6356.774719,   1 / 298.25),
              }


EARTH_RADIUS = 6371.009


def great_circle(ax, ay, bx, by, radius=EARTH_RADIUS):
    """
    calculate the great circle distance (km) between points

    Provide points (ax, ay) and (bx, by) as floats, or as
    vectors. If ax and ay are vectors or arrays of the
    same shape, the element-wise distance will be found
    between points in the vectors/arrays. If ax, ay are
    (Mx1) column vectors and (bx, by) are (1xN) row
    vectors, the vectors will be broadcast using numpy
    broadcasting rules and the distance between each pair
    of points will be returned as an (MxN) matrix.

    Parameters
    -----------
    ax : float or array
        x/long of point a
    ay : float or array
        y/lat of point a
    bx : float or array
        x/long of point b
    by : float or array
        y/lat of point b
    radius : float, optional
        Radius of the sphere on which to calculate the great
        circle distance (default is to use the Earth's radius in
        km, `6371.009`). Values returned will be in units of the
        radius provided.

    Returns
    --------
    distance : float or array
        great circle distance between points a and b. Units will
        match the radius provided (default km)
    """

    lat1, lng1 = np.radians(ay), np.radians(ax)
    lat2, lng2 = np.radians(by), np.radians(bx)

    sin_lat1, cos_lat1 = np.sin(lat1), np.cos(lat1)
    sin_lat2, cos_lat2 = np.sin(lat2), np.cos(lat2)

    delta_lng = lng2 - lng1
    cos_delta_lng, sin_delta_lng = np.cos(delta_lng), np.sin(delta_lng)

    d = np.arctan2(np.sqrt((cos_lat2 * sin_delta_lng) ** 2 +
                   (cos_lat1 * sin_lat2 -
                    sin_lat1 * cos_lat2 * cos_delta_lng) ** 2),
              sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lng)

    return radius * d
