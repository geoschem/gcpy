"""
Functions to apply grid-stretching transforms.
"""
import numpy as np


def rotate_vectors(x, y, z, k, theta):
    """
    Rotates 3D vectors about an arbitrary axis using Rodrigues'
    rotation formula.

    Parameters
    ----------
    x : array-like
        X-components of the vectors to be rotated.
    y : array-like
        Y-components of the vectors to be rotated.
    z : array-like
        Z-components of the vectors to be rotated.
    k : numpy.ndarray
        Unit vector defining the axis of rotation, shape (3,).
    theta : float
        Angle of rotation in radians.

    Returns
    -------
    x : numpy.ndarray
        X-components of the rotated vectors.
    y : numpy.ndarray
        Y-components of the rotated vectors.
    z : numpy.ndarray
        Z-components of the rotated vectors.
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    v = np.moveaxis(np.array([x, y, z]), 0, -1)  # shape: (..., 3)
    v = v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + \
        k[np.newaxis, :] * np.dot(v, k)[:, np.newaxis] * (1 - np.cos(theta))
    return v[..., 0], v[..., 1], v[..., 2]


def cartesian_to_spherical(x, y, z):
    """
    Converts Cartesian coordinates to spherical (longitude, latitude)
    coordinates in radians.

    Parameters
    ----------
    x : array-like
        X-components of the Cartesian coordinates.
    y : array-like
        Y-components of the Cartesian coordinates.
    z : array-like
        Z-components of the Cartesian coordinates.

    Returns
    -------
    x_sph : numpy.ndarray
        Longitude in radians, in the range [-pi, pi].
    y_sph : numpy.ndarray
        Latitude in radians, in the range [-pi/2, pi/2].
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    # Calculate x,y in spherical coordinates
    y_sph = np.arcsin(z)
    x_sph = np.arctan2(y, x)
    return x_sph, y_sph


def spherical_to_cartesian(x, y):
    """
    Converts spherical (longitude, latitude) coordinates in radians
    to Cartesian coordinates.

    Parameters
    ----------
    x : array-like
        Longitude in radians.
    y : array-like
        Latitude in radians.

    Returns
    -------
    x_car : numpy.ndarray
        X-components of the Cartesian coordinates.
    y_car : numpy.ndarray
        Y-components of the Cartesian coordinates.
    z_car : numpy.ndarray
        Z-components of the Cartesian coordinates.
    """
    x_car = np.cos(y) * np.cos(x)
    y_car = np.cos(y) * np.sin(x)
    z_car = np.sin(y)
    return x_car, y_car, z_car


def schmidt_transform(x, y, s):
    """
    Applies the Schmidt transform to stretch a spherical grid toward
    a pole.

    Parameters
    ----------
    x : array-like
        Longitude in radians.
    y : array-like
        Latitude in radians.
    s : float
        Stretch factor. Values greater than 1 compress the grid near
        the south pole and expand it near the north pole.

    Returns
    -------
    x : array-like
        Longitude in radians (unchanged by this transform).
    y : numpy.ndarray
        Stretched latitude in radians.

    Notes
    -----
    The Schmidt parameter D is computed as (1 - s**2) / (1 + s**2).
    """
    D = (1 - s ** 2) / (1 + s ** 2)
    y = np.arcsin((D + np.sin(y)) / (1 + D * np.sin(y)))
    return x, y


def scs_transform(x, y, s, tx, ty):
    """
    Applies the Stretched Cubed-Sphere (SCS) transform to a
    longitude/latitude grid.  This combines a Schmidt stretch with
    rotations about the x- and z-axes to relocate the high-resolution
    region to a target longitude and latitude.

    Parameters
    ----------
    x : array-like
        Longitude in degrees.
    y : array-like
        Latitude in degrees.
    s : float
        Stretch factor. Values greater than 1 increase the resolution
        near the target point.
    tx : float
        Target longitude in degrees. The center of the high-resolution
        region will be placed at this longitude.
    ty : float
        Target latitude in degrees. The center of the high-resolution
        region will be placed at this latitude.

    Returns
    -------
    x : numpy.ndarray
        Transformed longitude in degrees.
    y : numpy.ndarray
        Transformed latitude in degrees.

    Notes
    -----
    The transform proceeds in the following steps:

    1. Convert all angles from degrees to radians.
    2. Apply the Schmidt transform to stretch the grid.
    3. Convert to Cartesian coordinates.
    4. Rotate about the y-axis by (ty - y0) to shift the target latitude.
    5. Rotate about the z-axis by (tx - x0) to shift the target longitude.
    6. Convert back to spherical coordinates and then to degrees.
    """
    # Convert xy to radians
    x = x * np.pi / 180
    y = y * np.pi / 180
    tx = tx * np.pi / 180
    ty = ty * np.pi / 180
    # Calculate rotation about x, and z axes
    x0 = np.pi
    y0 = -np.pi / 2
    theta_x = ty - y0
    theta_z = tx - x0
    # Apply schmidt transform
    x, y = schmidt_transform(x, y, s)
    # Convert to cartesian coordinates
    x, y, z = spherical_to_cartesian(x, y)
    # Rotate about x axis
    xaxis = np.array([0, 1, 0])
    x, y, z = rotate_vectors(x, y, z, xaxis, theta_x)
    # Rotate about z axis
    zaxis = np.array([0, 0, 1])
    x, y, z = rotate_vectors(x, y, z, zaxis, theta_z)
    # Convert back to spherical coordinates
    x, y = cartesian_to_spherical(x, y, z)
    # Convert back to degrees and return
    x = x * 180 / np.pi
    y = y * 180 / np.pi
    return x, y
