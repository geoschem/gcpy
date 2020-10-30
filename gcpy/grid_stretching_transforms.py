import numpy as np


def rotate_vectors(x, y, z, k, theta):
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    v = np.moveaxis(np.array([x, y, z]), 0, -1)  # shape: (..., 3)
    v = v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + \
        k[np.newaxis, :] * np.dot(v, k)[:, np.newaxis] * (1 - np.cos(theta))
    return v[..., 0], v[..., 1], v[..., 2]


def cartesian_to_spherical(x, y, z):
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    # Calculate x,y in spherical coordinates
    y_sph = np.arcsin(z)
    x_sph = np.arctan2(y, x)
    return x_sph, y_sph


def spherical_to_cartesian(x, y):
    x_car = np.cos(y) * np.cos(x)
    y_car = np.cos(y) * np.sin(x)
    z_car = np.sin(y)
    return x_car, y_car, z_car


def schmidt_transform(x, y, s):
    D = (1 - s ** 2) / (1 + s ** 2)
    y = np.arcsin((D + np.sin(y)) / (1 + D * np.sin(y)))
    return x, y


def scs_transform(x, y, s, tx, ty):
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
