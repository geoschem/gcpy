"""
Internal utilities for managing datetime objects and strings
"""
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
from gcpy.util import verify_variable_type


def get_timestamp_string(date_array):
    """
    Convenience function returning the datetime timestamp based on
    the given input.

    Parameters
    ----------
    date_array : array of int
        Array of integers corresponding to [year, month, day, hour,
        minute, second]. Any integers not provided will be padded
        accordingly.

    Returns
    -------
    date_str : str
        String in datetime format (e.g. '2019-01-01T00:00:00Z').
    """
    # converts single integer to array for cases when only year is given
    date_array = [date_array] if isinstance(date_array, int) else date_array

    # datetime function must receive at least three arguments
    while len(date_array) < 3:
        date_array.append(None)

    # set default values for month and day if not present
    date_array[1] = date_array[1] or 1
    date_array[2] = date_array[2] or 1

    date_str = str(datetime(*date_array)).replace(" ", "T") + "Z"

    return date_str


def add_months(start_date, n_months):
    """
    Adds a given number of months to a numpy.datetime64 date.

    Parameters
    ----------
    start_date : numpy.datetime64
        The starting date.
    n_months : int
        Number of months to add to start_date.

    Returns
    -------
    new_date : numpy.datetime64
        Date with exactly n_months added to start_date.
    """
    new_date = start_date.astype(datetime) + relativedelta(months=n_months)
    return np.datetime64(new_date)


def is_full_year(start_date, end_date):
    """
    Verifies if two dates span a full year starting on January 1st.

    Parameters
    ----------
    start_date : numpy.datetime64
        The starting date.
    end_date : numpy.datetime64
        The ending date.

    Returns
    -------
    result : bool
        True if end_date is exactly 12 months after start_date and
        start_date falls on January 1st; False otherwise.
    """
    return (
        add_months(start_date, 12) == end_date
        and start_date.astype(datetime).month == 1
        and start_date.astype(datetime).day == 1
    )


def datetime64_to_str(timestamp, format_str="%Y-%m-%d"):
    """
    Convenience routine to convert a numpy.datetime64 object
    to a date/time string.

    Parameters
    ----------
    timestamp : numpy.datetime64
        Date and time to be converted.
    format_str : str, optional
        Format string for the output date/time string.
        Default value: '%Y-%m-%d'

    Returns
    -------
    date_str : str
        String representation of timestamp formatted according to
        format_str.
    """
    verify_variable_type(timestamp, np.datetime64)

    return timestamp.astype(
        'datetime64[s]').astype(datetime).strftime(format_str)
