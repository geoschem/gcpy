"""
Internal utilities for managing datetime objects and strings
"""
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np


def get_timestamp_string(date_array):
    """
    Convenience function returning the datetime timestamp based on the given input

    Args:
        date_array: array
            Array of integers corresponding to [year, month, day, hour, minute, second].
            Any integers not provided will be padded accordingly
    Returns:
        date_str: string
            string in datetime format (eg. 2019-01-01T00:00:00Z)
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

    Args:
        start_date: numpy.datetime64
            numpy datetime64 object
        n_months: integer
    Returns:
        new_date: numpy.datetime64
            numpy datetime64 object with exactly n_months added to the date
    """
    new_date = start_date.astype(datetime) + relativedelta(months=n_months)
    return np.datetime64(new_date)


def is_full_year(start_date, end_date):
    """
    Verifies if two dates are a full year starting Jan 1.

    Args:
        start_date: numpy.datetime64
            numpy datetime64 object
        end_date: numpy.datetime64
            numpy datetime64 object
    Returns: boolean
    """
    return (
        add_months(start_date, 12) == end_date
        and start_date.astype(datetime).month == 1
        and start_date.astype(datetime).day == 1
    )
