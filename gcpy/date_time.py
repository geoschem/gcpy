"""
Internal utilities for managing datetime objects and strings
"""
from datetime import datetime


def get_timestamp_string(
        date_array
):
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

    date_str = str(datetime(*date_array)).replace(' ', 'T') + 'Z'
    
    return date_str
