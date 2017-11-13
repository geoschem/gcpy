""" Specialized/auxiliary functions for dealing with timestamps and
timeseries of data. """


import datetime
import pandas as pd


def date_to_ymd(date):
    """ Given a date in YYYYMMDD format, return the 'year', 'month', and 'day'
    as separate variables in a dictionary. Also can be used to separate
    a time in HHMMSS format into 'hours', 'minutes', and 'seconds.

    TODO: It's generally bad form to have a single function have multiple
          return types, so see if we can't refactor this more elegantly.

    Parameters
    ----------
    date : str
        A str with the format 'YYYYMMDD' or 'HHmmss'


    Returns
    -------
    dict containing the timestamp components

    Notes
    -----
    N/A

    Examples
    --------

    >>> from gcpy.time import date_to_ymd

    Separate a date into year, month, and day

    >>> ymd = date_to_ymd("20170113")
    >>> ymd['year'], ymd['month'], ymd['day']
    2017 1 13

    Separate a time into hour, minute, and second
    >>> hms = date_to_ymd("134505")
    >>> hms['hour'], hms['minute'], hms['second']
    13 45 5

    """

    typ = 'ymd' if len(date) > 6 else 'hms'
    fmt = '%Y%m%d' if typ == 'ymd' else '%H%M%S'
    dt = datetime.datetime.strptime(date, fmt)

    if typ == 'ymd':
        return {'year': dt.year, 'month': dt.month, 'day': dt.day}
    else:
        return {'hour': dt.hour, 'minute': dt.minute, 'second': dt.second}



def isleap():
    # TODO
    pass


def strdate():
    # TODO
    pass


def tau_to_yymmdd():
    # TODO
    pass


def ymd_to_date():
    # TODO
    pass