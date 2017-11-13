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


def strptime(date_string, format):
    """ Parse a datetime.datetime object from a string with the given format.

    Parameters
    ----------
    date_string : str
        Any given timestamp represented as a string
    format : str
        Format to use when decoding `date_string`

    Returns
    -------
    datetime.datetime object corresponding to decoded string.

    """
    return datetime.datetime.strptime(date_string, format)


def tau_to_yymmdd(tau, nymd0="19850101", nhms0="000000"):
    """ Converts a tau value (elapsed hours between the current date/time and
    the beginning of an epoch) into a calendar date and time value.

    Parameters
    ----------
    tau : int
        The number of elapsed hours
    nymd0 : str

    nhms0

    Returns
    -------

    """

    delta_t = pd.Timedelta(str(int(tau)*60)+'M')
    dt0 = pd.Timestamp(str(nymd0) + " " + str(nhms0))
    result = delta_t + dt0
    result = result.to_pydatetime()

    return {'year': result.year, 'month': result.month, 'day': result.day,
            'hour': result.hour, 'minute': result.minute,
            'second': result.second}


def ymd_to_date(year, month, day):
    """ Given year, month, day (or hour, minute, second) values, returns a
    variable in YYYYMMDD (or HHMMSS) format.

    Parameters
    ----------
    year, month, day: int
        Bits corresponding to year (or hour), month (or minute), and day
        (or second) values

    Returns
    -------
    str containing the date components

    """

    return str(year) + str(month) + str(day)