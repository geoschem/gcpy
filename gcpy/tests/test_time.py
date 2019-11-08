'''Unit tests for methods in time.py'''

from __future__ import print_function

from gcpy.time import *


def test_date_to_ymd():

    # YMD -> dict
    assert date_to_ymd("20170113") == {'year': 2017, 'month': 1, 'day': 13}

    # HMS -> dict
    assert date_to_ymd("134505") == {'hour': 13, 'minute': 45, 'second': 5}


def test_tau_to_ymd():
    
    # 1970/01/01 00:00 GMT
    assert tau_to_yymmdd(-131496.0) == {'year': 1970, 'month': 1, 
                                        'day': 1, 'hour': 0,    
                                        'minute': 0, 'second': 0 }
    # 1985/01/01 00:00 GMT
    assert tau_to_yymmdd(0.0) == {'year': 1985, 'month': 1, 
                                  'day': 1, 'hour': 0,    
                                  'minute': 0, 'second': 0 }

    # 2020/02/29 15:00 GMT
    assert tau_to_yymmdd(308223.0) ==  {'year': 2020, 'month': 2, 
                                        'day': 29, 'hour': 15,    
                                        'minute': 0, 'second': 0 }


