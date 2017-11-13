
from __future__ import print_function

import pytest

import numpy as np
import pandas as pd

from xarray.testing import assert_equal

from gcpy.time import *


def test_date_to_ymd():

    # YMD -> dict
    assert date_to_ymd("20170113") == {'year': 2017, 'month': 1, 'day': 13}

    # HMS -> dict
    assert date_to_ymd("134505") == {'hour': 13, 'minute': 45, 'second': 5}