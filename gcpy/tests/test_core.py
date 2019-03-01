'''Unit tests for methods in core.py.'''

import pytest
from gcpy.core import *

def test_filter_names():
    '''Unit test for the filter_names routine.'''
    names = ['EmisCO', 'EmisO3', 'AREA']
    assert filter_names(names, 'Emis') == ['EmisCO', 'EmisO3']
