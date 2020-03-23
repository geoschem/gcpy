#!/usr/bin/env python

'''Unit tests for methods in core.py.'''

from gcpy.core import filter_names

def test_filter_names():
    '''Unit test for the filter_names routine.'''
    names = ['EmisCO', 'EmisO3', 'AREA']
    assert filter_names(names, 'Emis') == ['EmisCO', 'EmisO3']
