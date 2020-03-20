#!/usr/bin/env python

'''Unit tests for methods in units.py'''

from gcpy.units import adjust_units

def test_adjust_units():

    for v in ['kg/m2/s', 'kgm-2s-1', 'kgm^-2s^-1']:
        assert adjust_units(v) == 'kg/m2/s'

    for v in ['kgC/m2/s', 'kgCm-2s-1', 'kgCm^-2s^-1',
              'kgc/m2/s', 'kgcm-2s-1', 'kgcm^-2s^-1']:
        assert adjust_units(v) == 'kgC/m2/s'

    for v in ['molec/cm2/s', 'moleccm-2s-1', 'moleccm^-2s^-1']:
        assert adjust_units(v) == 'molec/cm2/s'
        
    for v in ['atoms C/cm2/s', 'atomsC/cm2/s']:
        assert adjust_units(v) == 'atomsC/cm2/s'
