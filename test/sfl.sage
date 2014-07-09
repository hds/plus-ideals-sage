# -*- coding: utf-8 -*-
import sys
import os

testpath = os.path.dirname(os.path.realpath(__file__))
libpath = os.path.realpath(os.path.join(testpath, '..', 'lib'))
sys.path += [ libpath ]

from plusideals import montes, path_of_precision, change_precision


p = 2
f = x^12 + p^2*x^6 + p^4*x^3 + p^6

K = NumberField(f, 'theta')

om_reps = montes(K, p)

tt = om_reps[0]

print "\n--------------------------------"
print "Before: {0} (slope: {1})"\
        .format(tt.rth_level().phi, tt.rth_level().slope)

tt.single_factor_lifting(20)

print "After: {0} (slope: {1})"\
        .format(tt.rth_level().phi, tt.rth_level().slope)

