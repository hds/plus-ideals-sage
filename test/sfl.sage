# -*- coding: utf-8 -*-
import sys
import os

testpath = os.path.dirname(os.path.realpath(__file__))
libpath = os.path.realpath(os.path.join(testpath, '..', 'lib'))
sys.path += [ libpath ]

from plusideals import montes, path_of_precision, change_precision


print path_of_precision(121, 6)

p = 3
zp = Zp(p, type='capped-abs', prec=5)
a = zp(67)
zpx.<x> = PolynomialRing(zp)

f = zpx([67, 7, 0, 1])

print f
g = change_precision(f, 2)
print g

print change_precision(g, 4)
