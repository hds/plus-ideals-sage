# -*- coding: utf-8 -*-
import sys
sys.path += ['../lib']

from plusideals import montes


Zx.<x> = PolynomialRing(ZZ)
# Simple example
p = 3
f = x^2 + 18*x + 20

# Article example
p = 2
f = x^12 + p^2*x^6 + p^4*x^3 + p^6

K = NumberField(f, 'K1')

montes(K, p)


