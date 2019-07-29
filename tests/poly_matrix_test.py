

import os,sys

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)


import poly_array as parr
import polynomial as ply
import numpy as np

poly = ply.Polynomial
pmatrix = parr.PolynomialMatrix
pvec = parr.PolynomialVector

ply.Polynomial.AUTO_SORT_VARIABLES = True

parr.PolynomialArray.AUTO_NUMPY = True

x,y,z = ply.init_poly_vars('xyz')

f1 = x**3+y**2+z**5
_rt =(3*x**2-2*x*y+z**3)
f2 = _rt**2

f3 = (x**3+y)**2+z+2


g1 = f1.gradient()
g2 = f2.gradient()
g3 = f3.gradient()
print g1
assert g1== pvec([3*x**2,2*y,5*z**4])
assert g2== 2*_rt*pvec([(6*x-2*y), -2*x, 3*z**2])

M = pmatrix([g1,g2,g3])

assert type(M)==pmatrix

print M

for i in [0,1,2]:
	assert type(M[i])== pvec
	assert M[i] == eval('g{}'.format(i+1))
	for j in [0,1,2]:
		assert type(M[i][j])==poly

assert M.dimensions() == (3,3)


