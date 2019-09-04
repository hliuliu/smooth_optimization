


import sys,os

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)


from polynomial import Polynomial as Poly
from fractions import Fraction as frac
from math import sqrt

from numerics import near_int, round_if_near
import numerics.error as error
from numerics.error import tol



def linear_root(a,b,tol=tol,use_rational=False):
	'''
		returns the root of ax+b
	'''
	if abs(a)<tol:
		return None if abs(b)>=tol else 'all real numbers'

	tmp = error.tol
	error.tol = tol

	a,b = map(round_if_near, [a,b])

	error.tol = tmp

	if use_rational and set(map(type,[a,b])).issubset({int,frac}):
		return -frac(b)/frac(a)

	return float(-b)/a


def plus_or_minus(a, b):
	return np.array([a+b,a-b])


def quadratic_roots(a,b,c,tol = tol):
	'''
		root of a*x^2+b*x+c
	'''
	tmp = error.tol
	
	if abs(a)<tol:
		return linear_root(b,c,tol)

	error.tol = tol
	a,b,c = map(round_if_near,[a,b,c])

	disc = b**2-4*a*c
	disc = round_if_near(disc)

	ans = -b

	if disc<0:
		ans = plus_or_minus(ans, sqrt(-disc)*1j)
	else:
		ans = plus_or_minus(ans, sqrt(disc))
	ans /= 2.
	ans /= a

	ans = round_if_near(ans)

	error.tol = tmp
	return ans









	 






