
import numpy as np

from numbers import Number as _num

import sys

class _err(object):
	def __init__(self,tol):
		self.tol = tol

	@property
	def tol(self):
		return self._tol

	@tol.setter
	def tol(self, tol):
		if not isinstance(tol,_num):
			raise TypeError('Cannot set tolerence to a non-numeric value')

		tol = abs(tol)
		self._tol = tol

DEFAULT_TOLERANCE = 1e-7

sys.modules['{}.error'.format(__name__)] = _err(DEFAULT_TOLERANCE)

del DEFAULT_TOLERANCE
del _err






def close(x,xappr):
	'''
		x and xapprr and both numbers OR both numpy arrays
			of the smae shape
	'''
	return np.abs(x-xappr)<error.tol

def round_if_near(x):
	test = near_int(x)
	if type(test) in  (bool,np.bool_):
		return int(np.round(x)) if test else x
	if np.all(test):
		return x.astype(int)
	arr = np.vectorize(lambda test, tfunc,x: tfunc(x) if test else x)
	arr = arr(test, np.round, x)
	return arr


def near_int(x):
	'''
		x: number or numpy array
	'''
	return np.abs(np.round(x)-x)<error.tol

def near_zero(x):

	return np.abs(np.round(x))<error.tol




