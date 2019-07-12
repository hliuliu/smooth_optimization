

import Polynomial as ply

import numpy as np


Poly = ply.Polynomial

zero_poly = ply.constant_poly(0)



class PolynomialArray(object):

	@staticmethod
	def _check_shape(seq):
		for i in seq:
			if type(i) not in [int,long]:
				raise TypeError('dimension not integral: '+str(i))
			if i<0:
				raise TypeError('dimension is negative: '+str(i))

	@staticmethod
	def _rec_array_constuct(shape,default):
		if not shape:
			return default
		rec_shape = shape[1:]
		return [PolynomialArray(rec_shape,default) for _ in [0]*shape[0]]


	def __new__(cls, shape, default = zero_poly):
		'''
			shape: a vector of nonnegative integers [n0, n1, ..., nk], or a single number, n0
			creates a n0 x n1 x ... x nk array of Polynomials, initially populated with default

			Note: for ideal purposes, polynomial entries should have the same list of variables in the same order.
		'''

		if not isinstance(default, Poly):
			default = ply.constant_poly(default)

		try: 
			iter(shape)
		except:
			shape = (shape,)

		cls._check_shape(shape)

		#if not shape:
		#	return  default

		obj = super(PolynomialArray, cls).__new__(cls)
		obj.shape = tuple(shape)


		obj.array = obj._rec_array_constuct(obj.shape,default)
		return obj


	def dimensions(self):
		return self.shape

	def to_numpy_array(self, indicies = ()):
		'''
			Convert to the numpy array with same shape and data.
			Only works if all polynomials are constant.
		'''
		if not self.shape:
			if not self.array.is_constant():
				raise ValueError('Polynomial {} is not constant: at index {}'.format(self.array, indices))
			num = float(self.array)
			if num==int(num):
				num = int(num)
			return num

		A = []
		for i, p in enumerate(self.array):
			A.append(p.to_numpy_array(indices+(i,)))

		return np.array(A)














