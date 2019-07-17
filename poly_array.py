

import polynomial as ply

import numpy as np

#from numbers import Number as _num




Poly = ply.Polynomial

zero_poly = ply.constant_poly(0)



class PolynomialArray(object):

	AUTO_NUMPY = True

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

	def to_numpy_array(self, indices = ()):
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

	def __call__(self,*args,**kwargs):
		'''
			Calls each polynomial entries, and replace that polynomial with the result.
			if all the entries are constant, then,
				if PoynomialArray.AUTO_NUMPY is True then a Numpy array with those numbers and same shape is returned.
				else, a PolynomialArray of constant polynomials is returned.
			else, a PolynomialArray with the returned Polynomials is returned. 
		'''

		if not self.shape:
			val = self.array(*args,**kwargs)
			if PolynomialArray.AUTO_NUMPY:
				if not isinstance(val, Poly) or val.is_constant():
					return int(val) if val == int(val) else float(val)
			if not isinstance(val,Poly):
				val = ply.constant_poly(val)
			return PolynomialArray((), val)



		new_array = [pa(*args,**kwargs) for pa in self.array]

		indices = [i for i in xrange(len(new_array)) if isinstance(new_array[i],PolynomialArray)]

		if not indices:
			return np.array(new_array)

		if len(indices)!=len(new_array):
			# disable numpy conversion
			PolynomialArray.AUTO_NUMPY = False
			for i,_ in enumerate(new_array):
				if not isinstance(new_array[i],PolynomialArray):
					new_array[i] = self.array[i](*args,**kwargs)
			PolynomialArray.AUTO_NUMPY = True

		npa = PolynomialArray(0)
		npa.array = new_array
		npa.shape = self.shape
		return npa

















