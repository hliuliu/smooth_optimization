
import numpy as np


class ARCModel(object):
	'''
		Model function used in ARC algorithm
	'''

	def __init__(self, f, x0, g, H=None, B0=None, sigma0=1):
		'''
			f: objective function to minimize
			x0: starting vector
			g: gradient of f
			H: exact hessian of f (set to None if unknown)
			B0: initial approximation to the hessian of f at iterate x0
			sigma0>0: initial regularization parameter
		'''
		self.f = f
		self.xk = x0
		self.g  =g
		self.H = H
		if B0 is None:
			if H is None:
				raise ValueError('Unknown Hessian: initial guess must be specified')
			B0 = H(*x0)
		self.Bk = B0
		self.sigma = sigma0
		self.sigma_update_fn = lambda x: x
		self.hessian_update_function = lambda B: B if H is None else H(*(self.xk))







