

import sys,os

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)

import numpy as np
from numpy.linalg import norm

from mat_helper import to_vector as tvec


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
		self.n = len(x0)
		if B0 is None:
			if H is None:
				raise ValueError('Unknown Hessian: initial guess must be specified')
			B0 = H(*tvec(x0))
		self.Bk = B0
		self.sigma = sigma0
		self.sigma_update_fn = lambda x: x
		self.hessian_update_function = lambda B,x,d: B if H is None else H(*tvec(x))

	def get_hessian_approx(self):
		return self.Bk if self.Bk is not None else H(*tvec(self.xk))
	
	def __call__(self, s):
		'''
			s: direction vector
		'''
		xk = self.xk
		Bk = self.get_hessian_approx()
		fk = self.f(*tvec(xk))
		gk = self.g(*tvec(xk))
		return fk+s.T*gk+self.sigma*norm(s)**3/3.



class ARCSubspaceModel(ARCModel):

	def __init__(self, f, x0, g, H=None, B0=None, sigma0=1, Q=None, T=None):
		'''
			Q is a rectangular matrix such that Q^T*Q = I, n-dimensional
			T = Q^T*B*Q, which will be computed if not specified,
				where B is whichever matrix we use for Hessian approximation
				(could be H(xk) itself)
		'''
		ARCModel.__init__(self, f, x0, g, H, B0, sigma0)
		# if Q is None:
		# 	Q = np.eye(self.n)
			#T = ARCModel.get_hessian_approx(self)
		self.T = T
		self.Q = Q
		#self._actual_g = self.g
		self.g = (
			lambda g: (lambda x: self.Q.T*g(*tvec(x)) ) 
			) (self.g)

	def get_hessian_approx(self):
		if self.T:
			return T
		Bk = ARCModel.get_hessian_approx(self)
		if self.Q is not None:
			return Bk

		return self.Q.T*Bk*self.Q








