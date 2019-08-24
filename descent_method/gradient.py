

import sys,os

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)


import numpy as np

from descent_method import DescentMethod
from mat_helper import to_vector as tvec


class GradientMethod(DescentMethod):
	'''
		Gradient (aka Steepest Descent) Method
	'''

	def __init__(self,x0,f,g,H=None, B0=None, tol = 1e-7, max_iter= 10**8):

		DescentMethod.__init__(self,x0,f,g, H, B0, tol, max_iter)
		self.search_dir_fn = lambda f,g,xk: -g(*tvec(xk))


class ScaledGradientMethod(GradientMethod):

	def __init__(self,x0,f,g,H=None, B0=None, tol = 1e-7, max_iter= 10**8):
		GradientMethod.__init__(self,x0,f,g, H, B0, tol, max_iter)

		self.scaled_matrix_fn = lambda f,g,H,xk: np.identity(self.n)

		#self.unscaled_step_length_fn = self.step_length_fn 

		self.scaled_matrix = None


	def pick_direction(self):
		self.scaled_matrix = self.scaled_matrix_fn(self.objective, self.gradient, self.hessian, self.xk)
		GradientMethod.pick_direction(self)
		self.direction = self.scaled_matrix* self.direction



def inverse_diag_hessian(f,g,H,xk):
	Hk = H(xk)
	dg = Hk.diagonal()
	return np.matrix(np.diag(1./dg))




