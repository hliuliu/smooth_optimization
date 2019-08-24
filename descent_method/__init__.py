

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

class DescentMethod(object):

	def __init__(self,x0,f,g,H=None, B0=None, tol = 1e-7, max_iter= 10**8):
		'''
			x0: starting vector
			f: smooth objective function
			g: gradint of f
			H: Hessian  of f
			B0: Approximation to H(x0)
		'''
		self.xk = x0
		self.n = len(x0)
		self.num_iter = 0
		self.objective = f
		self.gradient = g
		self.hessian = H
		self.Bk = B0
		self.direction = None
		self.step_length = 0
		self.tol = tol # tolerence

		self.max_iter = max_iter

		# set default functions

		## search direction: by default, pick a random vector with components in (0,1). Naive
		self.search_dir_fn = lambda f,g,xk: np.random.random(self.n)

		## Step length: always 1 by default
		self.step_length_fn = lambda f,g,xk,dk: 1

		## stopping crietrion: stop when g(xk) has length less then the tolerance
		self.stopping_fn = lambda g,xk: norm(g(*tvec(xk) ))<self.tol

		## Hessian approximation update: unchanged by default
		self.hessian_approx_update = lambda Bk,xk,dk,tk: Bk


	def iterate(self):

		if self.converged():
			return False # stopping criteron met. Unwilling to iterate further

		# pick a descent direction
		self.pick_direction()

		# pick step length
		self.pick_step_length()

		# update iterate
		self.update_iterate()

		# update Bk
		self.update_hessian_approximation()

		self.num_iter += 1

		return True

	def pick_direction(self):
		self.direction = self.search_dir_fn(self.objective,self.gradient,self.xk)

	def pick_step_length(self):
		self.step_length = self.step_length_fn(self.objective,self.gradient,self.xk,self.direction)

	def update_iterate(self):
		self.xk += self.step_length*self.direction

	def update_hessian_approximation(self):
		self.Bk = self.hessian_approx_update(self.Bk,self.xk,self.direction,self.step_length)

	def save_current_iterate(self,attrs, hists):
		'''
			attrs: list[str], names of arrtibutes of self
			hists: list[list[*]], storage space for attribute values
			The lists can also be numpy arrays.

			Example:
				xs, sds = [],[]

				Calling 
					self.save_current_iterate(['xk','direction'],[xs,sds])
				is equivalent to the code below

					xs.append(self.xk)
					sds.append(self.direction)
		'''

		for atname, store in zip(attrs,hists):
			store.append(getattr(self, atname))

	def optimize(self, save_iters= False, attrs=None,hists=None):
		'''
			Run the (reminder of) the algorithm to approximate a local minimizer
		'''

		while self.num_iter<self.max_iter and self.iterate():
			if save_iters:
				self.save_current_iterate(attrs,hists)


	def converged(self):
		return stopping_fn(self.gradient,self.xk)


	def report(self):
		if self.converged():
			print 'Iterative method successful:', 'Terminated after %d iterations'%self.num_iter
			print 'Optimal Point %s of value = %.5f'%(self.xk,self.objective(*tvec(self.xk) ))
			print 'Gradient error: %.5f'%self.gradient(*tvec(self.xk) )

		elif self.num_iter>=self.max_iter:
			print 'Reached a maximum of %d iterations'%self.max_iter
			print 'Method may not be convergent'
			print 'Gradient error: %.5f'%self.gradient(*tvec(self.xk) )

		else:
			print 'The optimizer has not yet finished.'
			print 'Current iterate', self.xk
			print 'Current value', self.objective(*tvec(self.xk) )









