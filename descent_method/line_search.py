
import sys,os

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)



import numpy as np
from numpy.random import random

from mat_helper import to_vector as tvec

class LineSearch(object):

	def __init__(self, routine=None):
		'''
			routine: a function (f,g,x,d) -> t where
			f,g: objective function,gradient
			x: Column vector. 
				also the current iterate of the associated optimization routine
			d: The search direction
			t: step length
		'''
		if routine is None:
			routine = lambda f,g,x,d: random()
		self.routine = routine

	def __call__(self, f,g,x,d):
		return self.routine(f,g,x,d)


class ConstantStepsize(LineSearch):

	def __init__(self, t):
		LineSearch.__init__(self, lambda f,g,x,d: t)
		self.step_length = t


class ExactLS(LineSearch):

	def __init__(self, min_routine, ray_fn=None):
		'''
			min_routine: A procedure that can finds and returns the minimum
				argument of t -> f(x+td) in the t parameter with x,d fixed.
				Should take in the function as above as its only argument.
			ray_fn: The functional argument of min_routine (except with extra x,d parameters). Optional
				Will be computed based on f on the fly if not provided.
				This funtion MUST assume f is known, and  (x,d,t) -> f(x+td)
				is impemented explicitly.

				For example, if f(s) = s^2 is known, one can simmply pass ray_fn as 
					(x,d,t) -> (x+td)^2, without referring to f directly
		'''
		self.ray_fn = ray_fn
		self.min_routine = min_routine

		def find_opt_length(f,g,x,d):
			if self.ray_fn is None:
				ray_fn = lambda x,d,t: f(x+td)
			else:
				ray_fn = self.ray_fn
			return self.min_routine(lambda t: ray_fn(x,d,t))

		LineSearch.__init__(self, find_opt_length)


def _backtracking_routine(f,g,x,d,s,a,b):
	fx = f(x)
	gtd = (g(x).T*d)[0,0]
	t=s
	while 1:
		td = t*d
		fxtd = f(x+td)
		
		if fx-fxtd>= -a*t*gtd:
			return t
		t = b*t



class Backtracking(LineSearch):

	def __init__(self, s, a, b):
		'''
			s>0: starting step length
			0<a<1: factor used in sufficient decrease lower bound
			0<b<1: factor to scale the step length by each iteration
		'''
		self.s = s
		self.a = a
		self.b = b

		LineSearch.__init__(
				lambda f,g,x,d: _backtracking_routine(
						f,g,x,d,self.s,self.a,self.b
					)
			)









		



