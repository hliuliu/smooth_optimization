

import numpy as np
import matplotlib.pyplot as plt



def product(seq):
	ans = 1
	for s in seq:
		ans *= s
	return ans


def zero_pad(seq, length):
	return list(seq)+[0]*(max(0, length-len(seq)))


def init_poly_vars(variables):
	return [Polynomial([v],[(1,[1])]) for v in variables]

def constant_poly(c):
	return Polynomial([],[(c,[])])


class Polynomial(object):

	def _check_exponents(self, seq):
		for i in seq:
			if type(i) not in [int,long]:
				raise TypeError('exponent not integral: '+str(i))
			if i<0:
				raise TypeError('exponent is negative: '+str(i))



	def __init__(self, variables, terms):
		'''
			terms = [ (coef, [exponents] ) OR (coef, { var: exponent, ... }), ... ]
			variables = ['x','y','z','x0','y1', ... ]
		'''

		self.variables = list(variables)
		self.variables_rlookup = {v:i for i,v in enumerate(variables)}
		self.terms = {}

		for coef, expcoll in terms:
			if coef==0:
				continue
			if type(expcoll) is list:
				expcoll = tuple(zero_pad(expcoll,len(variables)))
			else:
				seq = [0]*len(variables)
				for v,e in expcoll:
					seq[self.variables_rlookup[v]] = e
				expcoll = tuple(seq)
			self._check_exponents(expcoll)
			if coef == int(coef):
				coef = int(coef)
			self.terms[expcoll] = coef

	#@staticmethod
	def _monomial_str(self, exponents, coef):
		if coef == 0:
			return ''
		ans = '' if abs(coef) == 1 else str(abs(coef))
		if coef<0:
			ans = '-'+ans
		for e,v in zip(exponents,self.variables):
			if e:
				ans += '{}^{}'.format(v,e) if e>1 else str(v)

		if not ans:
			return '1'
		return ans



	def __str__(self):
		the_terms = map(self._monomial_str, self.terms.keys(), self.terms.values())
		if not the_terms:
			return '0'
		ans = ''
		for tm in the_terms:
			if ans and tm[0]!='-':
				ans += '+'
			ans += tm
		return ans


	def __call__(self, *args, **kwargs):
		'''
			If kwargs is used, unspecified params default to 0
		'''
		if not args and kwargs:
			args = [0]*len(self.variables)
			for k,v in kwargs.iteritems():
				args[self.variables_rlookup[k]] = v

		if len(args)!=len(self.variables):
			raise ValueError('number of arguments does not match number of variables')

		#assoc = dict(zip(self.variables,args))

		ans = 0

		for expl, coef in self.terms.iteritems():
			ans += coef * product([b**e for b,e in zip(args, expl)])

		return ans

	def varlist(self):
		return list(self.variables)

	def var_index(self, varname):
		return self.variables_rlookup[varname]

	def __neg__(self):
		newterms = []

		for expl, coef in self.terms.iteritems():
			newterms.append((-coef, list(expl)))

		return Polynomial(self.variables, newterms)

	@staticmethod
	def _to_var_dict(variables, exponents):
		vd = {}
		for v,e in zip(variables,exponents):
			if e:
				vd[v]=e
		return vd

	def __eq__(self,other):
		if not isinstance(other, Polynomial):
			other = constant_poly(other)

		svb = [(self._to_var_dict(self.variables, e),c) for e,c in self.terms.iteritems()]
		ovb = [(self._to_var_dict(other.variables, e),c) for e,c in other.terms.iteritems()]
		#print svb,ovb
		if len(svb)!=len(ovb):
			return False
		svb.sort()
		ovb.sort()

		for (vds, cs), (vdo, co) in zip(svb,ovb):
			#print vds,cs
			#print vdo,co
			#print 
			if vds!=vdo or cs!=co:
				return False
		return True















