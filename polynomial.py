

import numpy as np
# import matplotlib.pyplot as plt

from numbers import Number as _num



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

def split_list(seq, used):
	incl, excl = [],[]
	for i,j in enumerate(seq):
		(incl if i in used else excl).append(j)
	return incl,excl


def sorted_variables(poly, key=None):
	if key is None:
		key = lambda x: x
	if not isinstance(poly, Polynomial):
		raise TypeError('First argument must be a Polynomial instance')
	newvars = sorted(poly.variables, key=key)
	newterms = []
	for expl,coef in poly.terms.iteritems():
		expl = poly._to_var_dict(poly.variables,expl)
		newterms.append((coef, expl))

	return Polynomial(newvars,newterms)

def to_poly(poly):
	if not isinstance(poly,Polynomial):
		if isinstance(poly,_num):
			return constant_poly(poly)
		raise TypeError("Cannot convert object '{}' of type '{}' to polynomial".format(poly,type(poly)))
	return poly


class Polynomial(object):

	AUTO_SORT_VARIABLES = False


	@staticmethod
	def _check_exponents(seq):
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
		if Polynomial.AUTO_SORT_VARIABLES:
			self.variables.sort()
		self.variables_rlookup = {v:i for i,v in enumerate(self.variables)}
		self.terms = {}
		self.deg = None

		for coef, expcoll in terms:
			if coef==0:
				continue
			if type(expcoll) is list:
				expcoll = tuple(zero_pad(expcoll,len(variables)))
				expcoll = dict(zip(variables, expcoll))
				# print expcoll
			
			seq = [0]*len(variables)
			for v,e in expcoll.iteritems():
				seq[self.variables_rlookup[v]] = e
			expcoll = tuple(seq)

			self._check_exponents(expcoll)
			if coef == int(coef):
				coef = int(coef)
			self.terms[expcoll] = self.terms.get(expcoll,0)+coef
			if self.terms[expcoll] == int(self.terms[expcoll]):
				self.terms[expcoll] = int(self.terms[expcoll])
			if not self.terms[expcoll]:
				del self.terms[expcoll]


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

		if ans in ['','-']:
			ans += '1'
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

	def __ne__(self,other):
		return not self==other

	def is_constant(self):
		for expl in self.terms:
			if sum(expl):
				return False
		return True

	def __float__(self):
		if not self.is_constant():
			raise ValueError('Polynomial is nonconstant')
		if not self.terms:
			return 0.0
		[coef] = self.terms.values()
		return float(coef)

	def __int__(self):
		fl = float(self)
		if fl!=int(fl):
			raise ValueError('Constant polynomial is not an integer')
		return int(fl)

	def coef_of(self, monomial, partial = False):
		'''
			monomial = {x0:e0,x1:e1,...,xn:en}
			xi are the variables
			Returns the coefficient of x0^e0*x1^e1*...*xn^en

			If partial is set to True, then a Polynomial is returned.
			Otherwise, a number is returned.
				In this case, if m>n, and set(self.varlist()) == {x0,x1,...,xm},
					the numer returned is the coefficient of x0^e0*x1^e1*...*xn^en* x{n+1}^0* ... * xm^0

			If a variable, y not in self.varlist() appears with exponent e, then the following rule is applied:
				1. If e=0, this variable is simply ignored, hence implicitly discarded.
				2. If e>0, then 0 is returned (as a polynomial if partial is True, or a raw number otherwise)

		'''
		explist = [0]*len(self.variables)
		used = set()
		for x,e in monomial.iteritems():
			if x not in self.variables_rlookup:
				if e:
					return (constant_poly if partial else int)(0)
				continue
			explist[self.variables_rlookup[x]] = e
			used.add(self.variables_rlookup[x])

		if not partial:
			return self.terms.get(tuple(explist),0)

		varincl,varexcl = split_list(self.variables,used)
		expincl,_ = split_list(explist, used)
		newterms = []

		for texp, coef in self.terms.iteritems():
			incl,excl = split_list(texp,used)
			if incl == expincl:
				newterms.append((coef, excl))

		return Polynomial(varexcl,newterms)

	@staticmethod
	def get_sum(*polys):
		'''
			polys = [p1,p2,...,pn]
			each pi is a Polynomial
			returns s = p1+p2+...+pn, in the usual sense of Polynomial addition
			s.varlist() is the union of the pi.varlist(), in an arbitrary order
		'''
		newvars = set()
		for p in polys:
			newvars |= set(p.varlist())
		newvars = list(newvars)
		newterms = []
		for p in polys:
			for expl,coef in p.terms.iteritems():
				newterms.append((coef, Polynomial._to_var_dict(p.varlist(),expl)) )

		return Polynomial(newvars,newterms)


	def __add__(self, other):
		if not isinstance(other,Polynomial):
			other = constant_poly(other)

		return self.get_sum(self,other)

	def __radd__(self,other):
		return self+other

	def __sub__(self,other):
		return self+(-other)

	def __rsub__(self,other):
		return -(self-other)

	def is_zero(self):
		return not self.terms

	def __nonzero__(self):
		return bool(self.terms)

	def degree(self):
		if self.deg is None:
			self.deg = max([sum(ex) for ex in self.terms])
		return self.deg

	@staticmethod
	def _monomial_product(mon1,mon2):
		prod = {}

		for m in [mon1,mon2]:
			for v,e in m.iteritems():
				prod[v] = prod.get(v,0)+e

		return prod

	def __mul__(self,other):
		from poly_array import PolynomialArray as Parr
		if isinstance(other,Parr):
			return other*self
		if not isinstance(other,Polynomial):
			other = constant_poly(other)
		newvars = list(set(self.variables) | set(other.variables))
		newterms = []
		svb = [(self._to_var_dict(self.variables, e),c) for e,c in self.terms.iteritems()]
		ovb = [(self._to_var_dict(other.variables, e),c) for e,c in other.terms.iteritems()]

		for se, sc in svb:
			for oe,oc in ovb:
				newterms.append((sc*oc, self._monomial_product(se,oe)))

		return Polynomial(newvars,newterms)

	def __rmul__(self,other):
		return self*other

	def derivative(self,wrt):
		if wrt not in self.variables_rlookup:
			return constant_poly(0)
		newterms = []

		for expl,coef in self.terms.iteritems():
			if expl[self.variables_rlookup[wrt]]:
				expl = self._to_var_dict(self.variables, expl)
				coef *= expl[wrt]
				expl[wrt]-=1
				newterms.append((coef, expl))

		return Polynomial(self.variables, newterms)

	@staticmethod
	def _recpow(poly, n):
		if n==0:
			return constant_poly(1)
		if n==1:
			return poly
		psq = poly*poly
		n,r = divmod(n,2)
		ans = Polynomial._recpow(psq,n)
		if r:
			ans *= poly
		return ans

	def __pow__(self, n):
		self._check_exponents([n])
		return self._recpow(self,n)

	def gradient(self):
		from poly_array import PolynomialVector as pv
		gd = []

		for v in self.variables:
			gd.append(self.derivative(v))

		return pv(gd)





















