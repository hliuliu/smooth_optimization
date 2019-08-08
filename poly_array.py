

import polynomial as ply

import numpy as np

from math import ceil

#from numbers import Number as _num




Poly = ply.Polynomial

zero_poly = ply.constant_poly(0)


def _slice_len(n, sl):
	start,stop,step = sl.indices(n)
	if (stop-start)*step<=0:
		return 0
	return 1+(stop-start)//step


def _itermap(func, *args):
	args = map(iter,args)
	while 1:
		coll = []
		for it in args:
			try:
				nxt = next(it)
			except StopIteration:
				break
			else:
				coll.append(nxt)
		else:
			yield func(*coll)
			continue
		break

def _ident(*args):
	return args

def _iterzip(*args):
	return _itermap(_ident,*args)



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


	@staticmethod
	def from_numpy(arr):
		try:
			iter(arr)
		except:
			return PolynomialArray((),arr)
		polyarr = PolynomialArray(len(arr))

		for i in xrange(len(arr)):
			polyarr.array[i] = polyarr.from_numpy(arr[i])

		polyarr.shape = arr.shape
		return polyarr




	def __new__(cls, shape, default = zero_poly):
		'''
			shape: a vector of nonnegative integers [n0, n1, ..., nk], or a single number, n0
			creates a n0 x n1 x ... x nk array of Polynomials, initially populated with default

			Note: for ideal purposes, polynomial entries should have the same list of variables in the same order.
		'''

		default = ply.to_poly(default)

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


	def __len__(self):
		if not self.shape:
			raise ValueError('Cannot get length of a degenerate singleton array')
		return self.shape[0]


	def __iter__(self):
		if not self.shape:
			raise ValueError('No dimensions. not iterable')
		flag = len(self.shape)==1
		if flag:
			for item in self.array:
				yield item.array
		else:
			for item in self.array:
				yield item

	def to_nested_list(self,poly_to_str=False):
		if not self.shape:
			return self.array if not poly_to_str else str(self.array)
		return [entry.to_nested_list(poly_to_str) for entry in self.array]


	def __str__(self):
		return str(self.to_nested_list(True))

	@staticmethod
	def _get_yield_entry(index, entry, incl):
		return (index,entry) if incl else entry

	def entries(self, incl_indices=False):
		index = 0
		if not self.shape:
			yield self._get_yield_entry(index,self.array, incl_indices)
		else:
			for subarray in self.array:
				for entry in subarray.entries():
					yield self._get_yield_entry(index,entry,incl_indices)
					index +=1

	@staticmethod
	def populate(pa,iter_entries,iter_conv=True):
		'''
			Note:
				Default parameter, iter_conv, is for implemetation details,
					Users should ignore it 
			Populates the PolynomialArray instance with entries for an iterable object.
			Lex ordering is used.
			Population happens until the new entries of the array positions run out, whichever is first
			Returns True iff the iterator is not all used up BEFORE ALL the array positions are changed
			Assumes that all entries are polynomial objects
		'''
		if iter_conv:
			iter_entries = iter(iter_entries)
		
		if not pa.shape:
			nxt = next(iter_entries,None)
			if nxt is None:
				return False
			pa.array = nxt
			return True

		for spa in pa.array:
			if not pa.populate(spa,iter_entries,False):
				return False
		return True

	def derivative(self,wrt):
		der = self.transform(lambda f: f.derivative(wrt))
		if not self.shape:
			return der.array
		return der


	def reshape(self, shape, allow_size_diff= False):
		'''
			Returns an array of the new provided shape.
			Entries are populated in lex order of indices.
			Remaining entries of new array, if any, are 0 by default.
			If the new array is too small, some entries (lex order) of self will not appear in the new array
			If allow_size_diff is et to False, then equalent dimensions is forced. (via an exception otherwise)
		'''
		try:
			iter(shape)
		except:
			shape = (shape,)
		self._check_shape(shape)

		if not allow_size_diff and ply.product(shape)!=ply.product(self.shape):
			raise ValueError('Cannot reshape: number of Polynomial entries differ')

		newpa = PolynomialArray(shape)

		self.populate(newpa, self.entries())

		return newpa

	def __getitem__(self,inds):
		'''
			inds can be of type int,slice, of a tuple consisting of those types.
		'''
		if type(inds) is not tuple:
			inds = (inds,)

		pa = self._get_item_keep_as_array(inds)

		if not pa.shape:
			return pa.array
		return pa


	def __setitem__(self,inds,poly):
		'''
			self[inds] = poly
			currently: only supprots updating a single entry
		'''
		if not isinstance(poly,Poly):
			poly = ply.constant_poly(poly)
		if type(inds) is not tuple:
			inds = (inds,)
		if not self.shape:
			if not inds:
				self.array = poly
				return
		if len(inds)!=len(self.shape):
			raise IndexError('Indices Tuple does not match the dimensions in length')
		pa = self._get_item_keep_as_array(inds[:-1])
		pa.array[inds[-1]] = PolynomialArray((), poly)




	def _get_item_keep_as_array(self,inds):
		
		if inds and not self.shape:
			raise IndexError('Dimension overflow: %d'%inds[0])

		if not inds:
			return self

		currind,inds = inds[0],inds[1:]

		if type(currind) is int:
			return self.array[currind]._get_item_keep_as_array(inds)

		# currind is a slice
		subarray = self.array[currind]
		spa = ([entry._get_item_keep_as_array(inds) for entry in subarray])

		if spa:
			sspa = PolynomialArray(0)
			sspa.shape = (len(spa),)+(spa[0].shape)
			sspa.array = spa
			return sspa

		# slice currind yield an empty sequence
		# still need to retain the shape

		newshape = [0]

		for i,s in zip(inds,self.shape[1:]):
			if type(i) is slice:
				newshape.append(_slice_len(s,i))
			else:
				if -s<=i<s:
					newshape.append(1)
				else:
					raise IndexError('Index out of bounds, {}, of size {}'.format(i,s))

		return PolynomialArray(tuple(newshape))


	def __eq__(self,other):
		'''
			Two PolynomialArray instances are equal if their shape, and corresponding components, are equal
		'''
		return self.equals_ignore_shape(other) and self.shape==other.shape

	def equals_ignore_shape(self,other):
		'''
			returns True iff
				self can be reshaped, resp. lex order, to look like other
			iff
				The 'flattened' versions for self and other and identical

			return value logically equivalent to 
				self.reshape(other.dimensions(),True) == other
		'''
		if not isinstance(other, PolynomialArray):
			return False
		if ply.product(self.shape)!=ply.product(other.shape):
			return False

		for e1,e2 in _iterzip(self.entries(),other.entries()):
			if e1!=e2:
				return False

		return True

	def __add__(self,other):
		if not isinstance(other, PolynomialArray):
			return NotImplemented
		if self.shape!=other.shape:
			raise ValueError('Cannot add. Dimensions differ')
		pa = PolynomialArray(self.shape)
		def _func():
			for e1,e2 in _iterzip(self.entries(),other.entries()):
				yield e1+e2

		self.populate(pa,_func())
		return pa

	def hadamard(self,other):
		'''
			Componentwise multiplication
		'''
		if not isinstance(other, PolynomialArray):
			return NotImplemented
		if self.shape!=other.shape:
			raise ValueError('Cannot multiply componentwise. Dimensions differ')
		pa = PolynomialArray(self.shape)
		def _func():
			for e1,e2 in _iterzip(self.entries(),other.entries()):
				yield e1*e2

		self.populate(pa,_func())
		return pa

	def transform(self,func, in_place=False):
		'''
			replace each entry self[inds] with func(self[inds]),
				which must be a number or a polynomial object
		'''

		if not in_place:
			self = self.reshape(self.shape)

		self.populate(self,_itermap(func, self.entries()))

		if not in_place:
			return self

	def __neg__(self):
		return self.transform(lambda x: -x)

	def __sub__(self,other):
		if not isinstance(other, PolynomialArray):
			return NotImplemented
		return self+(-other)

	def __nonzero__(self):
		return 0 not in self.shape

	def __mul__(self,other):
		try:
			other = ply.to_poly(other)
		except:
			return NotImplemented
		return self.transform(lambda x: x*other)

	def __rmul__(self,other):
		return self*other

	

class PolynomialVector(PolynomialArray):


	def __new__(cls, seq=None, n=None, default = zero_poly):

		if seq is not None:
			seq = list(seq)
			n = len(seq)

		obj = super(PolynomialVector,cls).__new__(cls, n,default)


		if seq is not None:
			obj.populate(obj, _itermap(ply.to_poly,seq))

		return obj

	@staticmethod
	def from_numpy(arr):
		return PolynomialVector(PolynomialArray.from_numpy(arr))


	def __add__(self,other):
		sm = super(PolynomialVector,self).__add__(other)
		return PolynomialVector(sm)

	def __radd__(self,other):
		return self+other

	def transform(self,func,in_place=False):
		return PolynomialVector(super(PolynomialVector,self).transform(func,in_place))


	def __mul__(self,other):
		prod = super(PolynomialVector,self).__mul__(other)
		return PolynomialVector(prod)

	def hadamard(self,other):
		return PolynomialVector(super(PolynomialVector,self).hadamard(other))

	def dot(self,other):
		return sum(self.hadamard(other))

	def as_row(self):
		return PolynomialMatrix([self])

	def as_column(self):
		return PolynomialMatrix([[v] for v in self])

	def gradient(self):
		return PolynomialMatrix([f.gradient() for f in self])

	def __nonzero__(self):
		for entry in self:
			if entry:
				return True
		return False




class PolynomialMatrix(PolynomialArray):

	def __new__(cls, seq = None, nrows=None, ncols=None, default=zero_poly):
		if seq is not None:
			seq = map(PolynomialVector, seq)
			nrows = len(seq)
			ncols = len(seq[0]) if seq else 0
			if len(set(map(len,seq)))!= 1:
				raise ValueError('The row vectors are not all the same size')

		obj = super(PolynomialMatrix, cls).__new__(cls, (nrows,ncols),default)

		if seq is not None:
			for i in xrange(nrows):
				obj.array[i] = PolynomialVector(seq[i])

		return obj

	def square(self):
		return len(set(self.shape))==1

	def __add__(self,other):
		sm = super(PolynomialMatrix,self).__add__(other)
		return PolynomialMatrix(sm)

	def __radd__(self,other):
		return self+other

	@property
	def nrows(self):
		return self.shape[0]

	@property
	def ncols(self):
		return self.shape[1]

	def transpose(self):
		mat = PolynomialMatrix(nrows=self.ncols,ncols=self.nrows)
		for i in xrange(self.nrows):
			for j in xrange(self.ncols):
				mat[j][i] = self[i][j]
		return mat

	def is_row_vector(self):
		return self.nrows==1

	def is_column_vector(self):
		return self.ncols==1

	def symmetric(self):
		if not self.square():
			return False

		n = self.nrows

		for i in xrange(n-1):
			for j in xrange(i+1,n):
				if self[i,j]!=self[j,i]:
					return False
		return True

	def skew_symmetric(self):
		if not self.square():
			return False
		n = self.nrows
		for i in xrange(n-1):
			if self[i,i]:
				return False
			for j in xrange(i+1,n):
				if self[i,j]!=-self[j,i]:
					return False

		return True

	def transform(self, func,in_place):
		return PolynomialMatrix(
			super(PolynomialMatrix,self).transform(func,in_place)
			)

	def is_diagonal(self):
		for i in xrange(n-1):
			for j in xrange(i+1,n):
				if self[i,j] or self[j,i]:
					return False
		return True


	def __nonzero__(self):
		for row in self:
			if row:
				return True
		return False

	def __mul__(self,other):
		if not isinstance(other,PolynomialArray):
			return PolynomialMatrix(super(PolynomialMatrix, self).__mul__(other))

		try:
			other = PolynomialMatrix(other)
		except:
			return NotImplemented

		m,n = self.shape
		nn,p = other.shape
		if n!=nn:
			raise ValueError('Cannot multiply matrices of dimension {}x{} with {}x{}'.format(m,n,nn,p))

		res = PolynomialMatrix(nrows = m, ncols = p)

		for i in xrange(m):
			for j in xrange(n):
				for k in xrange(p):
					res[i][k] += self[i][j]*other[j][k]

		return res


	def __rmul__(self,other):
		if isinstance(other, PolynomialArray):
			try:
				other = PolynomialMatrix(other)
			except:
				return NotImplemented
			return other*self
		return self*other

	@staticmethod
	def from_numpy(arr):
		return PolynomialMatrix(PolynomialArray.from_numpy(arr))




























