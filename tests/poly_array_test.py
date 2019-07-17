
import os,sys

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)


import poly_array as parr
import polynomial as ply




######################

#print parr.ply.Polynomial==ply.Polynomial

######################







pv1 = parr.PolynomialArray(7)

assert(pv1.dimensions()==(7,))

for p in pv1.array:
	assert(p.dimensions()==())
	assert(isinstance(p,parr.PolynomialArray))
	assert(p.array==0)


pv1_vec = pv1.to_numpy_array()
print pv1_vec
print type(pv1_vec)
print pv1()

x,y = ply.init_poly_vars('xy')

print type(x**3+y)
print isinstance(x**3+y, parr.Poly)
print parr.Poly, ply.Polynomial
print issubclass(type(x**3+y),parr.Poly)

pv2 = parr.PolynomialArray((4,6),x**3+y
	)

for p in pv2.array:
	assert(p.dimensions()==(6,))
	for q in p.array:
		assert(q.dimensions()==())
		assert(q.array==x**3+y)

