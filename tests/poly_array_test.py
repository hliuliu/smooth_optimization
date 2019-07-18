
import os,sys

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)


import poly_array as parr
import polynomial as ply
import numpy as np




######################

#print parr.ply.Polynomial==ply.Polynomial

######################







pv1 = parr.PolynomialArray(7)

print 'pv1=',pv1

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

print 'pv2=',pv2
type(pv2).AUTO_NUMPY = False

pv2_0 = pv2(x=0,y=0)

print type(pv2_0)
assert(type(pv2_0)==parr.PolynomialArray)
assert((pv2_0.to_numpy_array()==np.zeros((4,6))).all())

print 'Trying to convert nonconstant Polynomial array to numpy array... Should raise an error'

try:
	pv2.to_numpy_array()
except:
	print 'PASS'
else:
	print 'FAIL'

pv3 = pv2.reshape(24)
pv4 = pv2.reshape((3,8))
print 'pv2 reshaped to 24: pv3=',pv3
print 'pv2 reshaped to (3,8): pv4=',pv4



for i in xrange(len(pv3)):
	item = pv3[i]
	assert(type(item)==ply.Polynomial)
	assert(item==x**3+y)
	assert(item == pv3.array[i].array)

for i, vec in enumerate(pv4):
	assert(vec is pv4[i])
	assert(vec.dimensions()==pv4.dimensions()[1:])
	for j in xrange(len(vec)):
		p = vec[j]
		assert(type(p)==ply.Polynomial)
		assert(p==pv4[i,j])

pv5 = pv4[:3:2, 4:9]

assert(type(pv5[0].array[0]) is parr.PolynomialArray)

assert(pv5.dimensions()==(2,4))

print pv5.to_nested_list()

ct=0
for entry in pv5.entries():
	assert(type(entry)==ply.Polynomial)
	assert(entry==x**3+y)
	ct+=1

assert(ct==2*4)
