
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


PolyVec = parr.PolynomialVector

ply.Polynomial.AUTO_SORT_VARIABLES = True

parr.PolynomialArray.AUTO_NUMPY = True

x,y,z = ply.init_poly_vars('xyz')

v = range(0,10,3)
pv = parr.PolynomialVector(v)
print pv



assert(type(pv)==parr.PolynomialVector)

for et in pv:
	assert(type(et)==ply.Polynomial)

p1 = x**2+y**2
p2 = (x-1)**3-(y+4)**2
p3 = (x+y)**2-(3*y+2)**3
p4 = x*y**4

print 'p1=',p1
print 'p2=',p2

for p in [p1,p2,p3]:
	assert(p.varlist()==['x','y'])

p1g,p2g,p3g,p4g = map(lambda f:f.gradient(), [p1,p2,p3,p4])

print 'p1g= D(p1)=',p1g

for gr in [p1g,p2g,p3g,p4g]:
	assert(type(gr)==PolyVec)

assert(p2g == PolyVec([3*(x-1)**2, -2*(y+4)]))

assert(p3g == PolyVec([2*(x+y), 2*(x+y)-9*(3*y+2)**2]))

print p4g

assert(list(p1g)==[p1g[0],p1g[1]])

p12 = p1+p2

print p12
assert(type(p1g+p2g)==PolyVec)
assert(p1g+p2g==p12.gradient())


pa0s = parr.PolynomialArray(2)
pa1s = parr.PolynomialArray(2,1)

assert(pa0s+p1g==p1g)

assert(type(pa0s+p1g)==type(p1g))

assert(1*p1g==p1g)
assert(type(1*p1g)==type(p1g))

vec = p1g(2,3)
print vec,type(vec)

print sum([p1,p2,p3])
assert type(sum([p1g,p3g],pa0s))==type(p1g)

print sum(p1g)

print p1g,p4g,p1g.dot(p4g)

print '-',p4g,'=',-p4g

assert type(-p1g)==type(p1g)
assert type(p1g-p2g)==type(p2g)



