

import os,sys

sys.path.append(
	os.path.join(
			os.path.dirname(__file__),
			os.pardir
		)
	)

import polynomial as ply


p = ply.Polynomial(['x','y'], [(1,[2,0]),(2,[0,2])])
q = ply.Polynomial(['x','y','z'], [(1,[2,1,0]),(2,[0,2,1])])

print 'p =',p
print 'q =',q

print ply.constant_poly(-1)

c = ply.Polynomial(['x','y','z'], [(5,[0,0,0])])

print 'c =',c

print 'p constant?',p.is_constant()
print 'c constant?',c.is_constant()


print 'coef of x^2 in p =', p.coef_of(dict(x=2))
print 'coef in z in q =', q.coef_of(dict(z=1))
print 'partial coef of z in q =', q.coef_of(dict(z=1), True)
print 'coef in y^2*z^0 in p =', p.coef_of(dict(y=2,z=0))

print 'coef in y^2*z^1 in p =', p.coef_of(dict(y=2,z=1))


x,y,z = ply.init_poly_vars('xyz')

print zip('xyz',map(str,[x,y,z]))


print 'p(1,1)=',p(1,1)
print 'q(-2,2,1)=',q(-2,2,1)

print 'p+q=',p+q

r = ply.Polynomial(['x','y','z'], [(4,[2,0,0]),(7,[2,3,1])])

print 'r=',r
print 'p+r=',p+r

print '-'*20

test_exprs =['p+q==p.get_sum(p,q)', 'p+q==q+p', '(p+q)+r==p+(q+r)', 'p+q+r==p.get_sum(p,q,r)', 
	'c==ply.constant_poly(5)', 'c==5', '5==c']

for comp in test_exprs:
	print comp,eval(comp)


print 'p+5 =', p+5
print '5+p =', 5+p

for func in 'pqrc':
	print '{}.degree() ='.format(func), eval(func).degree()

print 'p*q =', p*q
print 'p*r =', p*r
print 'q*c =', q*c


print 'x**2+2*y*x**3 =',x**2+2*y*x**3

f = x**3+2*y**4-z
g = x**2*y**5+z

print 'f(x,y,z) =',f
print 'g(x,y,z) =',g 

h = f*g

print 'h=f*g=',h

for ex in 'fgh':
	for wrt in 'wxyz':
		print 'd/d%s [%s(x,y,z)] ='%(wrt,ex), eval(ex).derivative(wrt)


p1 = (x+y)*(x-y)
p2 = x**2-y**2

assert(p1==p2)

p1 = (x-y)*(x**2+x*y+y**2)
p2 = x**3-y**3

assert(p1==p2)


def choose(n,k):
	if k>n//2:
		k=n-k
	ans = 1
	for i in xrange(k):
		ans *= n-i
		ans /= i+1
	return ans

def bin_thm_test(x,y, n):
	print 'expanding ({}+{})^{}'.format(x,y,n)
	t = (x+y)**n
	print 'result:',t
	for i in xrange(n+1):
		t -= choose(n,i)*x**i*y**(n-i)

	print 'PASS' if not t else 'FAIL'

	print 


for n in [0,1,2,3,5,7,9,12,20]:
	bin_thm_test(x,y,n)



test_exprs = ['h==g*f', 'p.derivative("x")+q.derivative("x")==(p+q).derivative("x")',
	'(p*q).derivative("y")==p.derivative("y")*q+p*q.derivative("y")']


for comp in test_exprs:
	print comp,eval(comp)

print 'f->', f.varlist()
print 'g->', g.varlist()
print 'h->', h.varlist()
fpg = f+g

print 'f+g=',fpg,'->', (fpg).varlist()
print 'sorting variables'
fpgs = ply.sorted_variables(fpg)

print 'f+g=',fpgs,'->', (fpgs).varlist()


fpg_xy = fpgs(x,y,0)
print '(f+g)(x,y)', fpg_xy, fpg_xy.varlist()

print 'sorting variables'
fpg_xys = ply.sorted_variables(fpg_xy)

print '(f+g)(x,y)', fpg_xys, fpg_xys.varlist()

def all_equal(seq):
	for i in xrange(len(seq)-1):
		if seq[i]!=seq[i+1]:
			print seq[i], seq[i+1], 'not equal'
			return False
	return True

assert(all_equal([fpg,fpgs,fpg_xy,fpg_xys]))

