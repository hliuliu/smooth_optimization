

import numpy as np 
from numpy.linalg import cholesky,LinAlgError
from numpy.random import random


def is_square(A):
	m,n = A.shape
	return m==n

def is_symmetric(A):
	return is_square(A) and np.all(A==A.T)


def is_SPD(A, assume_symmetric=False):
	if not assume_symmetric and not is_symmetric(A):
		return False

	try:
		cholesky(A)
	except LinAlgError:
		return False

	return True


def random_symmetric_matrix(n):
	m = n*(n+1)//2
	arr = random(m)
	A = np.zeros((n,n))

	k =0
	for i in xrange(n):
		for j in xrange(n):
			if i>j:
				A[i,j]=A[j,i]
			else:
				A[i,j]= arr[k]
				k+=1

	return A



def to_vector(x):
	if len(x.shape)==1:
		return x
	m,n = x.shape
	return x.reshape(m*n)







