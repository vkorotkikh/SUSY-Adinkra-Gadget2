# ******************************************************************************
# Name:    Adinkra NxN Constructor
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:	   February 2017
#
# Description: Algorithm for generating n-colour sized Adinkras of n-node
# matrix representations
#
# ******************************************************************************

# Library Imports
import time, itertools
# from array import array		# going to try and use this
from numpy.linalg import inv
import numpy as np

#>******************************************************************************
def unit_vector(n,i):
	""" Makes np.array rows that compose an identity matrix """
	vec		= [0 for j in range(n)]
	vec[i] 	= 1
	return np.array(vec)

#>******************************************************************************
def gen_permutations(n):
	"""
	Creates the 24 unsigned permutation matrices from permutation group S4
	"""
	# bracket_n = list(range(n))
	perms_list = list(itertools.permutations(range(n), n))
	ts = [unit_vector(n,i) for i in range(n)]
	temp = []
	for mi in perms_list:
		columns = [ts[mi[i]] for i in range(n)]
		temp_mat = np.column_stack(columns)
		matint = temp_mat.astype(int)
		temp.append(matint)
	return temp

#>******************************************************************************
def gen_signm(n_dim):
	"""
	Generate all possible sign permutations for an identiy matrix size n
	"""
	items		= [1] * int(n_dim)
	sign_mat 	= []
	for signs in itertools.product([-1,1], repeat=len(items)):
		temp = np.array([a*sign for a,sign in zip(items,signs)],dtype=int)
		# ptemp.append(temp)
		sign_mat.append(np.diag(temp))
	return sign_mat

#>******************************************************************************
def gen_product_matrices(n):
	"""
	Creates the L sign permutation matrices (384 for 4x4 Matrix)
	"""
	sigp_mats = []
	# from adinkra_tetrad_calc import gen_signm
	sign_pmat = gen_signm(n)
	uperm_mat = gen_permutations(n)
	for x in sign_pmat:
		for y in uperm_mat:
			# sigp_mats.append(np.matmul(t1,t2))
			sigp_mats.append(np.dot(x,y))
	return sigp_mats

#>******************************************************************************
def pairing(matli, matlj):
	"""
		Checks if two L sig-perm matrices are 'pairable'/satisfy Garden Algebra
		conditions
	"""
	ri = np.transpose(matli)
	rj = np.transpose(matlj)
	# tmat = np.dot(matli,rj) + np.dot(matlj,ri)
	rtmat = np.dot(ri,matlj) + np.dot(rj,matli)

	# if np.array_equal(ri, inv(matli)) and np.array_equal(rj, inv(matlj)):
	return (np.count_nonzero(rtmat) == 0)

#>******************************************************************************
def adinkra_gen(k, sigprod_mats):
	""" Make lists of number k matrices (adinkras), using sigprod_mats list
		k is the color number of Adinkra (or size), k=4 is a 4 color, 4 L Matrix
		Adinkra.
	"""
	# main_list	= [[] for i in range(len(sigprod_mats))]
	# print(len(main_list))
	xtest_pack	= []
	fourpack	= []

	if k == 1:
		return [list(l) for l in sigprod_mats]
	else:
		adinkra_list = []
		print("Length lmats", len(sigprod_mats))
		for i, mat in enumerate(sigprod_mats):			# Find all matrix pairs for mat
			# good_mats = [m for m in sigprod_mats if pairing(mat,m)]
			# test_mats = [ind[0] for ind in enumerate(sigprod_mats) if pairing(mat, ind[1])]
			xtest_pack = [ind for ind in enumerate(sigprod_mats) if pairing(mat, ind[1])]
			# main_list[i] = xtest_pack
			for val in xtest_pack:
				# main_list[val[0]] = [(i,mat)]
				fourpack	= [nmat for nmat in xtest_pack if pairing(val[1], nmat[1])]
				# print(len(fourpack))
				for xval in fourpack:
					for lastx in [nmat for nmat in fourpack if pairing(xval[1], nmat[1])]:
						# temp = [(i,mat), val, xval, lastx]
						# temp = [i, val[0], xval[0], lastx[0]]
						adinkra_list.append([(i,mat), val, xval, lastx])
		return adinkra_list


#>******************************************************************************
def create_adinkras(k,n):
	# main_tetrad = adinkra_gen(k,gen_product_matrices(n))
	return adinkra_gen(k,gen_product_matrices(n))
